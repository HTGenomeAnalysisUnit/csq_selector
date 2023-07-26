import hts
import times
import strformat
import strutils
import tables
import os
import streams
import re
import json
import std/sets
import csq_selector/impact_order
import csq_selector/utils
import csq_selector/tx_expression
import csq_selector/split_ann
import csq_selector/arg_parse

const VERSION = "0.1.1"
const TSV_HEADER = "#CHROM\tPOS\tID\tREF\tALT\tFILTER\tGENE_ID\tGENE_SYMBOL\tTRANSCRIPT\tCONSEQUENCE"

proc write_new_var(wrt:VCF, v:Variant): bool {.inline.} =
    result = wrt.write_variant(v)

proc var2string(v:Variant, csq: string): string {.inline.} =
    var alts: seq[string]
    for a in v.ALT: alts.add(a) 

    let outline = @[ $v.CHROM, $v.POS, $v.ID, v.REF, alts.join(","), v.FILTER, csq ]
    result = outline.join("\t")

proc write_new_var(wrt:FileStream, v:Variant, csq: string): bool {.inline.} =
    let outline = var2string(v, csq)
    try:
        wrt.writeLine(outline)
        result = true
    except:
        result = false

# Given a seq of impacts for a variant, generate a seq of strings representing variant anno in regenie format
proc write_anno_string(wrt:FileStream, v: Variant, csqs: seq[string], useid: bool = false): int =
  result = 0
  var var_id = [$v.CHROM, $v.POS, $v.ID, v.REF, v.ALT[0]].join(":")
  if useid: 
    var_id = $v.ID

  for c in csqs:
    try:
        wrt.writeLine([var_id,c].join("\t"))
        result += 1
    except:
        discard

# TODO: Future functions to write set file and annot file as per regenie format
# proc write_new_set(wrt:FileStream, v:Variant, csq: string): bool {.inline.} =
#     discard

iterator readvar(v: VCF, regions: seq[string]): Variant =
    if regions.len == 0:
        for variant in v: yield variant
    else:
        for r in regions:
            for variant in v.query(r): yield variant

proc main* () =
    log("INFO", fmt"CSQ SELECTOR v{VERSION}")
    
    let TX_VERS_RE = re"\.[0-9]+"
    
    #Load default impact order and show it if requested
    var impact_order = default_order
    var cmdline = commandLineParams()
    if cmdline[0] == "show_csq_order":
        log("INFO", "Default impact order:")
        echo show_impact_order(impact_order)
        quit QuitSuccess
    
    #Parse command line arguments
    var opts = parseCmdLine()
    opts.logArgs()

    #Check at least one filter is active for VCF output
    var most_severe = opts.most_severe
    if not (most_severe or opts.transcripts != "" or opts.min_impact != "" or opts.min_exp != ""):
        if opts.out_format == "tsv":
            log("WARN", "No filters active, just output TSV of all consequences")
        else:
            log("FATAL", "No filters active for a VCF output, nothing to do")
            quit "", QuitFailure

    # Check we have an output file if format is rarevar_set
    if opts.out_format == "rarevar_set" and opts.out == "":
        log("FATAL", "Output file must be specified for rarevar_set output")
        quit "", QuitFailure

    # Set consequences order
    if opts.impact_order != "":
        if not fileExists(opts.impact_order):
            log("FATAL", fmt"Impact order file {opts.impact_order} does not exist")
            quit "", QuitFailure
        impact_order = adjustOrder(opts.impact_order.readFile)

    # Load lists
    let
        transcripts_list = read_list(opts.transcripts, "transcript")
        regions = read_list(opts.region, "region")
        tissues = read_list(opts.tissues, "tissue")
    let allowed_transcripts = transcripts_list.cleanTxVersion(TX_VERS_RE).toHashSet

    # Check min impact is valid
    var min_impact: string

    if opts.min_impact != "":
        min_impact = ( if opts.min_impact.endsWith("CUTOFF"): opts.min_impact else: opts.min_impact.toLowerAscii)
        if min_impact.endsWith("_variant"):
            min_impact = min_impact[0..min_impact.high - 8]

        if not impact_order.hasKey(min_impact):
            log("FATAL", "Impact level " & opts.min_impact & " not found in known impacts order")
            quit "", QuitFailure


    # Load expression data if provided
    var 
        ranked_exp: HashSet[string]
        min_exp: float
    if opts.min_exp != "":
        min_exp = parseFloat(opts.min_exp)
        if opts.exp_data == "":
            log("FATAL", "--min_exp requires --exp_data to be set")
            quit "", QuitFailure
        if not fileExists(opts.exp_data):
            log("FATAL", "--exo_data file does not exist")
            quit "", QuitFailure
        if tissues.len == 0:
            log("FATAL", "No tissues can be loaded from file/list specified by --tissues")
            quit "", QuitFailure
        ranked_exp = read_exp(opts.exp_data,0,1,2,tissues, min_exp, TX_VERS_RE)
        if ranked_exp.len == 0:
            log("FATAL", fmt"No transcripts ranked above the minimum expression threshold ({min_exp})")
            quit "", QuitFailure

    # Load scores from JSON if provided
    var scores_json: JsonNode
    var scores_keys: seq[string]
    if opts.scores != "":
        if opts.out_format == "vcf":
            log("WARNING", "--scores can only be used with rarevar_set output and will be ignored")
        else:
            scores_json = parseFile(opts.scores)
            for k in scores_json.fields.keys:
                scores_keys.add(k)
    
    # Process variants
    log("INFO", "Variant processing started")
    var
        start_time = cpuTime()
        t0 = cpuTime()
        written_vars = 0
        interval = 10000
        n = 0
        n_noimpact = 0
        n_nocsqfield = 0
        by_gene = (if opts.mode == "by_gene": true else: false)

    var vcf:VCF
    if not open(vcf, opts.vcf):
        log("FATAL", "Could not open input VCF file " & opts.vcf)
        quit "", QuitFailure

    # Get GeneIndex from header
    var gene_fields: GeneIndexes
    let csq_columns = (if opts.csq_column != "": opts.csq_column.split(",") else: @[])
    vcf.set_csq_fields_idx(opts.csq_field, gene_fields, csq_columns)

    # Set header
    var out_vcf:VCF
    var out_tsv:FileStream
    var out_setlist:FileStream
    var out_annot:FileStream
    var gene_set: Table[string,GeneSet]
    
    case opts.out_format:
        of "vcf":
            var out_header = vcf.header
            if opts.keep_old_ann:
                if out_header.add_info("ORIGINAL_ANN", "1","String", "Original annotation before CSQ_SELECTOR") != Status.OK:
                    log("FATAL", "Could not add ORIGINAL_ANN to VCF header")
                    quit "", QuitFailure
            if opts.out != "":           
                doAssert(open(out_vcf, opts.out, mode="w"))
                out_vcf.header = out_header
                doAssert(out_vcf.write_header())
            else:
                stdout.writeLine($out_header)
        of "tsv":
            var header = TSV_HEADER
            for c in gene_fields.columns.keys():
                header &= &"\t{c}"  
            if opts.out != "":
                out_tsv = newFileStream(opts.out, fmWrite)
                if isNil(out_tsv):
                    log("FATAL", "Could not open output file " & opts.out)
                    quit "", QuitFailure
                else:
                    out_tsv.writeLine(header)
            else:
                stdout.writeLine(header)
        of "rarevar_set":
            out_setlist = newFileStream(fmt"{opts.out}.setlist", fmWrite)
            out_annot = newFileStream(fmt"{opts.out}.annot", fmWrite)

    #Process variants
    for v in vcf.readvar(regions):
        n = n + 1
        var (dolog, log_msg) = progress_counter(n, interval, t0)
        if dolog: log("INFO", log_msg)
        var (csqfield_missing, impacts) = v.split_csqs(opts.csq_field, gene_fields, impact_order, TX_VERS_RE, min_impact, allowed_transcripts, ranked_exp, scores_json)
        n_nocsqfield += csqfield_missing
        
        if impacts.len == 0: 
            n_noimpact += 1
        else:
            if opts.most_severe:
                impacts = impacts.get_most_severe(gene_fields, by_gene)
    
        if opts.keep_old_ann and csqfield_missing == 0 and opts.out_format == "vcf":
            var old_ann: string
            doAssert v.info.get(opts.csq_field, old_ann) == Status.OK
            doAssert v.info.set("ORIGINAL_ANN", old_ann) == Status.OK

        let selected_csqs = impacts.get_csq_string(gene_fields, opts.out_format)
        
        if opts.out_format == "vcf":
            if selected_csqs.len > 0:
                var new_ann = selected_csqs.join(",")
                doAssert v.info.set(opts.csq_field, new_ann) == Status.OK
            else:
                if csqfield_missing == 0:
                    doAssert v.info.delete(opts.csq_field) == Status.OK
            if opts.out != "":
                if out_vcf.write_new_var(v): written_vars += 1
            else:
                written_vars += 1
                stdout.writeLine(v.tostring)
        elif opts.out_format == "tsv":
            for x in selected_csqs:
                if opts.out != "":
                    if out_tsv.write_new_var(v, x): written_vars += 1
                else:
                    written_vars += 1
                    stdout.writeLine(v.var2string(x))
        elif opts.out_format == "rarevar_set":
            written_vars = written_vars + out_annot.write_anno_string(v, selected_csqs, opts.use_vcf_id)
            # update_gene_set*(gene_set: var Table[string, Gene_set], v: Variant, csqs: seq[Impact], useid: bool = false) 
            gene_set.update_gene_set(v, impacts, opts.use_vcf_id)

    close(vcf)
    if opts.out != "":
        case opts.out_format:
            of "vcf": close(out_vcf) 
            of "tsv": close(out_tsv)
            of "rarevar_set":
                for line in gene_set.make_set_string:
                    out_setlist.writeLine(line)
                close(out_setlist)
                close(out_annot)

    if n_noimpact > 0:
        log("WARN", fmt"{n_noimpact} variants had no impact after filters")
    if n_nocsqfield > 0:
        log("WARN", fmt"{n_nocsqfield} variants had no " & opts.csq_field & " field")
    log("INFO", fmt"{written_vars} records written in {elapsed_time(start_time)}")

when isMainModule:
  main()
