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

const VERSION = "0.4"
const TSV_HEADER = "#CHROM\tPOS\tID\tREF\tALT\tFILTER\tGENE_ID\tGENE_SYMBOL\tTRANSCRIPT\tCONSEQUENCE"

proc write_new_var(wrt:VCF, v:Variant): bool {.inline.} =
    result = wrt.write_variant(v)

proc var2string(v:Variant, csq: string, info_fields: seq[string]): string {.inline.} =
    var alts: seq[string]
    for a in v.ALT: alts.add(a) 

    var outline = @[ $v.CHROM, $v.POS, $v.ID, v.REF, alts.join(","), v.FILTER, csq ]
    var info_string: string
    for x in info_fields:
        if v.info.get(x, info_string) == Status.OK:
            outline.add(info_string)
        else:
            outline.add(".")
    result = outline.join("\t")

proc write_new_var(wrt:FileStream, v:Variant, csq: string, info_fields: seq[string]): bool {.inline.} =
    let outline = var2string(v, csq, info_fields)
    try:
        wrt.writeLine(outline)
        result = true
    except:
        result = false

# Given a seq of impacts for a variant, generate a seq of strings representing variant anno in regenie format
proc write_anno_string(wrt:FileStream, v: Variant, csqs: seq[string], useid: bool = false): int =
  result = 0
  var var_id = [$v.CHROM, $v.POS, v.REF, v.ALT[0]].join(":")
  if useid: 
    var_id = $v.ID

  for c in csqs:
    try:
        wrt.writeLine([var_id,c].join("\t"))
        result += 1
    except:
        log("WARN", fmt"Could not write annotation string for variant {var_id} and consequence {c}")

# Convenience iterator to read variants from a VCF with or without regions
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

    let 
        by_gene = (if opts.mode == "by_gene": true else: false)
        csq_columns = (if opts.csq_column != "": opts.csq_column.split(",") else: @[])
        info_fields = (if opts.info_column != "": opts.info_column.split(",") else: @[])
        most_severe = opts.most_severe

    #Check at least one filter is active for VCF output
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
        if tissues.len == 0:
            log("FATAL", "No tissues can be loaded from file/list specified by --tissues")
            quit "", QuitFailure
        log("INFO", fmt"Loading expression data from {opts.exp_data}")
        if not fileExists(opts.exp_data):
            log("FATAL", "--exp_data file does not exist")
            quit "", QuitFailure

        ranked_exp = read_exp(opts.exp_data,0,1,2,tissues, min_exp, TX_VERS_RE)
        if ranked_exp.len == 0:
            log("FATAL", fmt"No transcripts ranked above the minimum expression threshold ({min_exp})")
            quit "", QuitFailure

    # Load tag shema from JSON if provided
    var tag_config_json: JsonNode
    var tag_csq_keys: seq[string]
    var tag_info_fields: seq[string]
    var tag_csq_fields: seq[string]
    if opts.var_tagging_json != "":
        if opts.out_format == "vcf":
            log("WARNING", "--scores can only be used with rarevar_set output and will be ignored")
        else:
            log("INFO", fmt"Loading tagging schema from {opts.var_tagging_json}")
            tag_config_json = parseFile(opts.var_tagging_json)
            for k in tag_config_json["csq_classes"].keys:
                echo k
                tag_csq_keys.add(k)
                let tagging_config = (if tag_config_json["csq_classes"][k].hasKey("tagging"): tag_config_json["csq_classes"][k]["tagging"] else: %* {})
                let scoring_config = (if tag_config_json["csq_classes"][k].hasKey("scoring"): tag_config_json["csq_classes"][k]["scoring"] else: %* {})
                if scoring_config.hasKey("info"):
                    for s in scoring_config["info"].keys:
                        tag_info_fields.add(s)
                if scoring_config.hasKey("csq_field"):
                    for s in scoring_config["csq_field"].keys:
                        tag_csq_fields.add(s)
                for t in tagging_config.items:
                    if t.hasKey("info"):
                        for s in t["info"].keys:
                            tag_info_fields.add(s)
                    if t.hasKey("csq_field"):
                        for s in t["csq_field"].keys:
                            tag_csq_fields.add(s)
            log("INFO", fmt"Loaded scores for {$tag_csq_keys} consequences")
            log("INFO", fmt"Looking for the following fields from INFO {$tag_info_fields}")
            log("INFO", fmt"Looking for the following fields from CSQ field {$tag_csq_fields}")
    
    # Set CSQ config
    var csq_config: Config
    csq_config.csq_field_name = opts.csq_field
    csq_config.min_impact = opts.min_impact
    csq_config.allowed_tx = allowed_transcripts
    csq_config.ranked_exp = ranked_exp
    csq_config.most_severe = most_severe
    #csq_config.most_expressed: bool # not implemented yet
    csq_config.group_by_gene = by_gene
    csq_config.tagging_config = tag_config_json
    csq_config.csq_output_fields = csq_columns
    csq_config.tx_vers_re = TX_VERS_RE

    # Set output streams
    var out_vcf:VCF
    var out_tsv:FileStream
    var out_setlist:FileStream
    var out_annot:FileStream
    var gene_set: Table[string,GeneSet]

    # Set vars and open input VCF
    log("INFO", "Variant processing started")
    var
        start_time = cpuTime()
        t0 = cpuTime()
        written_vars = 0
        interval = 10000
        n = 0
        n_noimpact = 0
        n_nocsqfield = 0
        unknown_impacts: HashSet[string]

    var vcf:VCF
    if not open(vcf, opts.vcf):
        log("FATAL", "Could not open input VCF file " & opts.vcf)
        quit "", QuitFailure

    # Get GeneIndex from header
    let (csq_field, csq_fields_idx) = vcf.set_csq_fields_idx(opts.csq_field, csq_columns)
    csq_config.csq_field_idxs = csq_fields_idx
    csq_config.csq_field_name = csq_field
    
    # If we have scores are set, check they are defined in the header
    var desc: string
    for s in tag_info_fields:
        try:
            desc = vcf.header.get(s, BCF_HEADER_TYPE.BCF_HL_INFO)["Description"]
        except:
            log("FATAL", fmt"Could not find {s} in VCF header")
            quit "", QuitFailure

    # Open needed output streams and set headers
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
            for c in csq_config.csq_output_fields:
                header &= &"\tCSQ_{c}"
            for c in info_fields:
                header &= &"\tINFO_{c}"
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

    # Process variants
    for v in vcf.readvar(regions):
        n = n + 1
        var (dolog, log_msg) = progress_counter(n, interval, t0)
        if dolog: log("INFO", log_msg)
        
        var (csqfield_missing, impacts) = v.split_csqs(csq_config, impact_order, unknown_impacts)

        if opts.keep_old_ann and opts.out_format == "vcf" and not csqfield_missing:
            var old_ann: string
            doAssert v.info.get(csq_config.csq_field_name, old_ann) == Status.OK
            doAssert v.info.set("ORIGINAL_ANN", old_ann) == Status.OK

        n_nocsqfield += (if csqfield_missing: 1 else: 0)
        
        if impacts.len == 0: 
            n_noimpact += 1
    
        let selected_csqs = impacts.get_csq_string(csq_config.csq_output_fields, opts.out_format)
        
        if opts.out_format == "vcf":
            if selected_csqs.len > 0:
                var new_ann = selected_csqs.join(",")
                doAssert v.info.set(csq_config.csq_field_name, new_ann) == Status.OK
            else:
                if not csqfield_missing:
                    doAssert v.info.delete(csq_config.csq_field_name) == Status.OK
            if opts.out != "":
                if out_vcf.write_new_var(v): written_vars += 1
            else:
                written_vars += 1
                stdout.writeLine(v.tostring)
        elif opts.out_format == "tsv":
            for x in selected_csqs:
                if opts.out != "":
                    if out_tsv.write_new_var(v, x, info_fields): written_vars += 1
                else:
                    written_vars += 1
                    stdout.writeLine(v.var2string(x, info_fields))
        elif opts.out_format == "rarevar_set":
            written_vars = written_vars + out_annot.write_anno_string(v, selected_csqs, opts.use_vcf_id)
            gene_set.update_gene_set(v, impacts, opts.use_vcf_id)

    log("INFO", fmt"Finished processing {n} variants")
    close(vcf)

    # Close output streams and write setlist file if needed
    if opts.out != "":
        case opts.out_format:
            of "vcf": close(out_vcf) 
            of "tsv": close(out_tsv)
            of "rarevar_set":
                log("INFO", fmt"Writing {gene_set.len} gene sets to {opts.out}.setlist")
                for set_string in gene_set.make_set_string:
                    out_setlist.writeLine(set_string)
                close(out_setlist)
                close(out_annot)

    # Final log
    if n_noimpact > 0:
        log("WARN", fmt"{n_noimpact} variants had no impact after filters")
    if n_nocsqfield > 0:
        log("WARN", fmt"{n_nocsqfield} variants had no " & opts.csq_field & " field")
    log("INFO", fmt"{written_vars} records written in {elapsed_time(start_time)}")

when isMainModule:
  main()
