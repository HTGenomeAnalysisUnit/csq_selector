# This is just an example to get you started. A typical binary package
# uses this file as the main entry point of the application.

import hts
import times
import strformat
import strutils
import tables
import os
import streams
import re
import std/sets
import csq_selector/impact_order
import csq_selector/utils
import csq_selector/tx_expression
import csq_selector/split_ann
import csq_selector/arg_parse

const VERSION = "0.1.0"
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

iterator readvar(v: VCF, regions: seq[string]): Variant =
    if regions.len == 0:
        for variant in v: yield variant
    else:
        for r in regions:
            for variant in v.query(r): yield variant

proc main* () =
    let TX_VERS_RE = re"\.[0-9]+"
    
    var p: ParseSchema
    p.initParser("CSQ Selectot: select gene consequences from a VCF file"):
        p.addOption("-i", "--vcf", help="path to VCF/BCF", required=true)
        p.addOption("-o", "--out", help="Output file")
        p.addOption("-O", "--out_format", help="Output format", default = "vcf", choices = @["vcf", "tsv"])
        p.addOption("-c", "--csq_field", help="INFO field containing the gene name and impact. Usually CSQ, ANN or BCSQ", default = "ANN")
        p.addOption(long="--csq-column", help="CSQ sub-field(s) to extract (in addition to gene, impact, transcript) when output is TSV. A comma-separated list of sub-fields can be provided")
        p.addOption("-e", "--exp_data", help="path to expression file. Tab-separated table of transcript expression across conditions")
        p.addOption("-n", "--min_exp", help="Min expression value for a transcript to be considered expressed")
        p.addOption(long="--tissues", help="List of tissues to select in the expression file. Comma-separated list or file with 1 tissue per line.")        
        p.addOption("-r", "--region", help="Specify genomic region for processing. Format is chr[:start-end]. Comma-separated list of regions or file with 1 region per line.")                
        p.addOption("-t", "--transcripts", help="list of trancsripts to restrict consequences output. Comma-separated list or file with 1 tissue per line.")
        p.addOption("-m", "--mode", help="Mode to apply selected filters. Across all annotations or by gene", default="all", choices = @["all", "by_gene"])        
        p.addFlag(long="--most_severe", help="Only output the most severe consequence for each variant")
        p.addFlag(long="--keep_old_ann", help="Keep original annotation fields in the output VCF under ORIGINAL_ANN tag")
        p.addOption(long="--min_impact", help="Impact threashold from the consequences order")
        p.addOption(long="--impact_order", help="ordering of impacts to override the default. See default using --show-csq-order")
        p.addFlag(long="--show_csq_order", help="Show the default CSQ order and exit")

    var opts = p.parseOptions()
    #if len(argv) == 0: argv = @["--help"]
    log("INFO", fmt"Starting csq_selector {VERSION}")

    # Log parsed arguments
    for k, v in opts.pairs():
        log("ARG", fmt"{k} = {v}")

    # Check if some filtering is active
    var most_severe: bool 
    opts.getOpt("most_severe", most_severe)

    if not (most_severe or opts.isSet("transcripts") or opts.isSet("min_impact") or opts.isSet("min_exp")):
        if opts["out_format"] == "tsv":
            log("WARN", "No filters active, just output TSV of all consequences")
        else:
            log("FATAL", "No filters active for a VCF output, nothing to do")
            quit "", QuitFailure

    # Set consequences order
    var impact_order = default_order
    if opts.isSet("impact_order"):
        if not fileExists(opts["impact_order"]):
            log("FATAL", "Impact order file " & opts["impact_order"] & " does not exist")
            quit "", QuitFailure
        impact_order = adjustOrder(opts["impact_order"].readFile)

    # Load lists
    let
        transcripts_list = read_list(opts["transcripts"], "transcript")
        regions = read_list(opts["region"], "region")
        tissues = read_list(opts["tissues"], "tissue")
    let allowed_transcripts = transcripts_list.cleanTxVersion(TX_VERS_RE).toHashSet

    # Check min impact is valid
    var min_impact: string
    if opts.isSet("min_impact"):
        if opts["min_impact"] == "":
            if impact_order.hasKey("IMPACTFUL_CUTOFF"): min_impact = "IMPACTFUL_CUTOFF"
        else:
            min_impact = ( if opts["min_impact"].endsWith("CUTOFF"): opts["min_impact"] else: opts["min_impact"].toLowerAscii)
            if min_impact.endsWith("_variant"):
                min_impact = min_impact[0..min_impact.high - 8]

        if not impact_order.hasKey(min_impact):
            log("FATAL", "Impact level " & opts["min_impact"] & " not found in known impacts order")
            quit "", QuitFailure
    else:
        min_impact = ""

    # Load expression data if provided
    var 
        ranked_exp: HashSet[string]
        min_exp: float
    if opts.isSet("min_exp"):
        opts.getOpt("min_exp", min_exp)
        if not opts.isSet("exp_data"):
            log("FATAL", "--min_exp requires --exp_data to be set")
            quit "", QuitFailure
        if not fileExists(opts["exp_data"]):
            log("FATAL", "--exo_data file does not exist")
            quit "", QuitFailure
        if tissues.len == 0:
            log("FATAL", "No tissues can be loaded from file/list specified by --tissues")
            quit "", QuitFailure
        ranked_exp = read_exp(opts["exp_data"],0,1,2,tissues, min_exp, TX_VERS_RE)
        if ranked_exp.len == 0:
            log("FATAL", fmt"No transcripts ranked above the minimum expression threshold ({min_exp})")
            quit "", QuitFailure

    # Process variants
    log("INFO", "Starting variant processing")
    var
        start_time = cpuTime()
        t0 = cpuTime()
        written_vars = 0
        interval = 5000
        n = 0
        n_noimpact = 0
        n_nocsqfield = 0
        by_gene = (if opts["mode"] == "by_gene": true else: false)

    var vcf:VCF
    if not open(vcf, opts["vcf"]):
        log("FATAL", "Could not open input VCF file " & opts["vcf"])
        quit "", QuitFailure

    # Get GeneIndex from header
    var gene_fields: GeneIndexes
    let csq_columns = (if opts.isSet("csq_column"): opts["csq_column"].split(",") else: @[])
    vcf.set_csq_fields_idx(opts["csq_field"], gene_fields, csq_columns)

    # Set header
    var out_vcf:VCF
    var out_tsv:FileStream
    
    if opts["out_format"] == "vcf":
        var out_header = vcf.header
        if opts.isSet("keep_old_ann"):
            if out_header.add_info("ORIGINAL_ANN", "1","String", "Original annotation before CSQ_SELECTOR") != Status.OK:
                log("FATAL", "Could not add ORIGINAL_ANN to VCF header")
                quit "", QuitFailure
        if opts.isSet("out"):           
            doAssert(open(out_vcf, opts["out"], mode="w"))
            out_vcf.header = out_header
            doAssert(out_vcf.write_header())
        else:
            stdout.writeLine($out_header)
    else:
        var header = TSV_HEADER
        for c in gene_fields.columns.keys():
            header &= &"\t{c}"  
        if opts.isSet("out"):
            out_tsv = newFileStream(opts["out"], fmWrite)
            if isNil(out_tsv):
                log("FATAL", "Could not open output file " & opts["out"])
                quit "", QuitFailure
            else:
                out_tsv.writeLine(header)
        else:
            stdout.writeLine(header)

    #Process variants
    for v in vcf.readvar(regions):
        n = n + 1
        var (dolog, log_msg) = progress_counter(n, interval, t0)
        if dolog: log("INFO", log_msg)
        var (csqfield_missing, impacts) = v.split_csqs(opts["csq_field"], gene_fields, impact_order, TX_VERS_RE, min_impact, allowed_transcripts, ranked_exp)
        n_nocsqfield += csqfield_missing
        
        if impacts.len == 0: 
            n_noimpact += 1
        else:
            if opts.isSet("most_severe"):
                impacts = impacts.get_most_severe(gene_fields, by_gene)
    
        if opts.isSet("keep_old_ann") and csqfield_missing == 0 and opts["out_format"] == "vcf":
            var old_ann: string
            doAssert v.info.get(opts["csq_field"], old_ann) == Status.OK
            doAssert v.info.set("ORIGINAL_ANN", old_ann) == Status.OK

        let selected_csqs = impacts.get_csq_string(gene_fields, opts["out_format"])
        
        if opts["out_format"] == "vcf":
            if selected_csqs.len > 0:
                var new_ann = selected_csqs.join(",")
                doAssert v.info.set(opts["csq_field"], new_ann) == Status.OK
            else:
                if csqfield_missing == 0:
                    doAssert v.info.delete(opts["csq_field"]) == Status.OK
            if opts.isSet("out"):
                if out_vcf.write_new_var(v): written_vars += 1
            else:
                written_vars += 1
                stdout.writeLine(v.tostring)
        else:
            for x in selected_csqs:
                if opts.isSet("out"):
                    if out_tsv.write_new_var(v, x): written_vars += 1
                else:
                    written_vars += 1
                    stdout.writeLine(v.var2string(x))

    close(vcf)
    if opts.isSet("out"):
        case opts["out_format"]:
            of "vcf": close(out_vcf) 
            of "tsv": close(out_tsv)

    if n_noimpact > 0:
        log("WARN", fmt"{n_noimpact} variants had no impact after filters")
    if n_nocsqfield > 0:
        log("WARN", fmt"{n_nocsqfield} variants had no " & opts["csq_field"] & " field")
    log("INFO", fmt"{written_vars} variants written in {elapsed_time(start_time)}")

when isMainModule:
  main()
