import argparse
import strformat
from ./utils import log

var p = newParser("CSQ Selector"):
    help("Select gene consequences from snpEff, VEP or bcftools annotations in VCF file")
    option("-i", "--vcf", help="path to VCF/BCF", required=true)
    option("-o", "--out", help="Output file")
    option("-O", "--out_format", help="Output format", default = some("vcf"), choices = @["vcf", "tsv", "rarevar_set"])
    option("-c", "--csq_field", help="INFO field containing the gene name and impact. Usually CSQ, ANN or BCSQ", default = some("ANN"))
    option("--csq_column", help="CSQ sub-field(s) to extract (in addition to gene, impact, transcript) when output is TSV. A comma-separated list of sub-fields can be provided")
    option("--info_column", help="INFO field(s) to include in the TSV output. A comma-separated list can be provided")    
    option("-e", "--exp_data", help="path to expression file. Tab-separated table of transcript expression across conditions")
    option("-n", "--min_exp", help="Min expression value for a transcript to be considered expressed")
    option("--tissues", help="List of tissues to select in the expression file. Comma-separated list or file with 1 tissue per line.")        
    option("-r", "--region", help="Specify genomic region for processing. Format is chr[:start-end]. Comma-separated list of regions or file with 1 region per line.")                
    option("-t", "--transcripts", help="list of trancsripts to restrict consequences output. Comma-separated list or file with 1 tissue per line.")
    option("--mode", help="Mode to apply most severe filter. Across all annotations or by gene", default=some("all"), choices = @["all", "by_gene"])        
    flag("-m", "--most_severe", help="Only output the most severe consequence for each variant")
    flag("--keep_old_ann", help="Keep original annotation fields in the output VCF under ORIGINAL_ANN tag")
    option("--min_impact", help="Impact threshold from the consequences order")
    option("--impact_order", help="ordering of impacts to override the default. See default using show_csq_order")
    option("-j", "--var_tagging_json", help="A JSON file describing schema for variant tagging")
    flag("--use_vcf_id", help="Use the ID field from the VCF as the variant ID in the rarevar_set output")
        
proc parseCmdLine*(): ref =
    try:
        result = p.parse() 
    except ShortCircuit as e:
        if e.flag == "argparse_help":
            echo p.help
            quit QuitSuccess
    except UsageError:
        stderr.writeLine getCurrentExceptionMsg() 
        echo "Use --help for usage information"
        quit QuitSuccess

proc logArgs*(opts: ref) {.discardable.} =
    var active_filters: seq[string]
    log("ARG", fmt"Input VCF: {opts.vcf}")
    let out_stream = (if opts.out != "": opts.out else: "stdout")
    log("ARG", fmt"Output: {out_stream}")
    log("ARG", fmt"Output format: {opts.out_format}")
    log("ARG", fmt"CSQ field: {opts.csq_field}")
    log("ARG", fmt"Keep original CSQ field: {opts.keep_old_ann}")
    if opts.region != "":
        log("ARG", fmt"Regions: {opts.region}")
    if opts.impact_order != "":
        log("ARG", fmt"Custom impact order: {opts.impact_order}")
    if opts.csq_column != "":
        log("ARG", fmt"Additional CSQ columns: {opts.csq_column}")
    if opts.info_column != "":
        log("ARG", fmt"Additional INFO fields: {opts.info_column}")
    if opts.most_severe:
        active_filters.add("most_severe")
        log("ARG", fmt"Most severe active with mode {opts.mode}")
    if opts.min_exp != "":
        active_filters.add("expression")
        log("ARG", fmt"Min expression: {opts.min_exp}")
        log("ARG", fmt"Expression data: {opts.exp_data}")
        log("ARG", fmt"Tissues: {opts.tissues}")
    if opts.transcripts != "":
        active_filters.add("transcripts")
        log("ARG", fmt"Transcripts: {opts.transcripts}")
    if opts.min_impact != "":
        active_filters.add("min_impact")
        log("ARG", fmt"Min impact: {opts.min_impact}")
    log("ARG", fmt"Active filters: {active_filters}")
    log("ARG", fmt"Tagging congiguration: {opts.var_tagging_json}")