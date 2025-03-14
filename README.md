# CSQ SELECTOR

Small utility to manipulate gene consequences annotations from snpEff, VEP or bcftools.

Given a VCF file, it allows to filter consequences according to:

- min consequence (like missense, etc)
- transcripts with a min level of expression across specified conditions/tissues (expression values are loaded from a tab-separated file)
- a list of transcripts of interest

Additionally, for each variant one can decide to output only the most severe consequence or the most severe consequences per gene (in case a variant affect more than one gene).

The user can combine these filtering strategy to achieve precise control on the output consequences (see filter logic below). 

The tool can output variants in either VCF, TSV or regenie set variant format, thus allowing for high flexibility in the downstream analysis.

## Filter logic

Filters are applied in this order: 

1. consequences are filtered based on allowed transcripts (intersection of the [custom transcripts list](#set-transcripts-or-regions-of-interest) and transcripts [selected based on expression values](#expression-filter))
2. remaining consequences are filtered retaining only those more severe than the [specified min_impact](#min-impact-filter)
3. the most severe consequence is then selected if the ["most severe" flag is active](#most-severe-filter)

In this way, it is possible to accomodate sophisticated scenarios. For example one can select only the most severe consequence across transcripts expressed in the tissue of interest, but only when consequence is at least exonic.

The use of the `--mode` option can further refine the most severe filter. When set to `by_gene`, the most severe consequence is selected for each gene separately. This is useful when a variant affects more than one gene and the user is interested in the most severe consequence for each gene.

## Installation

Compiled binaries are provided with static linked htslib so they should work out of the box.

Otherwise, you can clone the repo and compile yourself. CSQ selector depends on

- nim >= 1.4.8
- htslib >= 1.10

And it uses the following nim packages

- hts-nim >= 0.3.21
- zip >= 0.3.1
- argparse >= 3.0.0

If you have nimble installed you can compile easily with:

```bash
git clone https://github.com/HTGenomeAnalysisUnit/csq_selector.git
cd csq_selector
nimble build 
```

If you have singularity installed, you can compile a static sel-contained build using the `nim_compile.sh` script

## Usage

A minimal example to select only the most severe consequence for each variant:

```bash
./csq_selector --vcf test/GRCh37_small_test.vcf.gz --out_format vcf --out test/out.vcf.gz --most_severe
```

You can use the following to see the default consequences order:

```bash
./csq_selector show_csq_order
```

Complete list of options:

```bash
Usage:
  CSQ Selector [options] 

Options:
  -h, --help
  -i, --vcf=VCF              path to VCF/BCF
  -o, --out=OUT              Output file
  -O, --out_format=OUT_FORMAT
                             Output format Possible values: [vcf, tsv, rarevar_set] (default: vcf)
  -c, --csq_field=CSQ_FIELD  INFO field containing the gene name and impact. Usually CSQ, ANN or BCSQ (default: ANN)
  -g, --gene_rarevar_identified=GENE_RAREVAR_IDENTIFIED
                             Set which gene identifier to use when building rarevar_set output Possible values: [gene_id, gene_symbol] (default: gene_id)
  --csq_column=CSQ_COLUMN    CSQ sub-field(s) to extract (in addition to gene, impact, transcript) when output is TSV. A comma-separated list of sub-fields can be provided
  --info_column=INFO_COLUMN  INFO field(s) to include in the TSV output. A comma-separated list can be provided
  -e, --exp_data=EXP_DATA    path to expression file. Tab-separated table of transcript expression across conditions
  -n, --min_exp=MIN_EXP      Min expression value for a transcript to be considered expressed
  --tissues=TISSUES          List of tissues to select in the expression file. Comma-separated list or file with 1 tissue per line.
  -r, --region=REGION        Specify genomic region for processing. Format is chr[:start-end]. Comma-separated list of regions or file with 1 region per line.
  -t, --transcripts=TRANSCRIPTS
                             list of trancsripts to restrict consequences output. Comma-separated list or file with 1 tissue per line.
  --mode=MODE                Mode to apply most severe filter. Across all annotations or by gene Possible values: [all, by_gene] (default: all)
  -m, --most_severe          Only output the most severe consequence for each variant
  --keep_old_ann             Keep original annotation fields in the output VCF under ORIGINAL_ANN tag
  --min_impact=MIN_IMPACT    Impact threshold from the consequences order
  --impact_order=IMPACT_ORDER
                             ordering of impacts to override the default. See default using show_csq_order
  -j, --var_tagging_json=VAR_TAGGING_JSON
                             A JSON file describing schema for variant tagging
  --use_vcf_id               Use the ID field from the VCF as the variant ID in the rarevar_set output
  --filter                   Only output variants where there is at least one consequence after filtering
```

### Most severe filter

When the `-m`, `--most_severe` flag is active, only the most severe consequence for each variant will be selected. Using the `--mode` option it is possible to control if this filter is applied across all consequences or by gene. In the latter case, if a variant has consequences on more than one gene, the most severe consequence for each gene is selected.

### Min impact filter

It is possible to filter by consequence level using `--min_impact`. In this case you should specify the minimum accepted consequence (like missense). The default order of known consequences can be found using `./csq_selector show_csq_order`

You can specify a new order providing a text file with a list of consequences, one per line. This custon file can be specified using `SELECTOR_CSQ_ORDER` environment variable or passing it to `--impact_order` option. Note that the latter take precedence in case both are specified.

You can add special cutoff values in this list in the form of string ending in `_CUTOFF` and these can then be passed to `--min_impact`. By default, 2 cutoffs are set: PROTCHANGE_CUTOFF and EXONIC_CUTOFF.

**NB.** Note that if `_variant` is present this is removed before matching so for example both `missense_variant` and `missense` should work.

### Set transcripts or regions of interest

User can specify a list of transcripts of interest using `-t, --transcripts` or a list of regions using `-r, --regions`. Regions are in the format chr:start-stop and both options take either a comma-separated list or a file with one value per line.

**NB.** Transcript version is removed from both transcripts list and annotations from the VCF file to ensure they can match. This is done by regex `\.[0-9]+$`

### Expression filter

You can provide transcripts expression data using `--exp_data` option. This data should be a tab-separated file with header. Column 1 must contain transcript ids and column 2 gene ids. Then you should have 1 column per tissue/condition. A file with median TPM expression from GTeX v8 is provided in the repository in exp_data folder.

Then you can use `--tissues` to specify tissues of interest. Here you can specify a comma-separated list or a file with one value per line. Values must correspond to column headers in the expression file. Use quotes if column names contains spaces, like `--tissue "Whole blood,Brain cortex"`.

The min expression threshold is set using `--min_exp`. Only consequences affecting a transcript with expression above this threshold in at least on of the specified tissues are kept.

### Refine variant consequence names

When the output format is set to `rarevar_set` or `tsv` it is possible to configure value thresholds based on annotations in the INFO fields and in consequence blocks to refine variants categories. You can pass a JSON file to the `-j, --var_tagging_json` option to configure the scoring/tagging strategy.

An example of the JSON file is in `test/configuration.json`

#### Scoring and tagging

Within the `csq_classes` key you can define for each possible consequence (e.g. `missense`) a scoring strategy under the `scoring` key. Here you can define threshold for `info` or `csq_field` fields. Similarly, you can define a tagging strategy under the `tagging` key to add a tag to the variant consequence name when a criteria is met.

For example you can ask to have missense variant scored based on CADD and REVEL annotations so that the variant annotation for that impact will reflect the number of scores thresholds that are passed. For example instead of just `missense`, the annotation could be `missense-1` or `missense-2` if the variant has a score above threshold for one or two values, respectively.

You can also ask to apply a tag `pathogenic` to the consequence name if the variant has a INFO pathogenic flag. In this case the final impact can become `missense-pathogenic-1`, or `missense-pathogenic-2`, combining the scoring and tagging strategies.

At the moment, this functionality only works with numeric or flag annotations. For numeric annotation, it is possible to configure a threashold value and the requested comparison (`>, >=, <, <=, ==, !=`).

#### Aliasing

Under the `aliases` key you can define aliases for given impacts so that they are all renamed to a single name. This is useful when you want to group together different impacts under a single name. For example, you can group all coding impacts or all loss-of-function ones.

```json
"aliases": {
		"LoF": ["stop_lost","stop_gained","frameshift_variant","splice_acceptor_variant","splice_donor_variant","start_lost"]
}
```

## Output formats

The output format is controlled by the `-O`, `--out_format` option. An output file/prefix can be specified with `--out`

When `--out` is not specified, the output is writted to stdout when using `tsv` or `vcf` format, while for `rarevar_set` format an output prefix is mandatory. 

### VCF

In VCF format, the tools write a standard VCF, updating the original ANN/BCSQ/CSQ field with the selected consequences. If no consequence is left after filtering, the field is removed. The old annotations can be kept in a separate field named ORIGINAL_ANN using `--keep_old_ann` option.

### TSV

When TSV is select the tool will write a tab-separated text file with header. By default the following fields are written in the output: CHROM,POS, ID, REF, ALT, FILTER, GENE_ID, GENE_SYMBOL, TRANSCRIPT, CONSEQUENCE.

In TSV format, the user can use the `--csq-column` option to provide a comma-separated list of additional fields that need to be extracted from the consequence string adn these will be included as additional columns. These names must match those described in the header description for the ANN/BCSQ/CSQ field. 

If no consequence is left after filtering, the line is omitted.

### Rarevar set

This output format generate 2 output files that can be used when performing rare variants analysis with [regenie](https://rgcgithub.github.io/regenie/). The first file (`.annot`) contains variants impact annotations (variant id, gene id, impact), while the second file (`.setlist`) contains the list of variants associated to each gene (gene id, chrom, position, variant id list).

By default the variant id created as `CHROM:POS:REF:ALT`, but one can set the `--use_vcf_id` flag to use the ID reported in the VCF instead.

The format of these files is described in more detail in the [regenie docs](https://rgcgithub.github.io/regenie/options/#annotation-input-files).

## Credits

Author: Edoardo Giacopuzzi

This tool uses [hts-nim](https://github.com/brentp/hts-nim) for VCF processing and the code to rank variant consequences is inspired from [a module in slivar](https://github.com/brentp/slivar/blob/master/src/slivarpkg/impact_order.nim) by Brent Pedersen.
