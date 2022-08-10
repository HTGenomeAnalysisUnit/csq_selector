# CSQ SELECTOR

Small utility to manipulate gene consequences annotations from snpEff, VEP or bcftools.

Given a VCF file, it allows to filter consequences according to:

- min consequence (like missense, etc)
- min level of expression of transcript across specified conditions/tissues (expression values are loaded from a tab-separated file)
- a list of transcripts of interest

Additionally, for each variant one can decide to output only the most severe consequence
or the most severe consequences per gene (in case a variant affect more than one gene).

Output is writted to file or stdout in TSV or VCF format. When TSV is select a list of additional CSQ columns to extract can be specified and variant with no remaining consequences after filtering are omitted.

## Installation

Compiled binaries are provided with static linked htslib so they should work out of the box.

Otherwise, you can clone the repo and compile yourself. CSQ selector depends on

- nim >= 1.4.8
- htslib >= 1.10

And it uses the following nim packages

- hts-nim >= 0.3.21
- zip >= 0.3.1

If you have nimble installed you can compile easily with:

```bash
git clone https://github.com/HTGenomeAnalysisUnit/csq_selector.git
cd csq_selector
nimble build 
```

If you have singularity installed, you can compile a static build using the `nim_compile.sh` script

## Usage

A minimal example to select only the most severe consequence for each variant:

```bash
./csq_selector --vcf test/GRCh37_small_test.vcf.gz --out_format vcf --out test/out.vcf.gz --most_severe
```

### Set transcripts or regions of interest

Additionally, you can specify a list of transcripts of interest using `-t, --transcripts` or a list of regions using `-r, --regions`. Regions are in the format chr:start-stop and both options take either a comma-separated list or a file with one value per line.

**NB.** Transcript version is removed from both transcripts list and annotations from the VCF file to ensure they can match. This is done by regex `\.[0-9]+$`

### Min impact filter

It is possible to filter by consequence level using `--min_impact`. In this case you should specify the minimum accepted consequence (like missense). The default order of known consequences can be found using `./csq_selector show_csq_order`

You can specify a new order providing a text file with a list of consequences, one per line. This custon file can be specified using `SELECTOR_CSQ_ORDER` environment variable or passing it to `--impact_order` option. Note that the latter take precedence in case both are specified.

You can add special cutoff values in this list in the form of string ending in `_CUTOFF` and these can then be passed to `--min_impact`. By default, 2 cutoffs are set: PROTCHANGE_CUTOFF and EXONIC_CUTOFF.

**NB.** Note that if `_variant` is present this is removed before matching so for example both `missense_variant` and `missense` should work.

### Expression filter

You can provide transcripts expression data using `--exp_data` option. This data should be a tab-separated file with header. Column 1 must contain transcript ids and column 2 gene ids. Then you should have 1 column per tissue/condition. A file with median TPM expression from GTeX v8 is provided in the repository in exp_data folder.

Then you can use `--tissues` to specify tissues of interest. Here you can specify a comma-separated list or a file with one value per line. Values must correspond to column headers in the expression file. Use quotes if column names contains spaces, like `--tissue "Whole blood,Brain cortex"`.

The min expression threshold is set using `--min_exp`. Only consequences affecting a transcript with expression above this threshold in at least on of the specified tissues are kept.

## Outputs

The output format can be set using `--out_format` option. You can use either `tsv` or `vcf` (default).

When using `tsv` format, the output is a tab-separated file with the following columns: CHROM,POS,ID,REF,ALT,FILTER,GENE_ID,GENE_SYMBOL,TRANSCRIPT,CONSEQUENCE + additional columns containig CSQ fields specified by `--csq-column`. Note that variants with no remaining consequences after filtering are omitted and the output will contain one consequence per line (so the number of output lines can be larger than the number of input variants)

When using `vcf` format, the output is a standard VCF file. All variants are omitted in this case and the CSQ field is removed if no consequence is left after filtering for a variant.

When not output file is specified (`--out`), the output is written to stdout.

## Credits

Author: Edoardo Giacopuzzi

This tool uses [hts-nim](https://github.com/brentp/hts-nim) for VCF processing and the most severe selection code is adapted from [a module in slivar](https://github.com/brentp/slivar/blob/master/src/slivarpkg/impact_order.nim) by Brent Pedersen.
