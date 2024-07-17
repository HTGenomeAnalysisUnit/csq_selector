# Split ANN from SnpEff and select most severe CSQ
import hts
import json
import strformat
import strutils
import system
import tables
import re
import std/sets
from math import floorMod
from ./utils import log, cleanTxVersion
from ./tx_expression import Ranked_exp
# import sequtils

#Store indexes for gene annotations in the CSQ field seq
type GeneIndexes* = object
  gene_id*: int
  gene_symbol*: int
  consequence*: int
  transcript*: int
  csq_field*: string

  columns*: OrderedTableRef[string, int]

#Store impact information
type Impact* = object
  gene_id: string
  gene_symbol: string
  transcript: string
  transcript_version: string
  impact: string
  order: int
  csq_string: string
  pass_scores: int

type Gene_set* = object 
  chrom: string
  position: int64
  vars: seq[string]

proc `$`*(x: Impact): string =
  result = fmt"[gene_id: {x.gene_id}, gene_symbol: {x.gene_symbol}, transcript: {x.transcript}, impact: {x.impact}, order: {x.order}]"

proc `$`*(x: seq[Impact]): string =
  var csqs: seq[string]
  for i in x:
    csqs.add($i)
  result = csqs.join("; ")

proc compare_values(x: float32, y: float, operator: string): bool =
  case operator:
    of "<": result = x < y
    of "<=": result = x <= y
    of ">": result = x > y
    of ">=": result = x >= y
    of "==": result = x == y
    of "!=": result = x != y
    else: raise newException(ValueError, fmt"unknown operator: {operator}")

proc compare_values(x: int32, y: int, operator: string): bool =
  case operator:
    of "<": result = x < y
    of "<=": result = x <= y
    of ">": result = x > y
    of ">=": result = x >= y
    of "==": result = x == y
    of "!=": result = x != y
    else: raise newException(ValueError, fmt"unknown operator: {operator}")

#Split consequences from ANN/CSQ/BCSQ and returns a list of csq as Impact object
proc split_csqs*(v:Variant, csq_field_name:string, gene_fields:GeneIndexes, impact_order: TableRef[string, int], tx_vers_re: Regex, max_impact: string = "", allowed_tx: HashSet[string], ranked_exp: HashSet[string], scores: JsonNode): (int, seq[Impact]) =
  var max_impact_order = 99
  if max_impact != "":
    max_impact_order = impact_order[max_impact]
  
  var 
    s = ""
    csqfield_missing = 0
    parsed_impacts: seq[Impact]
  
  if v.info.get(csq_field_name, s) != Status.OK: 
    csqfield_missing = 1
  else:
    let csqs = s.split(',')  
    for tr in csqs:
      var toks = tr.split('|')
      var tx = toks[gene_fields.transcript]
      tx = tx.cleanTxVersion(tx_vers_re)
      if allowed_tx.len > 0 and not allowed_tx.contains(tx): continue
      if ranked_exp.len > 0 and not ranked_exp.contains(tx): continue
      for i in toks[gene_fields.consequence].split('&'):
        var x: Impact

        var impact = i.toLowerAscii
        if impact.endsWith("_variant"):
          impact = impact[0..impact.high - 8]
        
        var val:int
        try:
          val = impact_order[impact]
        except:
          log("WARNING", fmt"unknown impact '{impact}' from csq '{tr}' please report the variant to the developers")
          val = impact_order.getOrDefault("PROTCHANGE_CUTOFF", 1)

        if val > max_impact_order: continue

        var n_pass_scores = 0
        try:
          var csq_scores = scores[impact]
          for score_tag in csq_scores.keys:
            let score_obj =  csq_scores[score_tag]
            if score_obj["value"].kind == JFloat or score_obj["value"].kind == JInt:
              let score_threshold = score_obj["value"].getFloat()
              var score_value: seq[float32] 
              if v.info.get(score_tag, score_value) == Status.OK:
                if compare_values(score_value[0], score_threshold, score_obj{"operator"}.getStr(">")):
                  n_pass_scores += 1
            elif score_obj["value"].kind == JBool:
              let flag_value = score_obj["value"].getBool()
              if info.has_flag(score_tag) == flag_value:
                n_pass_scores += 1
        except:
          discard

        x.gene_id = toks[gene_fields.gene_id].cleanTxVersion(tx_vers_re)
        x.gene_symbol = toks[gene_fields.gene_symbol]
        x.transcript = tx
        x.transcript_version = toks[gene_fields.transcript]
        x.impact = impact
        x.order = val
        x.csq_string = tr
        x.pass_scores = n_pass_scores
        
        parsed_impacts.add(x)

  result = (csqfield_missing, parsed_impacts)

#Given a seq of csq string returns csqs per gene in a table
proc get_csq_bygene(csqs: seq[Impact], gene_fields:GeneIndexes): Table[string, seq[Impact]] {.inline.} =
  for c in csqs:
    if result.hasKeyOrPut(c.gene_id, @[c]): 
      result[c.gene_id].add(c)

#Read header and set CSQ indexes for relevant fields
proc set_csq_fields_idx*(ivcf:VCF, field:string, gene_fields: var GeneIndexes, csq_columns: seq[string]= @[]): seq[string] {.discardable.} =
  gene_fields.gene_id = -1
  gene_fields.gene_symbol = -1
  gene_fields.csq_field = field
  gene_fields.consequence = -1
  gene_fields.transcript = -1
  gene_fields.columns = newOrderedTable[string, int]()

  # try to get the requested field, but iterate through other known csq fields
  # as a backup. sometimes, snpEff, for example will not add it's ANN or EFF
  # field to the header given an empty VCF.
  var desc: string
  for tryfield in [field, "CSQ", "BCSQ", "ANN"]:
    try:
      desc = ivcf.header.get(tryfield, BCF_HEADER_TYPE.BCF_HL_INFO)["Description"]
      break
    except:
      if tryfield == field:
        log("WARNING", fmt"Didn't find {field} in header in {ivcf.fname} trying other fields")

  if desc == "":
    raise newException(KeyError, fmt"consequence field: {field} not found in header")
  # snpEff ##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations: 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO' ">

  var spl = (if "Format: '" in desc: "Format: '" else: "Format: ")
  if spl notin desc:
    spl = ": '"
  var adesc:seq[string]
  try:
    adesc = desc.split(spl)[1].split("'")[0].strip().strip(chars={'"', '\''}).multiReplace(("[", ""), ("]", ""), ("'", ""), ("*", "")).split("|")
  except IndexDefect:
    # format field description not as expected. return emptyr result and don't fill gene fields
    return result

  for v in adesc.mitems: v = v.toUpperAscii.strip()
  result = adesc

  for cq in csq_columns:
    gene_fields.columns[cq] = adesc.find(cq.toUpperAscii)
    if gene_fields.columns[cq] == -1:
      raise newException(KeyError, fmt"requested csq column '{cq}' not found in {adesc}")

  #ANN: symbol=GENE_NAME, id=GENE_ID, transcript=FEATURE_ID, csq=ANNOTATION
  #BCSQ: symbol=GENE, id=N/A, transcript=TRANSCRIPT, csq=CONSEQUENCE
  #CSQ: symbol=SYMBOL, id=GENE, transcript=FEATURE, csq=CONSEQUENCE
  for check in ["GENE_ID", "GENE"]:
    gene_fields.gene_id = adesc.find(check)
    if gene_fields.gene_id != -1: break
  for check in ["SYMBOL", "GENE", "GENE_NAME"]:
    gene_fields.gene_symbol = adesc.find(check)
    if gene_fields.gene_symbol != -1: break
  for check in ["CONSEQUENCE", "ANNOTATION"]:
    gene_fields.consequence = adesc.find(check)
    if gene_fields.consequence != -1: break
  for check in ["FEATURE", "TRANSCRIPT", "FEATURE_ID"]:
    gene_fields.transcript = adesc.find(check)
    if gene_fields.transcript != -1: break

  if gene_fields.gene_id == -1:
    log("FATAL", fmt"unable to find gene id field in {field}")    
    quit "", QuitFailure
  if gene_fields.consequence == -1:
    log("FATAL", fmt"unable to find consequence field in {field}")
    quit "", QuitFailure
  if gene_fields.transcript == -1:
    log("FATAL", fmt"unable to find transcript field in {field}")
    quit "", QuitFailure

#Create csq output strings according to format
proc get_csq_string*(csqs: seq[Impact], gene_fields: GeneIndexes, format: string): seq[string] =
  ## get the gene_names and consequences for each transcript.
  ## Adapt this to be able to output TSV format

  for x in csqs:    
    case format:
      of "tsv":
        var line = @[
          x.gene_id,
          x.gene_symbol,
          x.transcript_version,
          x.impact
        ]
        var toks = x.csq_string.split('|')
        for c, ci in gene_fields.columns:
          line.add(if ci < toks.len: toks[ci] else: "")
        result.add(line.join("\t"))
      of "vcf":
        result.add(x.csq_string)
      of "rarevar_set":
        if x.pass_scores == 0:
          result.add([x.gene_id, x.impact].join("\t"))
        else:
          result.add([x.gene_id, fmt"{x.impact}-{x.pass_scores}"].join("\t"))
      else:
        raise newException(ValueError, fmt"unknown output format: {format}")

# Given a seq of impacts for a variant, update a gene_set with variants corresponding to each gene
proc update_gene_set*(gene_set: var Table[string, Gene_set], v: Variant, csqs: seq[Impact], useid: bool = false) =
  var var_id = [$v.CHROM, $v.POS, v.REF, v.ALT[0]].join("_")
  if useid: 
    var_id = $v.ID
  
  for c in csqs:
    var gene_values = gene_set.getOrDefault(c.gene_id)
    gene_values.chrom = $v.CHROM
    if gene_values.position > v.POS or gene_values.position == 0:
      gene_values.position = v.POS
    gene_values.vars.add(var_id)
    gene_set[c.gene_id] = gene_values

# Given a gene_set generate a seq of strings representing gene sets in regenie format
proc make_set_string*(gene_set: Table[string, Gene_set]): seq[string] =
  var n = 0
  let interval = 1000
  for gene_id, gene_values in gene_set.pairs():
    n += 1
    if n < 10:
      log("INFO", fmt"{gene_values.vars.len} variants in setlist for gene {gene_id}. Reported for the first 10 genes")
    result.add([gene_id, gene_values.chrom, $gene_values.position, gene_values.vars.join(",")].join("\t"))
    if floorMod(n, interval) == 0:
      log("INFO", fmt"{n} gene sets processed")

#Given a seq of csq strings and an impact order returns the highest severity consequence
proc get_highest_impact(csqs: seq[Impact], gene_fields:GeneIndexes): Impact =
  #Change this to process over a seq of strings as csqs
  result.gene_id = "unknown"
  result.gene_symbol = "unknown"
  result.transcript = "uknown"
  result.impact = "unknown"
  result.order = 99
  result.csq_string = "unknown"
  
  for c in csqs:
    if c.order < result.order:
      result = c

#Process variant and returns the most severe impact across all genes or by gene
proc get_most_severe*(csqs:seq[Impact], gene_fields:GeneIndexes, by_gene: bool): seq[Impact] =
  if by_gene:
    var by_gene_csqs = csqs.get_csq_bygene(gene_fields)
    for gene_id, gene_csqs in by_gene_csqs.pairs():
      var imp = gene_csqs.get_highest_impact(gene_fields)
      result.add(imp)

  else:
    var imp = csqs.get_highest_impact(gene_fields)
    result.add(imp)

#Given a seq of csq strings and an exp rank returns the highest severity consequence
proc get_highest_expressed(csqs: seq[Impact], gene_fields:GeneIndexes, exp_order: seq[string], allow_miss: bool = true): seq[Impact] =
  var imp: Impact
  var tx_seq = exp_order
  for c in csqs:
    let transcript_rank = tx_seq.find(c.transcript)
    if transcript_rank == -1:
      if allow_miss: result.add(c)
      continue
    else:
      tx_seq = tx_seq[0..transcript_rank]
      imp = c
  result.add(imp)

#Process variant and returns the csq for the most expressed transcript across all genes or by gene
proc get_most_expressed*(csqs:seq[Impact], gene_fields:GeneIndexes, exp: Ranked_exp, by_gene: bool, allow_miss: bool = true): seq[Impact] =  
  if by_gene:
    var by_gene_csqs = csqs.get_csq_bygene(gene_fields)
    for gene_id, gene_csqs in by_gene_csqs.pairs():
      var gene_tx_rank = exp.by_gene.getOrDefault(gene_id, @[])
      var imp = gene_csqs.get_highest_expressed(gene_fields, gene_tx_rank, allow_miss)
      result.add(imp)

  else:
    var imp = csqs.get_highest_expressed(gene_fields, exp.global, allow_miss)
    result.add(imp)

#Given a seq of csq strings and a max impact returns csqs more severe than max impact
proc filter_impact(csqs: seq[Impact], gene_fields:GeneIndexes, impact_order: TableRef[string, int], max_impact: string): seq[Impact] =
  #Max_impact is expected to be in the impact_order table
  let max_impact_order = impact_order[max_impact]
    
  for c in csqs:
    if c.order < max_impact_order:
      result.add(c)

#Process variant and returns the
proc get_filtered_csqs*(csqs:seq[Impact], gene_fields:GeneIndexes, impact_order: TableRef[string, int], max_impact: string, by_gene: bool): seq[Impact] = 
  if by_gene:
    var by_gene_csqs = csqs.get_csq_bygene(gene_fields)
    for gene_id, gene_csqs in by_gene_csqs.pairs():
      var imp = gene_csqs.filter_impact(gene_fields, impact_order, max_impact)
      result.add(imp)

  else:
    var imp = csqs.filter_impact(gene_fields, impact_order, max_impact)
    result.add(imp)