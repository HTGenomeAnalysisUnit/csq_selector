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
type CsqFieldIndexes* = object
  gene_id*: int
  gene_symbol*: int
  consequence*: int
  transcript*: int
  columns*: OrderedTableRef[string, int]

# Store configuration to parse and filter csqs
type Config* = object
  csq_field_name*: string
  csq_field_idxs*: CsqFieldIndexes
  min_impact*: string
  allowed_tx*: HashSet[string]
  ranked_exp*: HashSet[string]
  most_severe*: bool
  most_expressed*: bool
  group_by_gene*: bool
  tagging_config*: JsonNode
  csq_output_fields*: seq[string]
  tx_vers_re*: Regex

#Store impact information
type Impact* = object
  gene_id: string
  gene_symbol: string
  transcript: string
  transcript_version: string
  impact: string
  order: int
  csq_string: string
  csq_fields: TableRef[string, string]
  tag_suffix: HashSet[string]
  scores_suffix: int

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

#Given a seq of csq string returns csqs per gene in a table
proc get_csq_bygene(csqs: seq[Impact]): Table[string, seq[Impact]] {.inline.} =
  for c in csqs:
    if result.hasKeyOrPut(c.gene_id, @[c]): 
      result[c.gene_id].add(c)

#Given a seq of csq strings and an impact order returns the highest severity consequence
proc get_highest_impact(csqs: seq[Impact]): Impact =
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
proc get_most_severe*(csqs:seq[Impact], by_gene: bool): seq[Impact] =
  if by_gene:
    var by_gene_csqs = csqs.get_csq_bygene
    for gene_id, gene_csqs in by_gene_csqs.pairs():
      var imp = gene_csqs.get_highest_impact
      result.add(imp)

  else:
    var imp = csqs.get_highest_impact
    result.add(imp)

#Split consequences from ANN/CSQ/BCSQ and returns a list of csq as Impact object
proc split_csqs*(v:Variant, config: Config, impact_order: TableRef[string, int]): (int, seq[Impact]) =
  let field_indexes = config.csq_field_idxs
  var max_impact_order = 99
  if config.min_impact != "":
    max_impact_order = impact_order[config.min_impact]
  
  var 
    s = ""
    csqfield_missing = 0
    parsed_impacts: seq[Impact]
  
  if v.info.get(config.csq_field_name, s) != Status.OK: 
    csqfield_missing = 1
  else:
    let csqs = s.split(',')  
    for csq in csqs:
      var toks = csq.split('|')
      var tx = toks[field_indexes.transcript]
      tx = tx.cleanTxVersion(config.tx_vers_re)
      if config.allowed_tx.len > 0 and not config.allowed_tx.contains(tx): continue
      if config.ranked_exp.len > 0 and not config.ranked_exp.contains(tx): continue
      for i in toks[field_indexes.consequence].split('&'):
        var x: Impact

        var impact = i.toLowerAscii
        if impact.endsWith("_variant"):
          impact = impact[0..impact.high - 8]
        
        var val:int
        try:
          val = impact_order[impact]
        except:
          log("WARNING", fmt"unknown impact '{impact}' from csq '{csq}' please report the variant to the developers")
          val = impact_order.getOrDefault("PROTCHANGE_CUTOFF", 1)

        if val > max_impact_order: continue

        var scoring = 0
        var tags: HashSet[string]
        
        var csq_classes_config = %* {}
        if config.tagging_config.hasKey("csq_classes"): csq_classes_config = config.tagging_config["csq_classes"]
        if not csq_classes_config.hasKey(impact): continue

        let csq_classes_impact = csq_classes_config[impact]
        var tagging_config = %* {}
        var scoring_config = %* {}
        if csq_classes_impact.hasKey("tagging"): tagging_config = csq_classes_impact["tagging"]
        if csq_classes_impact.hasKey("scoring"): scoring_config = csq_classes_impact["scoring"]
        
        if tagging_config.len > 0:
          for tag_config in tagging_config.items:
            let tag = tag_config["tag"].getStr()
            if tag_config.hasKey("csq_field"):
              for k in tag_config["csq_field"].keys:
                let field_name = tag_config["csq_field"][k].getStr()
                let tag_obj = tag_config["csq_field"][k]
                if tag_obj["value"].kind == JFloat or tag_obj["value"].kind == JInt:
                  let field_value = toks[field_indexes.columns[field_name]].parseFloat()
                  let score_threshold = tag_obj["value"].getFloat()
                  if compare_values(field_value, score_threshold, tag_obj{"operator"}.getStr(">")):
                    tags.incl(tag)
                elif tag_obj["value"].kind == JString:
                  let field_value = toks[field_indexes.columns[field_name]]
                  let flag_value = tag_obj["value"].getStr()
                  if field_value == flag_value:
                    tags.incl(tag)
                elif tag_obj["value"].kind == JBool:
                  let field_value = toks[field_indexes.columns[field_name]]
                  let flag_value = tag_obj["value"].getBool()
                  if (field_value != "") == flag_value:
                    tags.incl(tag)
          
            if tag_config.hasKey("info"):
              for k in tag_config["info"].keys:
                let tag_obj = tag_config["csq_field"][k]
                if tag_obj["value"].kind == JFloat or tag_obj["value"].kind == JInt:
                  let score_threshold = tag_obj["value"].getFloat()
                  var info_value: seq[float32] 
                  if v.info.get(k, info_value) == Status.OK:
                    if compare_values(info_value[0], score_threshold, tag_obj{"operator"}.getStr(">")):
                      tags.incl(tag)
                elif tag_obj["value"].kind == JString:
                  let flag_value = tag_obj["value"].getStr()
                  var info_value: string
                  if v.info.get(k, info_value) == Status.OK:
                    if info_value == flag_value:
                      tags.incl(tag)
                elif tag_obj["value"].kind == JBool:
                  let flag_value = tag_obj["value"].getBool()
                  if v.info.has_flag(k) == flag_value:
                    tags.incl(tag)
              
        if scoring_config.len > 0:
          if scoring_config.hasKey("csq_field"):
            for k in scoring_config["csq_field"].keys:
              let field_name = scoring_config["csq_field"][k].getStr()
              let tag_obj = scoring_config["csq_field"][k]
              if tag_obj["value"].kind == JFloat or tag_obj["value"].kind == JInt:
                let field_value = toks[field_indexes.columns[field_name]].parseFloat()
                let score_threshold = tag_obj["value"].getFloat()
                if compare_values(field_value, score_threshold, tag_obj{"operator"}.getStr(">")):
                  scoring += 1
              elif tag_obj["value"].kind == JString:
                let field_value = toks[field_indexes.columns[field_name]]
                let flag_value = tag_obj["value"].getStr()
                if field_value == flag_value:
                  scoring += 1
              elif tag_obj["value"].kind == JBool:
                let field_value = toks[field_indexes.columns[field_name]]
                let flag_value = tag_obj["value"].getBool()
                if (field_value != "") == flag_value:
                  scoring += 1
          
          if scoring_config.hasKey("info"):
            for k in scoring_config["info"].keys:
              let tag_obj = scoring_config["csq_field"][k]
              if tag_obj["value"].kind == JFloat or tag_obj["value"].kind == JInt:
                let score_threshold = tag_obj["value"].getFloat()
                var info_value: seq[float32] 
                if v.info.get(k, info_value) == Status.OK:
                  if compare_values(info_value[0], score_threshold, tag_obj{"operator"}.getStr(">")):
                    scoring += 1
              elif tag_obj["value"].kind == JString:
                let flag_value = tag_obj["value"].getStr()
                var info_value: string
                if v.info.get(k, info_value) == Status.OK:
                  if info_value == flag_value:
                    scoring += 1
              elif tag_obj["value"].kind == JBool:
                let flag_value = tag_obj["value"].getBool()
                if v.info.has_flag(k) == flag_value:
                  scoring += 1

        x.gene_id = toks[config.csq_field_idxs.gene_id].cleanTxVersion(config.tx_vers_re)
        x.gene_symbol = toks[config.csq_field_idxs.gene_symbol]
        x.transcript = tx
        x.transcript_version = toks[config.csq_field_idxs.transcript]
        x.impact = impact
        x.order = val
        x.csq_string = csq
        x.tag_suffix = tags
        x.scores_suffix = scoring
        
        parsed_impacts.add(x)

        if config.most_severe:
          parsed_impacts = get_most_severe(parsed_impacts, config.group_by_gene)

  result = (csqfield_missing, parsed_impacts)

#Read header and set CSQ indexes for relevant fields
proc set_csq_fields_idx*(ivcf:VCF, field:string, csq_columns: seq[string]= @[]): (string, CsqFieldIndexes) =
  var field_indexes: CsqFieldIndexes
  var csq_field: string
  field_indexes.gene_id = -1
  field_indexes.gene_symbol = -1
  field_indexes.consequence = -1
  field_indexes.transcript = -1
  field_indexes.columns = newOrderedTable[string, int]()

  # try to get the requested field, but iterate through other known csq fields
  # as a backup. sometimes, snpEff, for example will not add it's ANN or EFF
  # field to the header given an empty VCF.
  var desc: string
  for tryfield in [field, "CSQ", "BCSQ", "ANN"]:
    try:
      desc = ivcf.header.get(tryfield, BCF_HEADER_TYPE.BCF_HL_INFO)["Description"]
      log("INFO", fmt"Reading gene consequences from {field}")
      csq_field = tryfield
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
  for (i, v) in adesc.pairs:
    field_indexes.columns[v] = i

  # check if all requested columns are present
  for cq in csq_columns:
    let idx = adesc.find(cq.toUpperAscii)
    if idx == -1:
      raise newException(KeyError, fmt"requested csq column '{cq}' not found in {adesc}")

  #ANN: symbol=GENE_NAME, id=GENE_ID, transcript=FEATURE_ID, csq=ANNOTATION
  #BCSQ: symbol=GENE, id=N/A, transcript=TRANSCRIPT, csq=CONSEQUENCE
  #CSQ: symbol=SYMBOL, id=GENE, transcript=FEATURE, csq=CONSEQUENCE
  for check in ["GENE_ID", "GENE"]:
    field_indexes.gene_id = adesc.find(check)
    if field_indexes.gene_id != -1: break
  for check in ["SYMBOL", "GENE", "GENE_NAME"]:
    field_indexes.gene_symbol = adesc.find(check)
    if field_indexes.gene_symbol != -1: break
  for check in ["CONSEQUENCE", "ANNOTATION"]:
    field_indexes.consequence = adesc.find(check)
    if field_indexes.consequence != -1: break
  for check in ["FEATURE", "TRANSCRIPT", "FEATURE_ID"]:
    field_indexes.transcript = adesc.find(check)
    if field_indexes.transcript != -1: break

  if field_indexes.gene_id == -1:
    log("FATAL", fmt"unable to find gene id field in {field}")    
    quit "", QuitFailure
  if field_indexes.consequence == -1:
    log("FATAL", fmt"unable to find consequence field in {field}")
    quit "", QuitFailure
  if field_indexes.transcript == -1:
    log("FATAL", fmt"unable to find transcript field in {field}")
    quit "", QuitFailure

  result = (csq_field, field_indexes)

#Create csq output strings according to format
proc get_csq_string*(csqs: seq[Impact], csq_columns: seq[string], format: string): seq[string] =
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
        for c in csq_columns:
          line.add(x.csq_fields.getOrDefault(c))
        result.add(line.join("\t"))
      of "vcf":
        result.add(x.csq_string)
      of "rarevar_set":
        var impact_str = x.impact
        if x.tag_suffix.len > 0:
          for t in x.tag_suffix: impact_str.add("-" & t)          
        if x.scores_suffix > 0:
          impact_str = fmt"{impact_str}-{x.scores_suffix}"
        result.add([x.gene_id, impact_str].join("\t"))
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
iterator make_set_string*(gene_set: Table[string, Gene_set]): string {.closure.} =
  var n = 0
  let interval = 1000
  for gene_id, gene_values in gene_set.pairs():
    n += 1
    if n < 10:
      log("INFO", fmt"{gene_values.vars.len} variants in setlist for gene {gene_id}. Reported for the first 10 genes")
    if floorMod(n, interval) == 0:
      log("INFO", fmt"{n} gene sets processed")
    yield [gene_id, gene_values.chrom, $gene_values.position, gene_values.vars.join(",")].join("\t")

#Given a seq of csq strings and an exp rank returns the highest severity consequence
proc get_highest_expressed(csqs: seq[Impact], gene_fields:CsqFieldIndexes, exp_order: seq[string], allow_miss: bool = true): seq[Impact] =
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
proc get_most_expressed*(csqs:seq[Impact], gene_fields:CsqFieldIndexes, exp: Ranked_exp, by_gene: bool, allow_miss: bool = true): seq[Impact] =  
  if by_gene:
    var by_gene_csqs = csqs.get_csq_bygene
    for gene_id, gene_csqs in by_gene_csqs.pairs():
      var gene_tx_rank = exp.by_gene.getOrDefault(gene_id, @[])
      var imp = gene_csqs.get_highest_expressed(gene_fields, gene_tx_rank, allow_miss)
      result.add(imp)

  else:
    var imp = csqs.get_highest_expressed(gene_fields, exp.global, allow_miss)
    result.add(imp)

#Given a seq of csq strings and a max impact returns csqs more severe than max impact
proc filter_impact(csqs: seq[Impact], impact_order: TableRef[string, int], max_impact: string): seq[Impact] =
  #Max_impact is expected to be in the impact_order table
  let max_impact_order = impact_order[max_impact]
    
  for c in csqs:
    if c.order < max_impact_order:
      result.add(c)

#Process variant and returns the
proc get_filtered_csqs*(csqs:seq[Impact], gene_fields:CsqFieldIndexes, impact_order: TableRef[string, int], max_impact: string, by_gene: bool): seq[Impact] = 
  if by_gene:
    var by_gene_csqs = csqs.get_csq_bygene
    for gene_id, gene_csqs in by_gene_csqs.pairs():
      var imp = gene_csqs.filter_impact(impact_order, max_impact)
      result.add(imp)

  else:
    var imp = csqs.filter_impact(impact_order, max_impact)
    result.add(imp)