import tables
import strutils
import os
import std/algorithm

const order_x = staticRead("./default-order.txt")

proc adjustOrder*(order: string): TableRef[string, int] =
  result = newTable[string, int](16)
  for o in order.strip().split("\n"):
    if o[0] == '#' or o.len == 0: continue
    var n: string
    if o.endsWith("CUTOFF"): n = o else: n = o.toLowerAscii.strip()
    # keep IMPACTFUL_CUTOFF as-is, all other impacts are lower-cased
    #if n[0] == 'i' and n in ["impact_cutoff", "impactful_cutoff"]: n = "IMPACTFUL_CUTOFF"
    #if o.len == 0: continue
    if n.endsWith("_variant"):
      n = n[0..n.high - 8]
    result[n] = result.len

var default_order* = adjustOrder(order_x)
if getEnv("SELECTOR_CSQ_ORDER") != "":
  if not fileExists(getEnv("SELECTOR_CSQ_ORDER")):
    raise newException(IOError, "[error] couldn't open file at:" & getEnv("SELECTOR_CSQ_ORDER") & " specified by env var 'SELECTOR_CSQ_ORDER'")
  default_order = adjustOrder(getEnv("SELECTOR_CSQ_ORDER").readFile)

proc show_impact_order*(order: TableRef[string, int]): string =
  var csq_list: seq[(int,string)]
  for k, i in order.pairs():
    csq_list.add((i,k))
  csq_list.sort()

  var result_seq: seq[string]
  for (i, val) in csq_list:
    result_seq.add($i & " - " & val)
  result = result_seq.join("\n")
