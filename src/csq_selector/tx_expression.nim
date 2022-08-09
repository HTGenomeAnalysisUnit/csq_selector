import zip/gzipfiles
import tables
import strutils
import strformat
import sequtils
import re
import std/enumerate
import std/algorithm
import std/sets
import ./utils

type Tx_exp* = tuple
    exp: float
    transcript: string

type Ranked_exp* = object
    global*: seq[string]
    by_gene*: Table[string, seq[string]]

#This take Tissue exp and return list of transcript id sorted by expression
proc sort_exp(exp: var seq[Tx_exp]): seq[string] =
    exp.sort()
    exp.reverse()
    result = exp.map(proc(x: Tx_exp): string = x.transcript)
    #for x in exp:
    #    if result.find(x.transcript) == -1:
    #        result.add(x.transcript)

#Read expression matrix and return a table with key = gene_id and seq containing transcript ids sorted by exp in relevant tissues
proc read_expression* (filename: string, transcript_col: int, gene_col: int, start_col: int, tissues: seq[string], min_exp: float, reg_exp_version: Regex): seq[string] =
    #var gene_exp: Table[string,seq[Tx_exp]]
    var global_exp: seq[Tx_exp]
    var n_tx, n_data_points = 0
    
    log("INFO", fmt"Reading expression matrix from {filename}")
    let fs = (if filename.endsWith(".gz"): newGzFileStream(filename) else: newFileStream(filename))

    let header = fs.readLine().split("\t")
    let start_idx = start_col
    log("INFO", fmt"{header.high - start_idx} tissues columns detected")
    
    while not fs.atEnd():
        let line = fs.readLine().split("\t")
        let tx = line[transcript_col]
        let transcript_id = tx.cleanTxVersion(reg_exp_version)
        n_tx += 1
        
        var exp_values = line[start_idx..line.high]

        for i, x in enumerate(exp_values):
            let t = header[start_idx + i]
            if tissues.find(t) == -1: continue
            var exp_value: float
            try:
                exp_value = parseFloat(x)
                n_data_points += 1
            except:
                log("WARNING", fmt"Could not parse '{x}' to float value for transcript {transcript_id} in tissue {t}")
            if exp_value > min_exp: global_exp.add((exp_value, transcript_id))
            #if gene_exp.hasKeyOrPut(line[gene_col], @[(exp_value, transcript_id)]):
            #    gene_exp[line[gene_col]].add((exp_value, transcript_id))
    
    log("INFO", fmt"{n_tx} transcripts and {n_data_points} exp values read. Sorting expression...")
    #for gene_id in gene_exp.keys():
    #    let id = gene_id.cleanTxVersion(reg_exp_version)
    #    var exp_data = gene_exp[gene_id]
    #    result.by_gene[id] = sort_exp(exp_data)
    
    #result.global = sort_exp(global_exp)
    result = sort_exp(global_exp)

    log("INFO", fmt"{result.len} transcripts selected with expesssion above {min_exp}")

proc read_exp* (filename: string, transcript_col: int, gene_col: int, start_col: int, tissues: seq[string], min_exp: float, reg_exp_version: Regex): HashSet[string] =
    #var gene_exp: Table[string,seq[Tx_exp]]
    #var global_exp: seq[Tx_exp]
    var n_tx, n_data_points = 0
    
    log("INFO", fmt"Reading expression matrix from {filename}")
    let fs = (if filename.endsWith(".gz"): newGzFileStream(filename) else: newFileStream(filename))

    let header = fs.readLine().split("\t")
    let start_idx = start_col
    log("INFO", fmt"{header.high - start_idx} tissues columns detected")
    
    while not fs.atEnd():
        let line = fs.readLine().split("\t")
        let tx = line[transcript_col]
        let transcript_id = tx.cleanTxVersion(reg_exp_version)
        n_tx += 1
        
        var exp_values = line[start_idx..line.high]

        for i, x in enumerate(exp_values):
            let t = header[start_idx + i]
            if tissues.find(t) == -1: continue
            var exp_value: float
            try:
                exp_value = parseFloat(x)
                n_data_points += 1
            except:
                log("WARNING", fmt"Could not parse '{x}' to float value for transcript {transcript_id} in tissue {t}")
            if exp_value > min_exp: 
                result.incl(transcript_id)
    
    log("INFO", fmt"{n_tx} transcripts and {n_data_points} exp values read")
    #for gene_id in gene_exp.keys():
    #    let id = gene_id.cleanTxVersion(reg_exp_version)
    #    var exp_data = gene_exp[gene_id]
    #    result.by_gene[id] = sort_exp(exp_data)
    
    #result.global = sort_exp(global_exp)
    #result = sort_exp(global_exp)

    log("INFO", fmt"{result.len} transcripts selected with expesssion above {min_exp}")