import std/parseopt
import strutils
import strformat
import sequtils
import tables
import parseutils
import os
import ./impact_order

type ArgDefinition = object 
    default: string
    required: bool
    choices: seq[string]
    help: string
    is_flag: bool

type ArgNames = tuple
    short: string
    long: string

type ParseSchema* = object
    main_name: string 
    options: Table[ArgNames, ArgDefinition]
    #args: seq[ArgDefinition]

type ParsedValues* = object
    args: seq[string]
    opts: Table[string, string]
    flags: Table[string, bool]

template initParser*(p: var ParseSchema, title: string = "", body: untyped): untyped =
    if title == "": 
        p.main_name = os.getAppFileName()
    else:
        p.main_name = title
    
    p.addFlag("-h", "--help", "Show help message and exit")

    body

proc printHelp*(args: ParseSchema) {.discardable.} =
    echo "=== HELP MESSAGE ==="
    echo &"{args.main_name}\n"

    var 
        tags: string
        expected_value: string
        is_required: string
        arguments: seq[string]
        flags: seq[string]

    for k, v in args.options.pairs():
        if k.long != "" and k.short != "": 
            tags = &"-{k.short}, --{k.long}"
        else:
            tags = (if k.short != "": &"-{k.short}" else: &"--{k.long}")
        
        if v.is_flag: 
            expected_value = ""
            flags.add(&"{tags}\n\t{v.help}\n")
        else:
            if v.choices.len > 0: expected_value = $v.choices else: expected_value = "value"
            if v.required: is_required = "required" else: is_required = "optional" 
            arguments.add(&"{tags} {expected_value}\t{is_required}\n\t{v.help}\n")

    if arguments.len > 0:    
        echo "Arguments description:"
        echo arguments.join("\n")

    if flags.len > 0:
        echo "Flags description:"
        echo flags.join("\n")

    echo "=== END HELP MESSAGE ==="
    quit QuitSuccess

proc add(args: var ParseSchema, short: string, long: string, values: ArgDefinition) {.discardable.} =
    var names: ArgNames
    names.short = short.replace("-", "")
    if names.short.len > 1:
        raise newException(IndexDefect, fmt"Short option name must be one character long '{names.short}'")
    names.long = long.replace("--", "")
    if names.long.len == 1:
        raise newException(IndexDefect, fmt"Long option name must longer than one character '{names.long}'")

    for (s, l) in args.options.keys():
        if (names.short == s and names.short != ""):
            raise newException(KeyError, fmt"Conflicting options {names.short}")
        if (names.long == l and names.long != ""):
            raise newException(KeyError, fmt"Conflicting options {names.long}")
    
    args.options[names] = values
      
proc addOption*(args: var ParseSchema, short: string = "", long: string = "", required: bool = false, default: string = "", help: string = "", choices: seq[string] = @[]) {.discardable.} =
    var values: ArgDefinition
    values.default = default
    values.required = required
    values.choices = choices
    values.help = help
    values.is_flag = false

    args.add(short, long, values)

proc addFlag*(args: var ParseSchema, short: string = "", long: string = "", help: string = "") {.discardable.} =
    var values: ArgDefinition
    values.default = ""
    values.required = false
    values.help = help
    values.is_flag = true

    args.add(short, long, values)
  
#proc addArg*(args: var ParseSchema, n: int = 1, datatype: string = "string", help: string = "") {.discardable.} =
#    var values: ArgDefinition
#    values.default = ""
#    values.help = help
#    values.is_flag = false
#    values.datatype = datatype
#
#    args.args.add(values)

proc parseOptions*(args: var ParseSchema, cmdline: string = "", ): ParsedValues =
    var
        longnoval: seq[string]
        shortnoval: set[char]
        #shortchar: char

    #Get flags so we can set which arguments expects no values
    for k, v in args.options.pairs():
        if v.is_flag: 
            if k.long != "": longnoval.add(k.long) 
            if k.short != "": #and k.short.parseChar(shortchar, 0) == 1:
                shortnoval.incl(cast[char](k.short)) 

    #Parse command line arguments
    var p = initOptParser(cmdline, shortNoVal=shortnoval, longNoVal=longnoval)
    var cmdArgsTable: Table[string, string]
    var cmdLineKeys: seq[string]
    while true:
        p.next()
        case p.kind
        of cmdEnd: break
        of cmdShortOption:
            cmdArgsTable[p.key] = p.val
            cmdLineKeys.add(p.key)
        of cmdLongOption:
            cmdArgsTable[p.key] = p.val
            cmdLineKeys.add(p.key)
        of cmdArgument:
            result.args.add(p.val)

    #Check if help is called and show help message
    if cmdArgsTable.contains("help") or cmdArgsTable.contains("h"):
        args.printHelp()

    # Show the default CSQ order and exit
    if cmdArgsTable.contains("show_csq_order"):
        quit show_impact_order(default_order), QuitSuccess

    #Check the same option is not seen twice
    let unique = cmdLineKeys.deduplicate()
    if unique.len != cmdLineKeys.len:
        raise newException(IndexDefect, fmt"Duplicate options found {$cmdLineKeys}")
    
    #Using the ParseSchema and table or args from command line, assign values
    for k, v in args.options.pairs():        
        var parsed_keys: seq[string]
        if cmdArgsTable.hasKey(k.short): parsed_keys.add(k.short)
        if cmdArgsTable.hasKey(k.long): parsed_keys.add(k.long)
        if parsed_keys.len > 1:
            raise newException(KeyError, fmt"Conflicting options. You specified both {parsed_keys[0]} and {parsed_keys[1]} for the same option")
        
        var parsed_key = (if cmdArgsTable.hasKey(k.short): k.short else: k.long)

        if cmdArgsTable.hasKey(parsed_key):
            if v.is_flag:
                if k.short != "": result.flags[k.short] = true
                if k.long != "": result.flags[k.long] = true
            else:
                if v.choices.len > 0 and cmdArgsTable[parsed_key] notin v.choices:
                    raise newException(ValueError, fmt"Invalid value {cmdArgsTable[parsed_key]} for option {parsed_key}. Allowed values: {v.choices}")
                if k.short != "": result.opts[k.short] = cmdArgsTable[parsed_key]
                if k.long != "": result.opts[k.long] = cmdArgsTable[parsed_key]

        else:
            if v.is_flag: 
                if k.long != "": result.flags[k.long] = false
                if k.short != "": result.flags[k.short] = false
                continue
            if v.default == "" and v.required:
                raise newException(ValueError, fmt"Required option -{k.short} / --{k.long} is not set")
            if v.default != "":
                if k.short != "": result.opts[k.short] = v.default
                if k.long != "": result.opts[k.long] = v.default
    
proc getOpt*[T](p: ParsedValues, y: string, x: var T) {.discardable.} =
    try:
        when $(x.typeof) == "string": 
            x = p.opts[y]
        elif $(x.typeof) == "int": 
            x = p.opts[y].parseInt()
        elif $(x.typeof) == "float": 
            x = p.opts[y].parseFloat()
        elif $(x.typeof) == "bool":
            x = p.flags[y]
        else:
            raise newException(TypeError, fmt"Invalid type of variable, only string, int, float or bool accepted")
    except KeyError:
        raise newException(KeyError, fmt"Argument {y} not found")

proc isSet*(p: ParsedValues, y: string): bool =
    try:
        result = p.flags[y]
    except KeyError:
        if y in p.opts: result = true else: result = false

proc `[]`*(p: ParsedValues, y: string): string =
    if y in p.opts:
        result = p.opts[y]
    elif y in p.flags:
        if p.flags[y]: result = "true" else: result = "false"
    else:
        result = ""

proc `[]`*(p: ParsedValues, y: int): string =
    try:
        result = p.args[y]
    except KeyError:
        raise newException(KeyError, fmt"Positional argument at index {y} not found")

proc getArgs*(p: ParsedValues): seq[string] =
    result = p.args

proc keys*(p: ParsedValues): seq[string] =
    for k in p.opts.keys(): result.add(k)
    for k in p.flags.keys(): result.add(k)

iterator pairs*(p: ParsedValues): (string, string) =
    for k in p.keys():
        yield (k, p[k])