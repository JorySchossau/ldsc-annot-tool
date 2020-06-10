import strscans
import glob # dependency
import streams
import seqUtils
import strUtils
import strformat
import tables
import intsets # 1 set represents 1 category file of rsids
import os
import std/decls # ref assignment {.byaddr.}
import algorithm # sort
import simple_parseopt # dependency
import memfiles
import parsecsv
import math # nextPowerOfTwo

# fix for slow hashing of higher-order bits (until next release of nim, or use BTrees)
import hashes
proc hash(x: int): Hash = Hash(x * 0x49249549)

# could make this modifiable by CLI
var verbosity:int = 1

# 1 intset is 1 category file listing rsids
# the string is the filename sans extension
type
  Category = IntSet
  CategoryWithAnnotations = Table[int,string]

proc all_wildcards_expanded(paths:seq[string]|string): seq[string] =
  when paths is string:
    result = toSeq(walkGlob(paths))
  else:
    for path in paths:
      let expansion = toSeq(walkGlob(path))
      result = concat(result, expansion)

# checks the contents of the file
# for appropriate formatting
# errors out and quits if not right
proc check_category_format(file:string) =
  var
    rsid = 0
    annotation = 0.0
    validFormat = true
  for line in file.lines:
    if line.scanf("rs$i$.",rsid): break
    if line.scanf("rs$i$s$f$.",rsid,annotation): break
    validFormat = false
    break
  if not validFormat:
    stderr.writeLine("")
    stderr.writeLine(&"Error: {file} is not a category file containing rsids")
    stderr.writeLine("use --examples to generate example input files")
    quit(1)

# checks the contents of the file
# for appropriate formatting
# errors out and quits if not right
proc check_bim_format(file:string) =
  var
    chr:int
    rsid:int
    cm:float
    bp:int
  for line in file.lines:
    if not line.scanf("$i\trs$i\t$f\t$i",chr,rsid,cm,bp):
      stderr.writeLine("")
      stderr.writeLine(&"Error: {file} is not a bim file of the expected format")
      stderr.writeLine("use --examples to generate example input files")
      quit(1)
    break

# determine which type of category file this is
# if it has no annotations (rsid-col) or has annotations (rsid-col annotation-col)
proc has_annotations(file:string):bool =
  var rsid:int
  var annotation:float
  result = false
  for line in file.lines:
    if line.scanf("rs$i$s$f$.",rsid,annotation):
      result = true
    break

# read annotation category file and return representative data structures
# Category:intSet, CategoryWithAnnotations:Table[int,float]
proc to_category(file:string): (Category,CategoryWithAnnotations) =
  var category:Category
  var categoryWithAnnots:CategoryWithAnnotations
  var fileLen = 0
  for line in lines(memfiles.open(file)): inc fileLen
  if verbosity >= 1:
    echo &"adding {file} ({fileLen} lines)"
  # init the category
  var rsid:int
  var annotation:string
  check_category_format file
  if file.has_annotations:
    # init the annotation table
    categoryWithAnnots = initTable[int,string](nextPowerOfTwo fileLen)
    # scrape both rsid and annotation
    for line in file.lines:
      if line.scanf("rs$i$s$*$.",rsid,annotation):
        category.incl rsid
        categoryWithAnnots[rsid] = annotation
  else:
    # scrape only the rsids
    for line in file.lines:
      if line.scanf("rs$i$.",rsid):
        category.incl rsid
  result = (category,categoryWithAnnots)

# the main work proc
# builds off of data in bim file
# adding extra columns, 1 for each category
# putting a 1 or 0 if rsid is in category
proc write_annot(bimfile:string, allsnpfiles:seq[string], outdir:string) =

  # prepare output file
  create_dir outdir
  let outfilename = bimfile.split_file.name & ".annot"
  let destinationPath = outdir / outfilename

  check_bim_format bimfile
  var
    chr:int
    rsid:int
    cm:string
    bp:int
    lineN:int

  if verbosity >= 1:
    echo &"reading {bimfile}"

  var fileLen = 0
  for line in lines(memfiles.open(bimfile)): inc fileLen

  var
    allbps = newSeqOfCap[int](fileLen)
    allrsids = newSeqOfCap[int](fileLen)
    allcms = newSeqOfCap[string](fileLen)
    allchr:int
    columnsWithoutAnnots = initOrderedTable[string,seq[char]](nextPowerOfTwo allsnpfiles.len)
    columnsWithAnnots = initOrderedTable[string,seq[string]](nextPowerOfTwo allsnpfiles.len)
    columnsAnnotStatus = newSeqOfCap[bool](nextPowerOfTwo allsnpfiles.len)

  # load bim data
  for line in lines(memfiles.open(bimfile)):
    if line.scanf("$i\trs$i\t$*\t$i",chr,rsid,cm,bp):
      allrsids.add rsid
      allcms.add cm
      allbps.add bp
  allchr = chr
  # load each snp category data
  for snpfile in allsnpfiles:
    let (rsids,annotsByRsid) = snpfile.to_category
    let columnName = snpfile.splitPath.tail
    let has_no_annotations = annotsByRsid.len == 0
    columnsAnnotStatus.add not has_no_annotations
    var column:ptr seq[char]
    var columnWithAnnots:ptr seq[string]
    if has_no_annotations:
      columnsWithoutAnnots[columnName] = newSeqOfCap[char](fileLen)
      column = cast[ptr seq[char]](unsafeAddr columnsWithoutAnnots[columnName])
    else:
      columnsWithAnnots[columnName] = newSeqOfCap[string](fileLen)
      columnWithAnnots = cast[ptr seq[string]](unsafeAddr columnsWithAnnots[columnName])
    if has_no_annotations: # if rsids only
      for rsid in allrsids:
        if rsids.contains rsid: column[].add '1'
        else:                   column[].add '0'
    else:
      for rsid in allrsids:
        if rsids.contains rsid: columnWithAnnots[].add annotsByRsid[rsid]
        else:                   columnWithAnnots[].add "0"

  # write initial content to annot
  if verbosity >= 1:
    echo &"saving to {destinationPath}"
  var outfile = newFileStream(destinationPath, fmWrite)
  var outstring:string
  # write column headers
  outstring = "CHR\tBP\tSNP\tCM\tbase"
  for name in allsnpfiles.mapIt(it.splitPath.tail):
    outstring.add "\t" & name
  outstring.add "\n"
  # write column contents for each row
  for i in 0 ..< allrsids.len:
    outstring.add &"{allchr}\t{allbps[i]}\trs{allrsids[i]}\t{allcms[i]}\t1"
    for coln,name in allsnpfiles.mapIt(it.splitPath.tail):
      if columnsAnnotStatus[coln]: outstring.add &"\t{columnsWithAnnots[name][i]}"
      else:                        outstring.add &"\t{columnsWithoutAnnots[name][i]}"
    outstring.add "\n"
  writeFile(destinationPath, outstring)

const
  baseline_filename = "example.22.bim"
  category_filename = "example.H3K27ac_PGC2"
proc generate_example_files() =
  const
    baseline_contents = static_read "examples" / baseline_filename
    category_contents = static_read "examples" / category_filename
  baseline_filename.write_file baseline_contents
  category_filename.write_file category_contents

# helper for multi-in "if ["a","b","c"] in alist"
proc contains[T](a:openArray[T], b:openArray[T]):bool =
  for each_item in a:
    if each_item in b: return true
  return false

proc has_continuous_values(col:var seq[string]):bool =
  for value in col:
    if value != "0" and value != "1":
      return true
  return false

proc unmake(files:seq[string],outdir:string) =
  # indexed [5+colid,[rsid,annot]]
  # 5+ so it's same as original annot file column index
  var colsrsids = initTable[int,seq[int]]()
  var colsvals = initTable[int,seq[string]]()
  var colnames:seq[string]
  var file: CsvParser
  let SNPCOL = 2
  let FIRSTUSERCOL = 5
  var totalFileLens = 0
  for annotfile in files:
    for line in memSlices(memfiles.open(annotfile)): inc totalFileLens
  for annotfile in files:
    if verbosity >= 1:
      echo &"reading {annotfile}"
    # count number of lines in the annot file
    # read data from annot file
    file.open(annotfile, separator='\t')
    file.readHeaderRow()
    # init columns if necessary
    for coln,column in file.headers[FIRSTUSERCOL .. ^1]:
      let adjustedColId = coln+FIRSTUSERCOL
      if column notin colnames:
        colnames.add column
        colsrsids[adjustedColId] = newSeqOfCap[int](totalFileLens)
        colsvals[adjustedColId] = newSeqOfCap[string](totalFileLens)
    while file.readRow():
      for coln in FIRSTUSERCOL ..< file.headers.len:
        # [2..^1] gets int part of of "rs1234567"
        if file.row[coln] != "0" and file.row[coln] != "0.0":
          colsrsids[coln].add file.row[SNPCOL][2..^1].parseInt
          colsvals[coln].add file.row[coln]
  # prepare output file
  create_dir outdir
  for coln,column in colnames:
    let outfilename = column
    var outfile = newFileStream(outdir / outfilename, fmWrite)
    if verbosity >= 1:
      echo &"creating " & outdir / outfilename
    if colsvals[coln+FIRSTUSERCOL].has_continuous_values:
      for i in 0 ..< colsrsids[coln+FIRSTUSERCOL].len:
        outfile.write "rs" & $colsrsids[coln+FIRSTUSERCOL][i] & "\t" & colsvals[coln+FIRSTUSERCOL][i] & "\n"
    else:
      for rsid in colsrsids[coln+FIRSTUSERCOL]:
        outfile.write "rs" & $rsid & "\n"

const EXE = "annot"
const HELP_MESSAGE = &"""

Converts a collection of rsid categories into formatted .annot files for ldsc regression.

Usage
  {EXE} examples
  {EXE} make -b|--bimprefix PATH -s|--snpfiles FILES... [-o|--outdir OUTDIR]
  {EXE} unmake -a|--annotprefix PATH [-o|--outdir OUTDIR]

Example
  ./{EXE} -b test_baseline/subset. -s test_annots/FetalDHS_Trynka test_annots/H3K27ac_PGC2

Commands
  examples   create example input files
  make       construct annot files from category annotations
  unmake     deconstruct annot files and create category annotations

  examples
  [no options]

  make
  -b --bimprefix PATH    nonunique prefix path to all bim files
  -s --snpfiles FILES... all snp category files
  -o --outdir OUTDIR     output directory, optional (default 'output')

  unmake
  -a --annotprefix PATH    nonunique prefix path to all annot files
  -o --outdir OUTDIR     output directory, optional (default 'output')
  -h --help              this help message

Note that the snp file names will be used as category headers (names) in the resulting file. So if the files have an extension (.txt) then that will show up in the resulting file and analysis output. We recommend naming files to not have an extension so the output is easier to read.

"""

proc main() =
  # CLI configuration
  simple_parseopt.config: dash_dash_parameters . can_name_bare . no_implicit_bare . no_slash . manual_help

  # command line arguments
  let (options,got) = get_options_and_supplied:
    command: string          {. positional .}
    bimprefix: string        {. alias("b") .}
    annotprefix: string      {. alias("a") .}
    snpfiles: seq[string]    {. alias("s") .}
    outdir = "output"        {. alias("o") .}
    examples: bool           {. alias("e") .}
    help: bool               {. alias("h") .}

  if options.help:
    case options.command
    of "examples":
      echo &"{EXE} examples\n  Creates example input files you can use to test the tool.\n"
    of "make":
      echo &"{EXE} make --bimprefix PREFIX --snpfiles FILES... [--output DIRNAME]\n  Constructs new annotation files based on the SNP locations in the bim files, and the categories with optional annotations in the snpfiles."
    of "unmake":
      echo &"{EXE} unmake --annotprefix PREFIX [--output DIRNAME]\n  Deconstructs annotation files into their categories, and adds annotation values if they are continuous."
    else:
      echo HELP_MESSAGE
    quit(0)
  
  # generate examples if requested
  if options.command == "examples":
    generate_example_files()
    echo &"created example files: {baseline_filename}, {category_filename}"
    echo ""
    echo "run using the example files:"
    echo &"  annot -b example. {category_filename}"
    echo ""
    quit(0)

  if options.command == "make":
    # ensure both bimprefix and snpfiles arguments provided
    if not (got.bimprefix and got.snpfiles) or got.annotprefix:
      echo "Error: please specify a bim prefix with -b or --bimprefix"
      echo "       and snp category files with with -s or --snpfiles"
      quit(1)
    # run normal program flow
    let prefix = if options.bimprefix.endswith ".bim": options.bimprefix[0..^5] else: options.bimprefix
    let
      allsnpfiles = options.snpfiles.all_wildcards_expanded.sorted
      allbimfiles   = (prefix & "*.bim").all_wildcards_expanded.sorted
    for bimfile in allbimfiles:
      write_annot(bimfile, allsnpfiles, options.outdir)
    quit(0)

  if options.command == "unmake":
    if not got.annotprefix or got.bimprefix or got.snpfiles:
      echo "Error: please specify an annot prefix with -a or --annotprefix"
      quit(1)
    let prefix = if options.annotprefix.endswith ".bim": options.annotprefix[0..^5] else: options.annotprefix
    let allannotfiles = (prefix & "*.annot").all_wildcards_expanded.sorted
    allannotfiles.unmake(outdir=options.outdir)
    quit(0)

  # otherwise, unrecognized command
  echo &"Error: unrecognized command '{options.command}'"
  quit(1)

when isMainModule:  # Preserve ability to `import api`/call from Nim
  main()

