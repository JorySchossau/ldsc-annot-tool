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
import parsecsv

# could make this modifiable by CLI
var verbosity:int = 1

# 1 intset is 1 category file listing rsids
# the string is the filename sans extension
type
  Categories = OrderedTable[string,IntSet]
  CategoriesWithAnnotations = Table[string,Table[int,float]]

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

# convert all files passed to
# table of intsets of rsids, grouped by filename
proc to_categories(files:openArray[string]): (Categories,CategoriesWithAnnotations) =
  var categories:Categories
  var categoriesWithAnnots:CategoriesWithAnnotations
  for file in files:
    if verbosity >= 1:
      echo &"reading {file}"
    # init the category
    categories[file.splitFile.name] = initIntSet()
    var category {.byaddr.} = categories[file.splitFile.name]
    var rsid:int
    var annotation:float
    check_category_format file
    if file.has_annotations:
      # init the annotation table
      categoriesWithAnnots[file.splitFile.name] = initTable[int,float]()
      var annottable {.byaddr.} = categoriesWithAnnots[file.splitFile.name]
      # scrape both rsid and annotation
      for line in file.lines:
        if line.scanf("rs$i$s$f$.",rsid,annotation):
          category.incl rsid
          annottable[rsid] = annotation
    else:
      # scrape only the rsids
      for line in file.lines:
        if line.scanf("rs$i$.",rsid):
          category.incl rsid
  result = (categories,categoriesWithAnnots)

# the main work proc
# builds off of data in bim file
# adding extra columns, 1 for each category
# putting a 1 or 0 if rsid is in category
proc write_annot(bimfile:string, cats:Categories, catswAnnots:CategoriesWithAnnotations, outdir:string) =
  check_bim_format bimfile
  var
    chr:int
    rsid:int
    cm:float
    bp:int

  # prepare output file
  create_dir outdir
  let outfilename = bimfile.split_file.name & ".annot"
  var outfile = newFileStream(outdir / outfilename, fmWrite)

  if verbosity >= 1:
    echo &"reading {bimfile}, creating " & outdir / outfilename

  # write header
  outfile.write("CHR\tBP\tSNP\tCM\tbase")
  for name,category in cats.pairs:
    outfile.write("\t" & name)
  outfile.write("\n")

  # write rest of contents
  for line in bimfile.lines:
    if line.scanf("$i\trs$i\t$f\t$i",chr,rsid,cm,bp):
      # write first part of row
      outfile.write &"{chr}\t{bp}\t{rsid}\t{cm}\t1"
      # write rest of row (all cats)
      for name,category in cats:
        # if this category has no annotations
        if not catswAnnots.hasKey name:
          if category.contains rsid: outfile.write "\t1"
          else:                      outfile.write "\t0"
        # if this category has annotations
        else:
          let catwannots = catswAnnots[name]
          if category.contains rsid: outfile.write "\t" & $catwannots[rsid]
          else:                      outfile.write "\t0"
      # finish row EOL
      outfile.write "\n"

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

proc has_continuous_values(col:var Table[int,float]):bool =
  for value in col.mvalues:
    if value != 0 and value != 1:
      return true
  return false

const HELP_MESSAGE = """

Converts a collection of rsid categories into formatted .annot files for ldsc regression.

Usage
  annot examples
  annot make -b|--bimprefix PATH -s|--snpfiles FILES... [-o|--output OUTDIR]
  annot unmake -b|--bimprefix PATH [-o|--output OUTDIR]

Example
  ./annot -b test_baseline/subset. -s test_annots/FetalDHS_Trynka test_annots/H3K27ac_PGC2

Commands
  examples   create example input files
  make       construct annot files from category annotations
  unmake     deconstruct annot files and create category annotations

  -b --bimprefix PATH    nonunique prefix path to all bim files
  -s --snpfiles FILES... all snp category files
  -o --output OUTDIR     output directory (default 'output')
  -h --help              this help message

Note that the snp file names will be used as category headers (names) in the resulting file. So if the files have an extension (.txt) then that will show up in the resulting file and analysis output. We recommend naming files to not have an extension so the output is easier to read.

"""
const EXE = "annot"

proc main() =
  # CLI configuration
  simple_parseopt.config: dash_dash_parameters . can_name_bare . no_implicit_bare . manual_help
  simple_parseopt.config: can_name_bare
  #help_text(customHelpHeader, customHelpFooter)

  # command line arguments
  let (options,got) = get_options_and_supplied:
    command: string          {. positional, info("Run mode: examples, make, unmake") .}
    bimprefix: string        {. alias("b"), info("Nonunique path and filename shared by all bim files") .}
    snpfiles: seq[string]    {. alias("s"), info("List of snp category files") .}
    outdir = "output"        {. alias("o"), info("Output directory (default 'output')") .}
    examples: bool           {. alias("e"), info("Create example input files") .}
    help: bool               {. alias("h"), info("Show this help") .}

  if options.help:
    case options.command
    of "examples":
      echo &"{EXE} examples\n  Creates example input files you can use to test the tool.\n"
    of "make":
      echo &"{EXE} make --bimprefix PREFIX --snpfiles FILES... [--output DIRNAME]\n  Constructs new annotation files based on the SNP locations in the bim files, and the categories with optional annotations in the snpfiles."
    of "unmake":
      echo &"{EXE} unmake --bimprefix PREFIX [--output DIRNAME]\n  Deconstructs annotation files into their categories, and adds annotation values if they are continuous."
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
    if not (got.bimprefix and got.snpfiles):
      echo "Error: please specify a bim prefix with -b or --bimprefix"
      echo "       and snp category files with with -s or --snpfiles"
      quit(1)
    # run normal program flow
    let
      allsnpfiles = options.snpfiles.all_wildcards_expanded.sorted
      allbimfiles   = (options.bimprefix & "*.bim").all_wildcards_expanded.sorted
      (categories,categoriesWithAnnots) = allsnpfiles.to_categories
    for bimfile in allbimfiles:
      write_annot(bimfile, categories, categoriesWithAnnots, options.outdir)
    quit(0)

  if options.command == "unmake":
    if not got.bimprefix:
      echo "Error: please specify a bim prefix with -b or --bimprefix"
      quit(1)
    let
      allannotfiles   = (options.bimprefix & "*.annot").all_wildcards_expanded.sorted
    # indexed [5+colid,[rsid,annot]]
    # 5+ so it's same as original annot file column index
    var colstable = initTable[int,Table[int,float]]()
    var colnames:seq[string]
    var file: CsvParser
    let SNPCOL = 2
    let FIRSTUSERCOL = 5
    for annotfile in allannotfiles:
      if verbosity >= 1:
        echo &"reading {annotfile}"
      file.open(annotfile, separator='\t')
      file.readHeaderRow()
      # init columns if necessary
      for column in file.headers[FIRSTUSERCOL .. ^1]:
        if column notin colnames:
          colnames.add column
          let adjustedColId = colnames.len-1+FIRSTUSERCOL
          colstable[adjustedColId] = initTable[int,float]()
      while file.readRow():
        for coln in FIRSTUSERCOL ..< file.headers.len:
          colstable[coln][file.row[SNPCOL].parseInt] = file.row[coln].parseFloat
    # prepare output file
    create_dir options.outdir
    for coln,column in colnames:
      let outfilename = column
      var outfile = newFileStream(options.outdir / outfilename, fmWrite)
      if verbosity >= 1:
        echo &"creating " & options.outdir / outfilename
      if colstable[coln+FIRSTUSERCOL].has_continuous_values:
        for rsid,value in colstable[coln+FIRSTUSERCOL].mpairs:
          outfile.write "rs" & $rsid & "\t" & $value & "\n"
      else:
        for rsid in colstable[coln+FIRSTUSERCOL].keys:
          outfile.write "rs" & $rsid & "\n"
    quit(0)

  # otherwise, unrecognized command
  echo "Error: unrecognized command"
  echo HELP_MESSAGE
  quit(1)

when isMainModule:  # Preserve ability to `import api`/call from Nim
  main()

