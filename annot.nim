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

# could make this modifiable by CLI
var verbosity:int = 1

# 1 intset is 1 category file listing rsids
# the string is the filename sans extension
type Categories = OrderedTable[string,IntSet]

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
  var rsid = 0
  for line in file.lines:
    if not line.scanf("rs$i$.",rsid):
      stderr.writeLine("")
      stderr.writeLine(&"Error: {file} is not a category file containing rsids")
      stderr.writeLine("use --examples to generate example input files")
      quit(1)
    break

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

# convert all files passed to
# table of intsets of rsids, grouped by filename
proc to_categories(files:openArray[string]): Categories =
  for file in files:
    if verbosity >= 1:
      echo &"reading {file}"
    result[file.splitFile.name] = initIntSet()
    var category {.byaddr.} = result[file.splitFile.name]
    var rsid:int
    check_category_format file
    for line in file.lines:
      if line.scanf("rs$i",rsid): category.incl rsid

# the main work proc
# builds off of data in bim file
# adding extra columns, 1 for each category
# putting a 1 or 0 if rsid is in category
proc write_annot(bimfile:string, categories:Categories, outdir:string) =
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
  for name,category in categories.pairs:
    outfile.write("\t" & name)
  outfile.write("\n")

  # write rest of contents
  for line in bimfile.lines:
    if line.scanf("$i\trs$i\t$f\t$i",chr,rsid,cm,bp):
      # write first part of row
      outfile.write &"{chr}\t{bp}\t{rsid}\t{cm}\t1"
      # write rest of row (all categories)
      for name,category in categories:
        if category.contains rsid: outfile.write "\t1"
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

proc showHelp =
  echo """

Converts a collection of rsid categories into formatted .annot files for ldsc regression.

Usage
  annot -b|--bimprefix PATH -s|--snpfiles FILES... [-o|--output OUTDIR]

Example
  ./annot -b test_baseline/subset. -s test_annots/FetalDHS_Trynka test_annots/H3K27ac_PGC2

  -b --bimprefix PATH    nonunique prefix path to all bim files
  -s --snpfiles FILES... all snp category files
  -o --output OUTDIR     output directory (default 'output')
  -h --help              this help message

Note that the snp file names will be used as category headers (names) in the resulting file. So if the files have an extension (.txt) then that will show up in the resulting file and analysis output. We recommend naming files to not have an extension so the output is easier to read.

"""

proc errorCheckCommandLineArgs(args:seq[string]) =
  # show help if requested
  if ["-h","--help"] in args or args.len == 0:
    showHelp()
    quit(0)

  # generate examples if requested
  if "--examples" in args:
    generate_example_files()
    echo &"created example files: {baseline_filename}, {category_filename}"
    echo ""
    echo "run using the example files:"
    echo &"  annot -b example. {category_filename}"
    echo ""
    quit(0)

  if ["-b","--bimprefix"] notin args:
    echo "Error: please specify a bim prefix with -b or --bimprefix"
    echo "       ex: bims/baseline.1.bim through bims/baseline.22.bim use:"
    echo "       annot -b bims/baseline."
    quit(1)

  if ["-s","--snpfiles"] notin args:
    echo "Error: please specify snp category files with -s or --snpfiles"
    echo "       ex: snps/FetalDHS_Trynka and snps/H3K27ac_PGC2 use:"
    echo "       annot -s snps/*"
    echo "       or"
    echo "       annot -s snps/FetalDHS_Trynka snps/H3K27ac_PGC2"
    quit(1)

  if args.len < 4:
    echo "Error: both --bimprefix (-b) and --snpfiles (-s) need arguments."
    echo ""
    showHelp()
    quit(1)

proc main() =
  var
    outdir:string = "output" # default output directory, user can override
    bimprefix:string
    snpfiles:seq[string]
    makeExamples:bool
    args = commandLineParams()
    loc = -1

  # remove possible optional output option and its argument
  if ["-o","--output"] in args:
    loc = -1
    loc = max(loc, args.find "-o")
    loc = max(loc, args.find "--output")
    if loc != -1:
      if args[loc+1].startsWith("-"):
        echo "Error: --output (-o) needs an argument."
        echo ""
        showHelp()
        quit(1)
      outdir = args[loc+1]
      # delete "--output" and following argument
      args.delete loc
      args.delete loc
  
  errorCheckCommandLineArgs args

  # find, load, and remove bimprefix and its argument
  loc = -1
  loc = max(loc, args.find "-b")
  loc = max(loc, args.find "--bimprefix")
  bimprefix = args[loc+1]

  args.delete loc
  args.delete loc

  # find and load snpfiles arguments
  loc = -1
  loc = max(loc, args.find "-s")
  loc = max(loc, args.find "--snpfiles")
  snpfiles = args[loc+1 .. ^1]

  # run normal program flow
  let
    allsnpfiles = snpfiles.all_wildcards_expanded.sorted
    allbimfiles   = (bimprefix & "*.bim").all_wildcards_expanded.sorted
    categories = allsnpfiles.to_categories

  for bimfile in allbimfiles:
    write_annot(bimfile, categories, outdir)

when isMainModule:  # Preserve ability to `import api`/call from Nim
  main()

