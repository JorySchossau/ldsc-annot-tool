import strscans
import glob # dependency
import cligen # dependency
import streams
import seqUtils
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
    result[file] = initIntSet()
    var category {.byaddr.} = result[file]
    var rsid:int
    check_category_format file
    for line in file.lines:
      if line.scanf("rs$i",rsid): category.incl rsid

# the main work proc
# builds off of data in bim file
# adding extra columns, 1 for each category
# putting a 1 or 0 if rsid is in category
proc write_annot(bimfile:string, categories:Categories) =
  check_bim_format bimfile
  var
    chr:int
    rsid:int
    cm:float
    bp:int

  # prepare output file
  create_dir "newAnnots"
  let outfilename = bimfile.split_file.name & ".annot"
  var outfile = newFileStream("newAnnots" / outfilename, fmWrite)

  if verbosity >= 1:
    echo &"reading {bimfile}, creating " & "newAnnots" / outfilename

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

proc main(bimprefix: string = "REQUIRED", annotfiles: seq[string], examples:bool = false) =

  # generate examples and quit, if requested
  if examples:
    generate_example_files()
    echo &"saved example files: {baseline_filename}, {category_filename}"
    echo ""
    echo "run using the example files:"
    echo &"  annot -b example. {category_filename}"
    echo ""
    quit(0)

  if bimprefix == "REQUIRED":
    echo "Error: please specify a bim prefix with -b or --bimprefix"
    echo "       ex: bims/baseline.1.bim through bims/baseline.22.bim use:"
    echo "       annot -b bims/baseline."
    quit(1)

  # otherwise, run normal program flow
  let
    allannotfiles = annotfiles.all_wildcards_expanded.sorted
    allbimfiles   = (bimprefix & "*.bim").all_wildcards_expanded.sorted
    categories = allannotfiles.to_categories
  for bimfile in allbimfiles:
    write_annot(bimfile, categories)

when isMainModule:  # Preserve ability to `import api`/call from Nim
  dispatch(main, help={"examples":"create and show example usage"})

