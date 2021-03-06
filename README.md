# ldsc-annot-tool
CLI tool to convert GWAS rsid-categories to the annot files that [ld-score-estimation (Alkes)](https://github.com/bulik/ldsc) requires.

[![Build status](https://ci.appveyor.com/api/projects/status/ej5xh7odg28rcwrt?svg=true)](https://ci.appveyor.com/project/JorySchossau/ldsc-annot-tool)

### Downloads
Latest releases listed here.
* [Linux](https://github.com/JorySchossau/ldsc-annot-tool/releases/latest/download/lin_annot.zip) (glibc 2.23+)
* [Mac](https://github.com/JorySchossau/ldsc-annot-tool/releases/latest/download/osx_annot.zip) (10.13+)
* [Windows](https://github.com/JorySchossau/ldsc-annot-tool/releases/latest/download/win_annot.zip)

[Older Releases](https://github.com/JorySchossau/ldsc-annot-tool/releases)

### About
The [LDSC tool](https://github.com/bulik/ldsc) includes a `makeannot.py` script, but sometimes you simply have a bunch of RSIDs and need to make annot files. This ldsc-annot-tool does just that.

```
Converts a collection of rsid categories into formatted .annot files for ldsc regression.

Usage
  annot examples
  annot make -b|--bimprefix PATH -s|--snpfiles FILES... [-o|--outdir OUTDIR]
  annot unmake -a|--annotprefix PATH [-o|--outdir OUTDIR]

Example
  ./annot -b test_baseline/subset. -s test_annots/FetalDHS_Trynka test_annots/H3K27ac_PGC2

Commands
  examples   create example input files
  make       construct annot files from category annotations
  unmake     deconstruct annot files and create category annotations

  examples
             create example input files
  [no options]

  make
             construct annot files from category annotations
  -b --bimprefix PATH    nonunique prefix path to all bim files
  -s --snpfiles FILES... all snp category files
  -o --outdir OUTDIR     output directory, optional (default 'output')

  unmake
             deconstruct annot files and create category annotations
  -a --annotprefix PATH    nonunique prefix path to all annot files
  -o --outdir OUTDIR     output directory, optional (default 'output')
  -h --help              this help message

Note that the snp file names will be used as category headers (names) in the resulting file. So if the files have an extension (.txt) then that will show up in the resulting file and analysis output. We recommend naming files to not have an extension so the output is easier to read.
```

### How to Use
Given a baseline (bim,bed,fam files) that are prepared for the [LDSC tool](https://github.com/bulik/ldsc), you can add arbitrary lists of RSIDs as new annotations.

For example, suppose you have an example baseline in `baseline/` as a set of files:

```sh
baseline/1kg_eur.1.bed
baseline/1kg_eur.1.bim
baseline/1kg_eur.1.fam
baseline/1kg_eur.2.bed
baseline/1kg_eur.2.bim
baseline/1kg_eur.2.fam
...
baseline/1kg_eur.22.fam
```

and suppose you have two new groups of RSIDs, each listed in their own files, like so:

(In the directory `categories` a file named `Exons`): `categories/Exons`
```sh
rs175643
rs15574
rs275644
...
```

and

(In the directory `categories` a file named `Introns`): `categories/Introns`
```sh
rs266758
rs211515
rs275230
...
```

then you can prepare the annotations for the Alkes ldsc tool, like this:

```sh
./annot make --bimprefix baseline/1kg_eur. --snpfiles categories/Introns categories/Exons
```

or

```sh
./annot make --bimprefix baseline/1kg_eur. --snpfiles categories/*
```

or with optional output destination

```sh
./annot make --outdir myoutput --bimprefix baseline/1kg_eur. --snpfiles categories/*
```

resulting in the new annotation files created in a directory `newAnnots` ready for ldsc regression, with each position annotated with a `0` or `1` for each category given, if the RSID is present (in this case, the `Exons` and `Introns` categories).

All command line options are repositionable (`--outdir myoutput` may be last, for example)

