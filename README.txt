PAWlib
======

Pedigree Analysis Workflow (rough)

INTRODUCTION
-------------------------------------------------------------------------------
This is not so much a library as it is a compilation of python scripts that are designed to do work in sequence. The sum purpose of these scripts is to facilitate simple pedigree analysis for a quartet: 

mother, father, proband, sibling

Certain methods in these scripts (notably variant_compare) will break if a trio or duet is passed.


OVERVIEW OF SCRIPTS
-------------------------------------------------------------------------------
There are 6 (7) scripts that do all the work in this pipeline. In order, they are:

1) RNA_processing.py (WES_processing.py)
2) batch_compare_and_merge.py
3) batch_filter_zygosity_common_pairwise_split.py
4) pythonrunVEP.py
5) variant_compare.py
6) sorting.py

Due to the variable nature of the input data, we request that all input files
be associated with a tab-delimited text file of the format

#FAMILYID ID  RNAFILENAME WESFILENAME GENDER
#<familyid>  <p1|s1|mo|fa> <RNA-seq file>  <WES file>  <m|f>
.
.
.

where p1 = proband, s1 = sibling, mo = mother, and fa = father.

EXECUTION OF SCRIPTS
-------------------------------------------------------------------------------
The first 5(6) scripts need to be in the same directory as the input files to
run. These scripts are automatically batched; if for some reason it is required
to run the script on a specific file(s), quarantine those files and move/copy the scripts to that new directory to run in isolation.

Regardless of where the working directory is located, it is optimal that it be
empty, or devoid of VCF files that are not intended to be processed. Ideally,
the directory should be empty.

***

DISCLAIMER: These scripts, on occasion, make calls to the terminal to execute
commands. These scripts were developed on Unix/Linux machines; there is no guarantee of success if used in a non Unix/Linux environment. 

Windows Powershell should be able to execute these scripts in theory, but extensive testing has not been performed. There are currently no plans to guarantee Windows compatibility.

These scripts were developed in a Python 2.x environment; no testing for compatibility with Python 3.x was performed, and thus is not supported.

***

To execute RNA_processing.py (or WES_processing.py), specify the key file as an argument:

python RNA_processing.py keyfile.txt

or something to that effect. The script will inform of erroneous execution with similar instructions.

To execute every other script except for sorting.py, just run them as-is, making sure they're in the same folder.

The output that these scripts generate (along with RNA/WES_processing) are all written to the same folder (the working directory), so again, it is strongly preferred that the working directory be otherwise empty when these scripts are executed. For 2 quartets, the final output of variant_compare puts the number of total files in the directory at 224, including original files and script files.

To execute sorting.py, which will move all of these new files to proper directories, call it in a similar manner to RNA/WES_processing.py:

python sorting.py keyfile.txt

This will let the script make directories for the families and move the files accordingly.

The variant_compare script generates a gene list of monoallelic genes present for proband, sibling, father and mother. These lists can be put into a venn diagram generating tool like venny (http://bioinfogp.cnb.csic.es/tools/venny/) to generate a venn diagram to see the intersection of the two gene lists.
