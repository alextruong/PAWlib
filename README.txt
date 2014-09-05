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
2) compare_merge.py
3) bin_by_zygosity.py
4) VEPcaller.py
5) variant_analysis.py
6) final_sorting.py

Due to the variable nature of the input data, we request that all input files
be associated with a tab-delimited text file of the format

#FAMILYID ID  RNAFILENAME WESFILENAME GENDER
<familyid>  <p1|s1|mo|fa> <RNA-seq file>  <WES file>  <m|f>
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

To execute every other script except for final_sorting.py, just run them as-is, making sure they're in the same folder.

The output that these scripts generate (along with RNA/WES_processing) are all written to the same folder (the working directory), so again, it is strongly preferred that the working directory be otherwise empty when these scripts are executed. The file
To execute final_sorting.py, which will move all of these new files to proper directories, call it in a similar manner to RNA/WES_processing.py:

python final_sorting.py keyfile.txt

This will let the script make directories for the families and move the files accordingly.

FILE BREAKDOWN
-------------------------------------------------------------------------------
Input files (per quartet) (8):
4x RNA files
4x WES files

After RNA/WES_processing.py (16):
4x <familyid>_<p1|s1|mo|fa>-<m|f>_RNA_snp.vcf
4x <familyid>_<p1|s1|mo|fa>-<m|f>_WES_snp.vcf
4x <familyid>_<p1|s1|mo|fa>-<m|f>_RNA_indel.vcf
4x <familyid>_<p1|s1|mo|fa>-<m|f>_RNA_indel.vcf

After compare_merge.py (8):
4x <familyid>_<p1|s1|mo|fa>-<m|f>_incomplete_match_coordinates.txt
4x <familyid>_<p1|s1|mo|fa>-<m|f>_RNA_WES_snp_shared_sorted_merged.vcf

After bin_by_zygosity.py (17):
4x <familyid>_<p1|s1|mo|fa>-<m|f>_RNA_WES_snp_shared_sorted_merged_hom-matching.vcf
4x <familyid>_<p1|s1|mo|fa>-<m|f>_RNA_WES_snp_shared_sorted_merged_het-matching.vcf
4x <familyid>_<p1|s1|mo|fa>-<m|f>_RNA_WES_snp_shared_sorted_merged_RNA-het_WES-hom.vcf
4x <familyid>_<p1|s1|mo|fa>-<m|f>_RNA_WES_snp_shared_sorted_merged_RNA-hom_WES-het.vcf
1x <familyid>_RNA_WES_cross_pairwise_matrix_counts.txt

After VEPcaller.py (12):
4x <familyid>_<p1|s1|mo|fa>-<m|f>_RNA_snp_shared_sorted_merged_RNA-hom_WES-het_annotated.vcf
4x <familyid>_<p1|s1|mo|fa>-<m|f>_RNA_snp_shared_sorted_merged_RNA-hom_WES-het_annotated.vcf_summary.html
2x <familyid>_<p1|s1>-<m|f>_RNA_snp_shared_sorted_merged_annotated.vcf
2x <familyid>_<p1|s1>-<m|f>_RNA_snp_shared_sorted_merged_annotated.vcf_summary.html

After variant_analysis.py (18):
4x <familyid>_<p1|s1|mo|fa>-<m|f>_RNA_WES_snp_shared_sorted_merged_RNA-hom_WES-het_annotated_genelist.txt
4x <familyid>_<p1|s1|mo|fa>-<m|f>_RNA_WES_snp_shared_sorted_merged_RNA-hom_WES-het_annotated_monoallelic_gene_dictionary.txt
2x <familyid>_<p1|s1>-<m|f>_RNA_WES_snp_shared_sorted_merged_RNA-hom_WES-het_annotated_monoallelic_common_gene_dictionary.txt
2x <familyid>_<p1|s1>-<m|f>_RNA_WES_snp_shared_sorted_merged_RNA-hom_WES-het_annotated_monoallelic_common_sorted.vcf
2x <familyid>_<p1|s1>-<m|f>_RNA_WES_snp_shared_sorted_merged_RNA-hom_WES-het_annotated_monoallelic_unique_sorted.vcf
2x <familyid>_<p1|s1>-<m|f>_RNA_WES_snp_shared_sorted_merged_RNA-hom_WES-het_annotated_monoallelic_common_in_<s1|p1>_sorted.vcf
2x <familyid>_<p1|s1>-<m|f>_RNA_WES_snp_shared_sorted_merged_RNA-hom_WES-het_traced_WES.txt

General documents (9):
1x key file (*.txt)
7x scripts (*.py)
1x README (README.txt)

TOTAL: 88 files (79 files per family + 9 base files)

DESCRIPTION OF FILE BREAKDOWN (where necessary)
-----------------------------------------------------------------------------------
Input files (per quartet) (8):
4x RNA files
4x WES files

After RNA/WES_processing.py (16):
4x <familyid>_<p1|s1|mo|fa>-<m|f>_RNA_snp.vcf
4x <familyid>_<p1|s1|mo|fa>-<m|f>_WES_snp.vcf
4x <familyid>_<p1|s1|mo|fa>-<m|f>_RNA_indel.vcf
4x <familyid>_<p1|s1|mo|fa>-<m|f>_RNA_indel.vcf

These files are filtered by QUAL and GQ (dependent on zygosity), with default values of 20, 20 if heterozygous, 40 if homozygous.

After compare_merge.py (8):
4x <familyid>_<p1|s1|mo|fa>-<m|f>_incomplete_match_coordinates.txt

These files contain the locations of variants where a multiallelic variant shares a common alternate allele with a monoallelic variant (other tools such as bcftools and gsearch do not handle this appropriately, and it may or may not be desired).

4x <familyid>_<p1|s1|mo|fa>-<m|f>_RNA_WES_snp_shared_sorted_merged.vcf

These files contain the positions that are shared between an individual's RNA and WES files with both "extra" columns appended (using bcftools merge).

After bin_by_zygosity.py (17):
4x <familyid>_<p1|s1|mo|fa>-<m|f>_RNA_WES_snp_shared_sorted_merged_hom-matching.vcf

These files contain the variants that are homozygous in both WES and RNA files.

4x <familyid>_<p1|s1|mo|fa>-<m|f>_RNA_WES_snp_shared_sorted_merged_het-matching.vcf

These files contain the variants that are heterozygous in both WES and RNA files.

4x <familyid>_<p1|s1|mo|fa>-<m|f>_RNA_WES_snp_shared_sorted_merged_RNA-het_WES-hom.vcf

These files contain the variants that are heterozygous in RNA and homozygous in WES files.

4x <familyid>_<p1|s1|mo|fa>-<m|f>_RNA_WES_snp_shared_sorted_merged_RNA-hom_WES-het.vcf
1x <familyid>_RNA_WES_cross_pairwise_matrix_counts.txt

These files contain the variants that are homozygous in RNA and heterozygous in WES files.

After VEPcaller.py (12):
4x <familyid>_<p1|s1|mo|fa>-<m|f>_RNA_snp_shared_sorted_merged_RNA-hom_WES-het_annotated.vcf

These files are the annotated versions of the corresponding files above.

4x <familyid>_<p1|s1|mo|fa>-<m|f>_RNA_snp_shared_sorted_merged_RNA-hom_WES-het_annotated.vcf_summary.html

These files are the html summary files for the annotated files above.

2x <familyid>_<p1|s1>-<m|f>_RNA_snp_shared_sorted_merged_annotated.vcf

These files are the annotated versions of the corresponding files above.

2x <familyid>_<p1|s1>-<m|f>_RNA_snp_shared_sorted_merged_annotated.vcf_summary.html

These files are the html summary files for the annotated files above.

After variant_analysis.py (18):
4x <familyid>_<p1|s1|mo|fa>-<m|f>_RNA_WES_snp_shared_sorted_merged_RNA-hom_WES-het_annotated_genelist.txt

These files are the lists of genes derived from the corrresponding annotated files above.

4x <familyid>_<p1|s1|mo|fa>-<m|f>_RNA_WES_snp_shared_sorted_merged_RNA-hom_WES-het_annotated_monoallelic_gene_dictionary.txt

These files are the lists of genes derived from the corresponding files above, with associated representative variant positions.

2x <familyid>_<p1|s1>-<m|f>_RNA_WES_snp_shared_sorted_merged_RNA-hom_WES-het_annotated_monoallelic_common_gene_dictionary.txt

These files contain the common genes and positions derived from the comparison between the RNA-hom_WES-het file of the query (p1|s1) compared to the _merged file of the reference(s1|p1).

2x <familyid>_<p1|s1>-<m|f>_RNA_WES_snp_shared_sorted_merged_RNA-hom_WES-het_annotated_monoallelic_common_sorted.vcf

These files contain the common variants derived from the comparison between the RNA-hom_WES-het file of the query (p1|s1) compared to the _merged file of the reference(s1|p1).

2x <familyid>_<p1|s1>-<m|f>_RNA_WES_snp_shared_sorted_merged_RNA-hom_WES-het_annotated_monoallelic_unique_sorted.vcf

These files contain the unique variants (query-side) derived from the comparison between the RNA-hom_WES-het file of the query (p1|s1) compared to the _merged file of the reference(s1|p1).

2x <familyid>_<p1|s1>-<m|f>_RNA_WES_snp_shared_sorted_merged_RNA-hom_WES-het_annotated_monoallelic_common_in_<s1|p1>_sorted.vcf

These files contain the common variants (reference-side) derived from the comparison between the RNA-hom_WES-het file of the query (p1|s1) compared to the _merged file of the reference(s1|p1).

2x <familyid>_<p1|s1>-<m|f>_RNA_WES_snp_shared_sorted_merged_RNA-hom_WES-het_traced_WES.txt

These files contain the zygosity information from the parents given an individual's RNA-hom_WES-het file.







The variant_compare script generates a gene list of monoallelic genes present for proband, sibling, father and mother. These lists can be put into a venn diagram generating tool to generate a venn diagram to see the intersection of the two gene lists.
