PAWlib: Pedigree Analysis Workflow (rough)
-------------------------------------------------------------------------------
INTRODUCTION
-------------------------------------------------------------------------------
This is not so much a library as it is a compilation of python scripts that are designed to do work in sequence. The sum purpose of these scripts is to facilitate simple pedigree analysis for a quartet: 

mother, father, proband, sibling

Certain methods in these scripts (notably variant_compare) will break if a trio or duet is passed.

-------------------------------------------------------------------------------
OVERVIEW OF SCRIPTS
-------------------------------------------------------------------------------
There are 6 scripts that do all the work in this pipeline. In order, they are:

1) VCF_processing.py - Input: All vcf files (RNA and WES), Output: Filtered VCF files by gq/quality, and snp/indel
2) compare_merge.py - Input: Output SNP files from step 1, Output: Merged files
3) bin_by_zygosity.py - Input: Merged files from step 2, Output: Counts of variants by zygosity (.txt), VCFs for 2x2 table 
4) VEPcaller.py - Input: All RNA-homozygous/WES-heterozygous VCFs, and proband/sibling full merges, Output: Annotated input files
5) variant_analysis.py - Input: Annotated files from step 4, original WES files from parents, RNA-homozygous/WES-heterozygous files from step 3, and proband/sibling full merges from step 2, Output: 
6) final_sorting.py - Input: All files in directory, Output: Directories by family id, SNPS/INDELS, and script outputs

Due to the variable nature of the input data, we request that all input files
be associated with a tab-delimited text file of the format

#FAMILYID ID  RNAFILENAME WESFILENAME GENDER
<familyid>  <p1|s1|mo|fa> <RNA-seq file>  <WES file>  <m|f>
.
.
.

where p1 = proband, s1 = sibling, mo = mother, and fa = father.


-------------------------------------------------------------------------------
DEPENDENCIES
-------------------------------------------------------------------------------
In order for these scripts to function properly, you need to have bcftools (http://samtools.github.io/bcftools/), vcftools (http://vcftools.sourceforge.net/), and the standalone perl version of VEP (http://www.ensembl.org/info/docs/tools/vep/index.html) installed on the machine and callable from commandline.

These scripts were tested on VCF v.4.1 files; see below disclaimer for further caveats.


-------------------------------------------------------------------------------
EXECUTION OF SCRIPTS
-------------------------------------------------------------------------------
***
DISCLAIMER: These scripts, on occasion, make calls to the terminal to execute
commands. These scripts were developed on Unix/Linux machines; there is no guarantee of success if used in a non Unix/Linux environment. 

Windows Powershell should be able to execute these scripts in theory, but extensive testing has not been performed. There are currently no plans to guarantee Windows compatibility.

These scripts were developed in a Python 2.x environment; no testing for compatibility with Python 3.x was performed, and thus is not supported.
***

IN ORDER, EXECUTE AS FOLLOWS (or similarly depending on environment):
1) VCF_processing.py

python VCF_processing.py <keyfile.txt>

2) compare_merge.py

python compare_merge.py

3) bin_by_zygosity

python bin_by_zygosity.py

4) VEPcaller.py

python VEPcaller.py

5) variant_comparison.py

python variant_comparison.py

6) final_sorting.py

python final_sorting.py <keyfile.txt>


-------------------------------------------------------------------------------
FILE BREAKDOWN
-------------------------------------------------------------------------------
For a single set of input files (per quartet) (8):
4x RNA files
4x WES files

-----

After VCF_processing.py (16):
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

After VEPcaller.py (16):
4x <familyid>_<p1|s1|mo|fa>-<m|f>_RNA_snp_shared_sorted_merged_RNA-hom_WES-het_annotated.vcf
4x <familyid>_<p1|s1|mo|fa>-<m|f>_RNA_snp_shared_sorted_merged_RNA-hom_WES-het_annotated.vcf_summary.html
4x <familyid>_<p1|s1|mo|fa>-<m|f>_RNA_snp_shared_sorted_merged_RNA-het_WES-hom_annotated.vcf
4x <familyid>_<p1|s1|mo|fa>-<m|f>_RNA_snp_shared_sorted_merged_RNA-het_WES-hom_annotated.vcf_summary.html

After variant_comparison.py (9):
1x <familyid>_common_p1_s1_monoallelic_genelist_variants.txt
2x <familyid>_<p1|s1>-<m|f>_RNA_WES_snp_shared_sorted_merged_RNA-hom_WES-het_annotated_unique_monoallelic_genelist_variants.txt
2x <familyid>_<mo|fa>-<m|f>_RNA_WES_snp_shared_sorted_merged_RNA-hom_WES-het_annotated_full_monoallelic_genelist_variants.txt
2x <familyid>_<p1|s1>-<m|f>_RNA_WES_snp_shared_sorted_merged_RNA-hom_WES-het_unique_lineage_parent_WES.txt
2x <familyid>_<p1|s1>-<m|f>_RNA_WES_snp_shared_sorted_merged_RNA-hom_WES-het_unique_sorted.vcf

TOTAL: 66 output files per family


-------------------------------------------------------------------------------
SCRIPT FUNCTIONALITY
-------------------------------------------------------------------------------
VCF_processing - 
Takes input files, filters by QUAL/GQ, and separates SNPs from indels.

compare_merge.py - 
Takes RNA SNP and WES SNP files and grabs all the shared variants, then consolidates the information into a single "merge" file (one variant, two sets of information: RNA and WES, e.g. zygosity, scores, etc.) per individual. Matches at multiallelic sites are noted in a separate file.

bin_by_zygosity.py - 
Takes the merge files and generates a set of counts based on zygosity combinations (e.g. RNA heterozygous/WES homozygous, etc.). Files containing the corresponding variants are also genrated.

VEPcaller.py - executes a command-line version of VEP to annotate all fles from the previous step that are RNA homozygous/WES heterozygous (monoallelic), as well as both the proband and sibling merge files from the step before. Generates the annotated files and corresponding HTML summary files.

variant_analysis.py - does a number of things: 
1) compare proband monoallelic variant files to sibling merge files, and vice versa. Record common variants and information from both sides; output in files appropriately named. 
1.5) generate gene lists based on step 1.
2) Take monoallelic child variants (proband|sibling) and trace them in mother and father WES files. Note zygosities to infer inheritance of variants.
3) take parental monoallelic variants and generate lists of interesting genes (for external analysis).


---------------------------------------------------------------------------------
CLOSING REMARKS
---------------------------------------------------------------------------------
This was the culmination of a summer's worth of work at Boston Children's Hospital in an internship under Dr. Sek Won Kong.
Please contact Joey Orofino, Alex Truong, or Deep Shah (in that order) if there are any questions regarding the code.

Joey Orofino: jorofino@bu.edu
Alex Truong: atruong@bu.edu
Deep Shah: deepsh@bu.edu
