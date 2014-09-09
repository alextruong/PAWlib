#!/usr/bin/env python


from __future__ import division
import os
import sys
import glob

def read_data(filename):
	"""reads in data from file, creates list of list, returns list of list and headers"""

	data = []
	headers = []

	with open(filename, 'r') as f:
		for line in f:
			if line.startswith('#'):
				headers.append(line.strip())
				continue
			else:
				data.append(line.strip())
	
				
	data = [line.split('\t') for line in data]


	return data


def main():
	if len(sys.argv) != 2:
		print "Please call script: {0} <key file>".format(sys.argv[0])

	key_file = sys.argv[1]
	key_data = read_data(key_file)

	all_files = glob.glob('*')


	#outputs from processing script
	processed_files = glob.glob('*_snp.vcf')
	unprocessed_indels = glob.glob('*_indel.vcf')

	#outputs from gsearch and merge script
	gsearch_files = glob.glob('*_shared.vcf')
	merge_files = glob.glob('*_merged.vcf')
	match_coordinates = glob.glob('*_match_coordinates.txt')
	compare_merge_files = gsearch_files + merge_files + match_coordinates
	

	#output from VEP calling script
	annotated_files = glob.glob('*_annotated.vcf')
	annotated_html = glob.glob('*.html')
	VEP_calling_files = annotated_files + annotated_html

	#all outputs from binning script
	monoallelic_files = glob.glob('*_RNA-hom_WES-het.vcf')
	other_mismatch = glob.glob('*_RNA-het_WES-hom.vcf')
	matching_het = glob.glob('*_het-matching.vcf')
	matching_hom = glob.glob('*_hom-matching.vcf')
	matrix_counts = glob.glob('*_matrix_counts.txt')
	binning_files = monoallelic_files + other_mismatch + matching_het + matching_hom + matrix_counts

	#all outputs from variant compare script
	gene_lists = glob.glob('*_genelist*.txt')
	lineage_files = glob.glob('*lineage*.txt')
	unique_vcfs = glob.glob('*_unique_sorted.vcf')
	annotated_lists = glob.glob('*annotated*.txt')
	variant_comparison_files = gene_lists + lineage_files + unique_vcfs + annotated_lists

	#this creates a unique set of all of the family IDs
	family_ids = []
	for file_name in all_files:
		if file_name[0] not in family_ids:
			family_ids.append(file_name[0])
		else:
			continue

	#this grabs the original files, and the scripts
	scripts = [files for files in all_files if files.endswith('.py')]
	RNA_original_files = [i[2].strip() for i in key_data if i in all_files]
	WES_original_files = [i[3].strip() for i in key_data if i in all_files]
	original_files = RNA_original_files + WES_original_files + scripts

	print len(all_files)
	print len(unprocessed_indels), len(processed_files), len(gsearch_files), len(merge_files), len(match_coordinates), len(annotated_files), len(annotated_html), len(monoallelic_files), len(other_mismatch), len(matching_het), len(matching_hom), len(matrix_counts), len(gene_lists), len(lineage_files), len(unique_vcfs), len(annotated_lists)

	if len(all_files) != 2 + len(unprocessed_indels) + len(processed_files) + len(compare_merge_files) + len(VEP_calling_files) + len(binning_files) + len(variant_comparison_files) + len(original_files):
		print 'Script does not account for all file types'
		sys.exit(1)
	else:
		print 'Initalizing directories and moving files'

	os.system('mkdir original_files_and_scripts')

	directories = ['processed_files', 'compare_merge_files', 'binned_files', 'annotated_files', 'variant_comparison_files']

	#initializes all needed directories
	for family_id in family_ids:
		os.system('mkdir %s' % family_id)
		os.system('mkdir %s/snp' % family_id)
		for directory in directories:
			os.system('mkdir %s/snp/%s' % (family_id, directory))
			if directory == "variant_comparison_files":
				os.system('mkdir %s/snp/%s/gene_lists; mkdir %s/snp/%s/lineage_files; mkdir %s/snp/%s/unique_vcfs' % (family_id, directory, family_id, directory, family_id, directory))
		os.system('mkdir %s/indel' % family_id)

	#moves original data files to original folder
	for vcf in all_files:
		if vcf in original_files or (vcf == 'sampleKey_truncated.txt') or (vcf == "README.txt"):
			os.system('mv %s %s/' % (vcf, 'original_files_and_scripts'))

	#moves files into directories by what script produced them
		else:
			family_id = vcf.split('_')[0]

			if vcf in processed_files:
				os.system('mv %s %s/snp/%s' % (vcf, family_id, directories[0]))

			elif vcf in compare_merge_files:
				os.system('mv %s %s/snp/%s' % (vcf, family_id, directories[1]))

			elif vcf in binning_files:
				os.system('mv %s %s/snp/%s' % (vcf, family_id, directories[2]))

			elif vcf in VEP_calling_files:
				os.system('mv %s %s/snp/%s' % (vcf, family_id, directories[3]))

			elif vcf in variant_comparison_files:
				if vcf in gene_lists or annotated_lists:
					os.system('mv %s %s/snp/%s/%s' % (vcf, family_id, directories[4], 'gene_lists'))
				elif vcf in lineage_files:
					os.system('mv %s %s/snp/%s/%s' % (vcf, family_id, directories[4], 'lineage_files'))
				else:
					os.system('mv %s %s/snp/%s/%s' % (vcf, family_id, directories[4], 'unique_vcfs'))
			elif 'indel' in vcf:
				os.system('mv %s %s/indel' % (vcf, family_id))
		
			else:
				print vcf, 'File unaccounted for by script'
				sys.exit(1)

			

if __name__ == "__main__":
	main()
