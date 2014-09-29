#!/usr/bin/env python

from __future__ import division

import os
import sys
import itertools
import glob

def read_data(filename):
        data = []
        headers = []
        
        with open(filename, 'r') as f:
                for line in f:
                        if line.startswith('#'):
                                headers.append(line.strip())
                        else:
                                data.append(line.strip())
        
        data = [line.split('\t') for line in data]
        
        
        return data, headers
        
        
def parse_combinations(data, file_name):

        #initializes counts
        matching_het_count = 0
        matching_hom_count = 0
        rna_hom_wes_het_count = 0
        rna_het_wes_hom_count = 0
        
        name = file_name
        file_name = []
        master_file = []
        
        #initializes arrays
        monoallelic_variants = []
        other_mismatch = []
        matching_het = []
        matching_hom = []

        #for each element in the data, check what pattern of zygosity each variant displays and bin it into the appropriate group
        for i in xrange(len(data)):
                wes_genotype = data[i][-1].split(':')[0]
                rna_genotype = data[i][-2].split(':')[0]
        
                if (wes_genotype == rna_genotype) and (wes_genotype == '0/1' and rna_genotype == '0/1'):        
                        matching_het_count += 1
                        matching_het.append(data[i])
                
                elif (wes_genotype == rna_genotype) and (wes_genotype == '1/1' and rna_genotype == '1/1'):
                        matching_hom_count += 1
                        matching_hom.append(data[i])
                
                elif (wes_genotype != rna_genotype) and (wes_genotype == '1/1' and rna_genotype == '0/1'):
                        rna_het_wes_hom_count += 1      
                        other_mismatch.append(data[i])
                
                elif (wes_genotype != rna_genotype) and (wes_genotype == '0/1' and rna_genotype == '1/1'):
                        rna_hom_wes_het_count += 1
                        monoallelic_variants.append(data[i])
                #this assumes 0/1, 1/1 are the only combinations, remove the next three lines if the data has diff. zygosities
                else:
                        print 'Parsing failed, columns in wrong position'
                        sys.exit(1)
        
        file_name.append(name)
        file_name.append(matching_het_count)
        file_name.append(matching_hom_count)
        file_name.append(rna_het_wes_hom_count)
        file_name.append(rna_hom_wes_het_count)
        print 'Counts written'
        
        
        return file_name, monoallelic_variants, other_mismatch, matching_het, matching_hom
        
def write_to_file(master_file, family_number):
        
        file_name = family_number
        print 'Writing File...'
        
        #writes the counts from each individual
        with open('%s_RNA_WES_cross_pairwise_matrix_counts.txt' % file_name, 'w') as w:
                w.write('\t')
                for x in xrange(len(master_file)):
                        w.write('%s\t' % master_file[x][0])
                        
                w.write('\n')
                w.write('heterozygous\t')
                for y in xrange(len(master_file)):
                        w.write('%s\t' % master_file[y][1])
                
                w.write('\n')
                w.write('homozygous\t')
                for z in xrange(len(master_file)):
                        w.write('%s\t' % master_file[z][2])     
                        
                w.write('\n')
                w.write('RNA-heterozygous, WES-homozygous\t')
                for a in xrange(len(master_file)):
                        w.write('%s\t' % master_file[a][3])
                        
                w.write('\n')
                w.write('RNA-homozygous, WES-heterozygous\t')
                for b in xrange(len(master_file)):
                        w.write('%s\t' % master_file[b][4])
        
        print 'Done!'


def write_vcf_files(monoallelic_variants, other_mismatch, matching_het, matching_hom, file_name, headers):
        
        
        monoallelic_file = file_name[0:-4] + '_RNA-hom_WES-het.vcf'
        
        print 'Writing monoallelic file...'
        
        with open(monoallelic_file, 'w') as w:
                for line in headers:
                        w.write(line + '\n')
                
                for line in monoallelic_variants:
                        for index, row in enumerate(line):
                                if index != len(line) - 1:
                                        w.write('%s\t' % row)
                                else:
                                        w.write('%s\n' % row)
        
        other_mismatch_file = file_name[0:-4] + '_RNA-het_WES-hom.vcf'
        
        print 'Writing other mismatch file...'
        
        with open(other_mismatch_file, 'w') as w:
                for line in headers:
                        w.write(line + '\n')
                
                for line in other_mismatch:
                        for index, row in enumerate(line):
                                if index != len(line) - 1:
                                        w.write('%s\t' % row)
                                else:
                                        w.write('%s\n' % row)                
                        
        matching_het_file = file_name[0:-4] + '_het-matching.vcf'
        
        print 'Writing matching heterozygous file...'
        
        with open(matching_het_file, 'w') as w:
                for line in headers:
                        w.write(line + '\n')
                
                for line in matching_het:
                        for index, row in enumerate(line):
                                if index != len(line) - 1:
                                        w.write('%s\t' % row)
                                else:
                                        w.write('%s\n' % row)        
                                        
        matching_hom_file = file_name[0:-4] + '_hom-matching.vcf'

        print 'Writing matching homozygous file...'
        
        with open(matching_hom_file, 'w') as w:
                for line in headers:
                        w.write(line + '\n')
                
                for line in matching_hom:
                        for index, row in enumerate(line):
                                if index != len(line) - 1:
                                        w.write('%s\t' % row)
                                else:
                                        w.write('%s\n' % row)
                                        

def main():

        file_names = glob.glob('*_merged.vcf')

        family_list = set([i.split('_')[0] for i in file_names])
        
        master_file = []
        
        #loops through on a family level, concatenates the information from each individual into the count file
        for family_id in family_list:
                family_files = glob.glob("%s*_merged.vcf" % family_id)

                family_number = family_files[0].split('_')[0]

                for individual in family_files:
                        data, headers = read_data(individual)
                        new_file, monoallelic_variants, other_mismatch, matching_het, matching_hom = parse_combinations(data, individual)
                        write_vcf_files(monoallelic_variants, other_mismatch, matching_het, matching_hom, individual, headers)
                        master_file.append(new_file)
        
                write_to_file(master_file, family_number)
                master_file = []
                                
if __name__ == "__main__":
        main()
