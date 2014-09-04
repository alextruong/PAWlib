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

        matching_het_count = 0
        matching_hom_count = 0
        rna_hom_wes_het_count = 0
        rna_het_wes_hom_count = 0
        
        name = file_name
        file_name = []
        master_file = []
        
        monoallelic_variants = []
        other_mismatch = []
        matching_het = []
        matching_hom = []

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
                
                else:
                        print 'There are indels present in data set'
                        sys.exit(1)
        
        file_name.append(name)
        file_name.append(matching_het_count)
        file_name.append(matching_hom_count)
        file_name.append(rna_het_wes_hom_count)
        file_name.append(rna_hom_wes_het_count)
        print 'Counts written'
        
        
        return file_name, monoallelic_variants, other_mismatch, matching_het, matching_hom
        
def write_to_file(master_file, filename):
        
        file_name = filename.split('_')[0]
        print 'Writing File...'
        
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


def write_monoallelic_file(monoallelic_variants, file_name, headers):
        
        
        new_file_name = file_name[0:-4] + '_RNA-hom_WES-het.vcf'
        
        print 'Writing monoallelic file...'
        
        with open(new_file_name, 'w') as w:
                for line in headers:
                        w.write(line + '\n')
                
                for line in monoallelic_variants:
                        for row in line:
                                w.write('%s\t' % row)
                        w.write('\n')
                        
                        
                        
def write_other_mismatch(other_mismatch, file_name, headers):

        
        new_file_name_other = file_name[0:-4] + '_RNA-het_WES-hom.vcf'
        
        print 'Writing other mismatch file...'
        
        with open(new_file_name_other, 'w') as w:
                for line in headers:
                        w.write(line + '\n')
                
                for line in other_mismatch:
                        for row in line:
                                w.write('%s\t' % row)
                        w.write('\n')

def write_matching_het(matching_het, file_name, headers):

        
        new_file_name = file_name[0:-4] + '_het-matching.vcf'
        
        print 'Writing matching heterozygous file...'
        
        with open(new_file_name, 'w') as w:
                for line in headers:
                        w.write(line + '\n')
                
                for line in matching_het:
                        for row in line:
                                w.write('%s\t' % row)
                        w.write('\n')

def write_matching_hom(matching_hom, file_name, headers):

        
        new_file_name = file_name[0:-4] + '_hom-matching.vcf'

        print 'Writing matching homozygous file...'
        
        with open(new_file_name, 'w') as w:
                for line in headers:
                        w.write(line + '\n')
                
                for line in matching_hom:
                        for row in line:
                                w.write('%s\t' % row)
                        w.write('\n')
        
def main():

        file_names = glob.glob('*_merged.vcf')


        family_list = set([i.split('_')[0] for i in file_names])
        
        master_file = []
        
        for family_id in family_list:
                family_files = glob.glob("%s*_merged.vcf" % family_id)

                for individual in family_files:
                        data, headers = read_data(individual)
                        new_file, monoallelic_variants, other_mismatch, matching_het, matching_hom = parse_combinations(data, individual)
                        write_monoallelic_file(monoallelic_variants, individual, headers)
                        write_other_mismatch(other_mismatch, individual, headers)
                        write_matching_het(matching_het, individual, headers)
                        write_matching_hom(matching_hom, individual, headers)


                        master_file.append(new_file)
        
                        write_to_file(master_file)
                        master_file = []
                                
if __name__ == "__main__":
        main()
