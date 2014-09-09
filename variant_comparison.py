#!/usr/bin/env python


from __future__ import division
import os
import sys
import commands
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


        return headers, data


def create_gene_list(annotated_query_file_data):
        """makes a gene list based on VEP annotation"""

        monoallelic_genes = set([str(i[-1].split(';')[-1].split('=')[-1]) for i in annotated_query_file_data if i[-1].split(';')[-1].split('=')[0] == 'SYMBOL'])

        return monoallelic_genes


def create_gene_variant_dictionary(annotated_file, filename):
        """Function takes annotated file and returns dictionary with genes and chrom positions"""

        gene_var_dict = {}      
                
        for line in annotated_file:
                
                if line[-1].split(';')[-1].split('=')[0] == 'SYMBOL':
                        key = line[-1].split(';')[-1].split('=')[-1]
                        value = line[1]
        
                        if key not in gene_var_dict:
                                gene_var_dict[key] = []
                                gene_var_dict[key].append(value)
                
                        else:
                                if value not in gene_var_dict[key]:
                                        gene_var_dict[key].append(value)
                                else:
                                        continue

        return gene_var_dict

def unique_and_common_monoallelic_genes(proband_monoallelic_genes, proband_variant_dict, proband_file_name, sibling_monoallelic_genes, sibling_variant_dict, sibling_file_name):

        proband_unique_monoallelic_genes = [unique_gene for unique_gene in proband_monoallelic_genes if unique_gene not in sibling_monoallelic_genes]
        sibling_unique_monoallelic_genes = [unique_gene for unique_gene in sibling_monoallelic_genes if unique_gene not in proband_monoallelic_genes]

        #proband_unique_dict = {key:value for key, value in proband_variant_dict.iteritems() if key in proband_unique_monoallelic_genes}
        #sibling_unique_dict = {key:value for key, value in sibling_variant_dict.iteritems() if key in sibling_unique_monoallelic_genes}
        proband_unique_dict = dict((k, v) for (k, v) in proband_variant_dict.iteritems() if k in proband_unique_monoallelic_genes)
        sibling_unique_dict = dict((k, v) for (k, v) in sibling_variant_dict.iteritems() if k in sibling_unique_monoallelic_genes)
        
        common_genes = set(proband_monoallelic_genes).intersection(sibling_monoallelic_genes)

        proband_unique_name = proband_file_name[0:-4] + '_unique_monoallelic_genelist_variants.txt'
        sibling_unique_name = sibling_file_name[0:-4] + '_unique_monoallelic_genelist_variants.txt'
        common_file_name = proband_file_name.split('_')[0] + '_common_p1_s1_monoallelic_genelist_variants.txt'

        with open(proband_unique_name, 'w') as w:
                w.write('#gene' + '\t' + "variants\n")

                for unique_gene in proband_unique_monoallelic_genes:
                        variants = ','.join(proband_variant_dict[unique_gene])
                        w.write(unique_gene + '\t' + variants + '\n')

        with open(sibling_unique_name, 'w') as w:
                w.write('#gene' + '\t' + "variants\n")

                for unique_gene in sibling_unique_monoallelic_genes:
                        variants = ','.join(sibling_variant_dict[unique_gene])
                        w.write(unique_gene + '\t' + variants + '\n')
        

        with open(common_file_name, 'w') as w:
                headers = ['#gene', 'p1 variants', 's1 variants']
                w.writelines('\t'.join(headers))
                w.write('\n')
                for common_gene in common_genes:
                        proband_positions = ','.join(proband_variant_dict[common_gene])
                        sibling_positions = ','.join(sibling_variant_dict[common_gene])

                        w.write(str(common_gene) + '\t' + proband_positions + '\t' + sibling_positions + '\n')

        return proband_unique_monoallelic_genes, proband_unique_dict, sibling_unique_monoallelic_genes, sibling_unique_dict

def make_unique_monoallelic_vcf(monoallelic_file_data, monoallelic_file_headers, positions, monoallelic_file_name):
        lreference_monoallelic_data = list(monoallelic_file_data)
        reference_key = [tuple((i[0], i[1])) for i in lreference_monoallelic_data]

        positions_key = [tuple((position.split(':')[0], position.split(':')[1])) for position in positions]
        
        unique_monoallelic_rows = []
        for position in positions_key:
                if position in reference_key:
                        reference_index = reference_key.index(position)
                        unique_monoallelic_rows.append(lreference_monoallelic_data[reference_index])
                else:
                        continue

        new_vcf = monoallelic_file_name[0:-4] + '_unique.vcf'
        with open(new_vcf, 'w') as w:
                for line in monoallelic_file_headers:
                        w.write(line + '\n')

               
                for index, unique_row in enumerate(unique_monoallelic_rows):
                        if index != len(unique_monoallelic_rows) - 1:
                                for column_index, column in enumerate(unique_row):
                                        if column_index != len(unique_row) - 1:
                                                w.write(column + '\t')
                                        else:
                                                w.write(column + '\n')
                        else:
                                for column_index, column in enumerate(unique_row):
                                        if column_index != len(unique_row) - 1:
                                                w.write(column + '\t')
                                        else:
                                                w.write(column)   
                         

def trace_lineage(unique_variants_dict, reference_WES):
        """same function as gsearch compare merge"""

        positions = [single_key for key in unique_variants_dict for single_key in unique_variants_dict[key]]
        query_key = [tuple((position.split(':')[0], position.split(':')[1])) for position in positions]

        lreference_WES = list(reference_WES)
        WES_variants_key = [tuple((i[0], i[1])) for i in lreference_WES]

        zygosity = []
        for position in query_key:
                if position in WES_variants_key:
                        reference_index = WES_variants_key.index(position)
                        zygosity.append(lreference_WES[reference_index][-1].split(':')[0])

                else:
                        zygosity.append('N/A')                                  

        return zygosity, positions

def write_lineage_files(unique_dict, individual, positions, mother_zygosity, father_zygosity, filename):
        """writes out file with variant, and zygosity in mother and father"""
        header_lines = [('#' + individual) + ' genes', 'variant', 'mother zygosity', 'father zygosity']

        with open(filename[0:-4] + '_unique_lineage_parent_WES.txt', 'w') as w:
                for header in header_lines:
                        w.write(header + '\t')
                w.write('\n')

                for index, zygosity in enumerate(mother_zygosity):
                        for key in unique_dict:
                                for values in unique_dict[key]:
                                        if positions[index] in values:
                                                gene = key
                        w.write(gene + '\t' + positions[index] + '\t')
                        w.write(mother_zygosity[index] + '\t')
                        w.write(father_zygosity[index] + '\n')

def write_gene_dict(gene_var_dict, filename):
        new_gene_var_file = filename[0:-4] + '_full_monoallelic_genelist_variants.txt'
        
        with open(new_gene_var_file, 'w') as w:
                w.write('#gene\t' + 'variants\n')
                
                for index, key in enumerate(gene_var_dict):
                        positions = ','.join(gene_var_dict[key])
                        if index != len(gene_var_dict) - 1:
                                w.write(key + '\t' + positions + '\n')
                        else:
                                w.write(key + '\t' + positions)

def main():
        
        #script loops through on a family level
        file_names = glob.glob('*_merged.vcf')
        family_list = set([i.split('_')[0] for i in file_names])
        
        master_file = []

        all_monoallelic_files = glob.glob("*_RNA-hom_WES-het.vcf")
        all_monoallelic_annotated_files = glob.glob("*_RNA-hom_WES-het_annotated.vcf")
        
        for family_id in family_list:
                
                family_specific_monoallelic_annotated_files = [i for i in all_monoallelic_annotated_files if i.split('_')[0] == family_id]

                for monoallelic_annotated_file in family_specific_monoallelic_annotated_files:
                        
                        individual = monoallelic_annotated_file.split('-')[0][-2:]

                        if individual == 'p1':
                                headers, data = read_data(monoallelic_annotated_file)
                                proband_file_name = monoallelic_annotated_file
                                proband_monoallelic_genes = create_gene_list(data)
                                proband_monoallelic_genes_variants = create_gene_variant_dictionary(data, proband_file_name)


                        elif individual == 's1':
                                headers, data = read_data(monoallelic_annotated_file)
                                sibling_file_name = monoallelic_annotated_file
                                sibling_monoallelic_genes = create_gene_list(data)
                                sibling_monoallelic_genes_variants = create_gene_variant_dictionary(data, sibling_file_name)

                        else:
                                headers, data = read_data(monoallelic_annotated_file)
                                gene_var_dict = create_gene_variant_dictionary(data, monoallelic_annotated_file)
                                write_gene_dict(gene_var_dict, monoallelic_annotated_file)
                                
                proband_unique_monoallelic_genes, proband_unique_dict, sibling_unique_monoallelic_genes, sibling_unique_dict = unique_and_common_monoallelic_genes(proband_monoallelic_genes, proband_monoallelic_genes_variants, proband_file_name, sibling_monoallelic_genes, sibling_monoallelic_genes_variants, sibling_file_name)


                family_specific_monoallelic_files = [i for i in all_monoallelic_files if i.split('_')[0] == family_id]
                mother_WES = "%s_mo-f_WES_snp.vcf" % family_id
                father_WES = "%s_fa-m_WES_snp.vcf" % family_id

                headers, mother_WES_data = read_data(mother_WES)
                headers, father_WES_data = read_data(father_WES)

                for monoallelic_file in family_specific_monoallelic_files:           #traces lineage from spawn to parent WES
                        
                        individual = monoallelic_file.split('-')[0][-2:]
                        
                        if individual == 'p1':
                                monoallelic_file_headers, monoallelic_file_data = read_data(monoallelic_file)
                                mother_zygosity, positions = trace_lineage(proband_unique_dict, mother_WES_data)
                                father_zygosity, positions = trace_lineage(proband_unique_dict, father_WES_data)
                                write_lineage_files(proband_unique_dict, individual, positions, mother_zygosity, father_zygosity, monoallelic_file)
                                make_unique_monoallelic_vcf(monoallelic_file_data, monoallelic_file_headers, positions, monoallelic_file)

                        elif individual == 's1':
                                monoallelic_file_headers, monoallelic_file_data = read_data(monoallelic_file)
                                mother_zygosity, positions = trace_lineage(sibling_unique_dict, mother_WES_data)
                                father_zygosity, positions = trace_lineage(sibling_unique_dict, father_WES_data)
                                write_lineage_files(sibling_unique_dict, individual, positions, mother_zygosity, father_zygosity, monoallelic_file)
                                make_unique_monoallelic_vcf(monoallelic_file_data, monoallelic_file_headers, positions, monoallelic_file)
                       
                        else:
                                continue


                unique_vcfs = glob.glob('*_unique.vcf')

                sort = 'vcf-sort %s > %s'
        
                for vcf in unique_vcfs:
                        new_name = str(vcf[0:-4]) + '_sorted.vcf' 
                        os.system(sort % (vcf, new_name))
                
                for vcf in unique_vcfs:
                        os.system('rm %s' % vcf)


if __name__ == "__main__":
        main()
