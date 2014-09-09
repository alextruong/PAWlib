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

        with open(filename[0:-4] + '_monoallelic_genelist_variants.txt', 'w') as w:
                for key in gene_var_dict:
                        positions = ','.join(gene_var_dict[key])
                        w.write(str(key) + '\t' + positions + '\n')

        return gene_var_dict



def unique_and_common_monoallelic_genes(proband_monoallelic_genes, proband_variant_dict, proband_file_name, sibling_monoallelic_genes, sibling_variant_dict, sibling_file_name):

        proband_unique_monoallelic_genes = [unique_gene for unique_gene in proband_monoallelic_genes if unique_gene not in sibling_monoallelic_genes]
        sibling_unique_monoallelic_genes = [unique_gene for unique_gene in sibling_monoallelic_genes if unique_gene not in proband_monoallelic_genes]

        common_genes = set(proband_monoallelic_genes).intersection(sibling_monoallelic_genes)

        proband_unique_name = proband_file_name[0:-4] + '_unique_monoallelic_genelist.txt'
        sibling_unique_name = sibling_file_name[0:-4] + '_unique_monoallelic_genelist.txt'
        common_file_name = proband_file_name.split('_')[0] + '_common_proband_sibling_monoallelic_genelist.txt'

        with open(proband_unique_name, 'w') as w:

                w.writelines('\n'.join(proband_unique_monoallelic_genes))

        with open(sibling_unique_name, 'w') as w:
                
                w.writelines('\n'.join(sibling_unique_monoallelic_genes))

        with open(common_file_name, 'w') as w:

                w.writelines('\n'.join(common_genes))

        with open(common_variants_name, 'w') as w:
                headers = ['Gene', 'p1 variants', 's1 variants']
                w.writelines('\t'.join(headers))
                w.write('\n')
                for common_gene in common_genes:
                        proband_positions = ','.join(proband_variant_dict[common_gene])
                        sibling_positions = ','.join(sibling_variant_dict[common_gene])

                        w.write(str(common_gene) + '\t' + proband_positions + '\t' + sibling_positions + '\n'))



def trace_lineage(query_variants, reference_WES):
        """same function as gsearch compare merge"""

        lquery_variants = list(query_variants)
        lreference_WES = list(reference_WES)
        query_variants_alt = [i[4] for i in lquery_variants]
        WES_alt = [i[4] for i in lreference_WES]

        query_variants_key = [tuple((i[0], i[1], i[3])) for i in lquery_variants]
        WES_variants_key = [tuple((i[0], i[1], i[3])) for i in lreference_WES]

        zygosity = []
        for i in query_variants_key:
                if i in WES_variants_key:
                        query_index = query_variants_key.index(i)
                        reference_index = WES_variants_key.index(i)
                        
                        if query_variants_alt[query_index] == WES_alt[reference_index]:
                                zygosity.append(lreference_WES[reference_index][-1].split(':')[0])

                        elif query_variants_alt[query_index] != WES_alt[reference_index]:
                                query_temp = query_variants_alt[query_index].split(',')
                                ref_temp = WES_alt[reference_index].split(',')
                        
                                if any(j in query_temp for j in ref_temp):
                                        zygosity.append(lreference_WES[reference_index][-1].split(':')[0])
                else:
                        zygosity.append('N/A')                                  

        return zygosity

def write_lineage_files(individual, query_variants, mother_zygosity, father_zygosity, filename):
        """writes out file with variant, and zygosity in mother and father"""
        header_lines = [individual, 'Mother', 'Father']

        with open(filename[0:-4] + '_traced_WES.txt', 'w') as w:
                for header in header_lines:
                        w.write(header + '\t')
                w.write('\n')

                for index, zygosity in enumerate(mother_zygosity):
                        query_variant_info = query_variants[index][0] +  '_' + query_variants[index][1] + '_'  + query_variants[index][3]+ '_'  + query_variants[index][4]
                        w.write(query_variant_info + '\t')
                        w.write(mother_zygosity[index] + '\t')
                        w.write(father_zygosity[index] + '\n')

def write_gene_list(gene_file_name, gene_list):
        monoallelic_file = gene_file_name[0:-4] + '_full_monoallelic_genelist.txt'

        with open(monoallelic_file, 'w') as w:
                for monoallelic_gene in gene_list:
                        w.write(monoallelic_gene + '\n')

def main():
        
        #script loops through on a family level
        file_names = glob.glob('*_merged.vcf')
        family_list = set([i.split('_')[0] for i in file_names])
        
        master_file = []

        monoallelic_files = glob.glob("*_RNA-hom_WES-het.vcf")
        monoallelic_annotated_files = glob.glob("*_RNA-hom_WES-het_annotated.vcf")
        
        for family_id in family_list:
                
                family_specific_monoallelic_files = [i for i in monoallelic_files if i.split('_')[0] == family_id]
                mother_WES = "%s_mo-f_WES_snp.vcf" % family_id
                father_WES = "%s_fa-m_WES_snp.vcf" % family_id

                headers, mother_WES_data = read_data(mother_WES)
                headers, father_WES_data = read_data(father_WES)

                for monoallelic_file in family_specific_monoallelic_files:           #traces lineage from spawn to parent WES
                        
                        individual = monoallelic_file.split('-')[0][-2:]
                        
                        if individual == 'p1':
                                headers, monoallelic_proband_data = read_data(monoallelic_file)
                                mother_zygosity = trace_lineage(monoallelic_proband_data, mother_WES_data)
                                father_zygosity = trace_lineage(monoallelic_proband_data, father_WES_data)
                                write_lineage_files(individual, monoallelic_proband_data, mother_zygosity, father_zygosity, monoallelic_file)

                        elif individual == 's1':
                                headers, monoallelic_sibling_data = read_data(monoallelic_file)
                                mother_zygosity = trace_lineage(monoallelic_sibling_data, mother_WES_data)
                                father_zygosity = trace_lineage(monoallelic_sibling_data, father_WES_data)
                                write_lineage_files(individual, monoallelic_sibling_data, mother_zygosity, father_zygosity, monoallelic_file)
                       
                        else:
                                continue

                for monoallelic_annotated_file in monoallelic_annotated_files:
                        
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
                                monoallelic_genes = create_gene_list(data)
                                write_gene_list(monoallelic_annotated_file, monoallelic_genes)


                unique_and_common_monoallelic_genes(proband_monoallelic_genes, proband_monoallelic_genes_variants, proband_file_name, sibling_monoallelic_genes, sibling_monoallelic_genes_variants sibling_file_name)



        
        # monoallelic_vcfs = glob.glob('*monoallelic*.vcf')
        
        # sort = 'vcf-sort %s > %s'
        
        # for vcf in monoallelic_vcfs:
        #         new_name = str(vcf[0:-4]) + '_sorted.vcf' 
        #         os.system(sort % (vcf, new_name))
        
        # updated_monoallelic_vcfs = glob.glob('*monoallelic*.vcf')
        # for vcf in monoallelic_vcfs:
        #         if not vcf.split('_')[-1][0:-4] == 'sorted':
        #                 os.system('rm %s' % vcf)
        #         else:
        #                 continue
if __name__ == "__main__":
        main()
