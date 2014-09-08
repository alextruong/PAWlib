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


def monoallelic_comparison(monoallelic_query, full_merge_reference):
        """same method as gsearchcompare, it compares tuples of positional information to match"""
        
        lmonoallelic_query = list(monoallelic_query)
        lfull_merge_reference = list(full_merge_reference)

        monoallelic_alt_key = [i[4] for i in lmonoallelic_query]
        full_merge_alt_key = [i[4] for i in lfull_merge_reference]

        monoallelic_variant_key = [tuple((i[0], i[1], i[3])) for i in lmonoallelic_query]
        full_merge_variant_key = [tuple((i[0], i[1], i[3])) for i in lfull_merge_reference]

        common = set(full_merge_variant_key).intersection( set(monoallelic_variant_key) )

        unique = set(monoallelic_variant_key).difference( set(full_merge_variant_key) )
        

        
        query_common_monoallelic = []
        ref_common_monoallelic = []
        query_unique = []

        for i in common:

                query_index = monoallelic_variant_key.index(i)
                reference_index = full_merge_variant_key.index(i)
                
                if monoallelic_alt_key[query_index] == full_merge_alt_key[reference_index]:
                        query_common_monoallelic.append(lmonoallelic_query[query_index])
                        ref_common_monoallelic.append(lfull_merge_reference[reference_index])

                else:
                        query_temp = monoallelic_alt_key[query_index].split(',')
                        ref_temp = full_merge_alt_key[reference_index].split(',')
                        
                        if any(j in query_temp for j in ref_temp):
                                query_common_monoallelic.append(lmonoallelic_query[query_index])
                                ref_common_monoallelic.append(lfull_merge_reference[reference_index])

        for i in unique:
                query_index = monoallelic_variant_key.index(i)
                query_unique.append(lmonoallelic_query[query_index])
                        

        return query_common_monoallelic, ref_common_monoallelic, query_unique

def write_monoallelic_to_file(query_common_headers, query_common_monoallelic, ref_common_monoallelic, ref_common_headers, query_unique, other_child, filename):
        """standard writing to file"""
        name_split = filename[0:-4]
        query_common_full_name = name_split + '_monoallelic_common.vcf'
        ref_common_full_name = name_split + '_monoallelic_common_in_%s.vcf' % other_child
        query_unique_full_name = name_split + '_monoallelic_unique.vcf'

        with open(query_common_full_name, 'w') as w:
                for line in query_common_headers:
                        w.write(line + '\n')

                for row in query_common_monoallelic:
                        for index, column in enumerate(row):
                                if index != len(row) - 1:
                                        w.write('%s\t' % column)
                                else:
                                        w.write('%s\n' % column)


        with open(ref_common_full_name, 'w') as w:
                for line in ref_common_headers:
                        w.write(line + '\n')

                for row in ref_common_monoallelic:
                        for index, column in enumerate(row):
                                if index != len(row) - 1:
                                        w.write('%s\t' % column)
                                else:
                                        w.write('%s\n' % column)

        with open(query_unique_full_name, 'w') as w:
                for line in query_common_headers:
                        w.write(line + '\n')

                for row in query_unique:
                        for index, column in enumerate(row):
                                if index != len(row) - 1:
                                        w.write('%s\t' % column)
                                else:
                                        w.write('%s\n' % column)

def make_gene_list(annotated_query_file, monoallelic_file):
        """makes a gene list based on VEP annotation"""

        query_genes = set([str(i[-1].split(';')[-1].split('=')[-1]) for i in annotated_query_file if i[-1].split(';')[-1].split('=')[0] == 'SYMBOL'])

        filename = monoallelic_file[0:-4] + '_genelist.txt'

        with open(filename, 'w') as w:
                for query_gene in query_genes:
                        w.write(query_gene + '\n')

        return query_genes


def gene_comparison(*args):
        """Function takes annotated file and/or query input gene list and returns dictionary with genes and chrom positions"""

        annotated_reference_merge = args[0]
        ref_gene_dict = {}      
        common_genes = {}
                
        if len(args) >= 1:
                for line in annotated_reference_merge:
                
                        if line[-1].split(';')[-1].split('=')[0] == 'SYMBOL':
                                key = line[-1].split(';')[-1].split('=')[-1]
                                value = line[1]
        
                                if key not in ref_gene_dict:
                                        ref_gene_dict[key] = []
                                        ref_gene_dict[key].append(value)
                
                                else:
                                        if value not in ref_gene_dict[key]:
                                                ref_gene_dict[key].append(value)
                                        else:
                                                continue

        if len(args) == 2:
                query_genes = args[1]
                for gene in query_genes:
                        if gene in ref_gene_dict:
                                common_genes[gene] = ref_gene_dict[gene]        

                        else:
                                continue

        return ref_gene_dict, common_genes

def write_gene_comparison_dictionary(gene_dictionary, filename):
        """standard write out"""
        with open(filename[0:-4] + '_monoallelic_gene_dictionary.txt', 'w') as w:
                for key in gene_dictionary:
                        positions = ','.join(gene_dictionary[key])
                        w.write(str(key) + '-' + positions + '\n')

def write_gene_comparison_dictionary_common(common_genes, filename):
        """standard write out"""
        with open(filename[0:-4] + '_monoallelic_common_gene_dictionary.txt', 'w') as w:
                for key in common_genes:
                        positions = ','.join(common_genes[key])
                        w.write(str(key) + '-' + positions + '\n')

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

def main():
        
        #script loops through on a family level
        file_names = glob.glob('*_merged.vcf')
        family_list = set([i.split('_')[0] for i in file_names])
        
        master_file = []

        family_files_monoallelic = glob.glob("*_RNA-hom_WES-het.vcf")
        family_files_monoallelic_annotated = glob.glob("*_RNA-hom_WES-het_annotated.vcf")
        
        for family_id in family_list:
                
                family_files = [i for i in family_files_monoallelic if i.split('_')[0] == family_id]
                
                for i in family_files:
                        if i.split('-')[0][-2:] == 's1':
                                sibling_annotated = i
                        elif i.split('-')[0][-2:] == 'p1':
                                proband_annotated = i
                        else:
                                continue


                for monoallelic_file in family_files:           #sibling-proband monoallelic comparison

                        individual = monoallelic_file.split('-')[0][-2:]
                        if individual  == 'p1':

                                full_reference_file = sibling_annotated
                                proband_common_headers, monoallelic_data = read_data(monoallelic_file)
                                ref_headers, full_reference_file_data = read_data(full_reference_file)

                                other_child = 's1'
                                proband_common_monoallelic, sibling_common_monoallelic, proband_unique = monoallelic_comparison(monoallelic_data, full_reference_file_data)
                                write_monoallelic_to_file(proband_common_headers, proband_common_monoallelic, sibling_common_monoallelic, ref_headers, proband_unique, other_child, monoallelic_file)
        
                        elif individual == 's1':
                                
                                full_reference_file = proband_annotated
                                sibling_common_headers, monoallelic_data = read_data(monoallelic_file)
                                ref_headers, full_reference_file_data = read_data(full_reference_file)
                                
                                other_child = 'p1'
                                sibling_common_monoallelic, proband_common_monoallelic, sibling_unique = monoallelic_comparison(monoallelic_data, full_reference_file_data)
                                write_monoallelic_to_file(sibling_common_headers, sibling_common_monoallelic, proband_common_monoallelic, ref_headers, sibling_unique, other_child, monoallelic_file)
                        else:
                                continue
                
                family_files_annotated = [i for i in family_files_monoallelic_annotated if i.split('_')[0] == family_id]

                for i in family_files_annotated:
                        if i.split('-')[0][-2:] == 's1':
                                sibling_annotated = i
                        elif i.split('-')[0][-2:] == 'p1':
                                proband_annotated = i
                        else:
                                continue

                for monoallelic_file in family_files_monoallelic_annotated:             #write gene lists from monoallelic and gene comparison
                
                       
                        individual = monoallelic_file.split('-')[0][-2:]
                        
                        #individual = monoallelic_file.split('-')[0][-2:], gene dictionary called with two arguments (query genes and gene list)
                        if individual == 'p1':
                                headers, monoallelic_data = read_data(monoallelic_file) 
                                query_genes = make_gene_list(monoallelic_data, monoallelic_file)
                                annotated_reference_merge = sibling_annotated
                                headers, sibling_annotated_data = read_data(annotated_reference_merge)
                                ref_gene_dict, common_genes = gene_comparison(sibling_annotated_data, query_genes)
                                write_gene_comparison_dictionary(ref_gene_dict, monoallelic_file)
                                write_gene_comparison_dictionary_common(common_genes, monoallelic_file)

                        elif individual == 's1':
                                headers, monoallelic_data = read_data(monoallelic_file) 
                                query_genes = make_gene_list(monoallelic_data, monoallelic_file)
                                annotated_reference_merge = proband_annotated
                                headers, proband_annotated_data = read_data(annotated_reference_merge)
                                ref_gene_dict, common_genes = gene_comparison(proband_annotated_data, query_genes)
                                write_gene_comparison_dictionary(ref_gene_dict, monoallelic_file)
                                write_gene_comparison_dictionary_common(common_genes, monoallelic_file)
                        #gene_dictionary called with one argument (gene list)
                        else:
                                headers, monoallelic_data = read_data(monoallelic_file)
                                ref_gene_dict, common_genes = gene_comparison(monoallelic_data)
                                write_gene_comparison_dictionary(ref_gene_dict, monoallelic_file)
                
                mother_WES = "%s_mo-f_WES_snp.vcf" % family_id
                father_WES = "%s_fa-m_WES_snp.vcf" % family_id
                headers, mother_WES_data = read_data(mother_WES)
                headers, father_WES_data = read_data(father_WES)
                for monoallelic_file in family_files:           #traces lineage from spawn to parent WES
                        
                        individual = monoallelic_file.split('-')[0][-2:]

                        for i in family_files:
                                if i.split('-')[0][-2:] == 's1':
                                        unique_monoallelic_sibling = i
                                elif i.split('-')[0][-2:] == 'p1':
                                        unique_monoallelic_proband = i
                                
                                else:
                                        continue
                        
                        if individual == 'p1':
                                headers, unique_monoallelic_proband_data = read_data(unique_monoallelic_proband)
                                mother_zygosity = trace_lineage(unique_monoallelic_proband_data, mother_WES_data)
                                father_zygosity = trace_lineage(unique_monoallelic_proband_data, father_WES_data)
                                write_lineage_files(individual, unique_monoallelic_proband_data, mother_zygosity, father_zygosity, monoallelic_file)
                        elif individual == 's1':
                                headers, unique_monoallelic_sibling_data = read_data(unique_monoallelic_sibling)
                                mother_zygosity = trace_lineage(unique_monoallelic_sibling_data, mother_WES_data)
                                father_zygosity = trace_lineage(unique_monoallelic_sibling_data, father_WES_data)
                                write_lineage_files(individual, unique_monoallelic_sibling_data, mother_zygosity, father_zygosity, monoallelic_file)
                        else:
                                continue
        
        monoallelic_vcfs = glob.glob('*monoallelic*.vcf')
        
        sort = 'vcf-sort %s > %s'
        
        for vcf in monoallelic_vcfs:
                new_name = str(vcf[0:-4]) + '_sorted.vcf' 
                os.system(sort % (vcf, new_name))
        
        updated_monoallelic_vcfs = glob.glob('*monoallelic*.vcf')
        for vcf in monoallelic_vcfs:
                if not vcf.split('_')[-1][0:-4] == 'sorted':
                        os.system('rm %s' % vcf)
                else:
                        continue
if __name__ == "__main__":
        main()
