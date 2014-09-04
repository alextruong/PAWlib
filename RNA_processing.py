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

        print 'Data read into buffer'

        return headers, data


def filter_qual(data, input_qual, input_het, input_hom):
        """reads in gt_pl_gq, and qual_score, returns rows with qual_score >= 20 and gq >= 20 (het), gq >= 40 (hom) or values defined by user"""

        gq_score = [i[-1].split(':')[-1] for i in data]

        genotype = [i[-1].split(':')[0] for i in data]

        qual_score = [i[5] for i in data]


        qual_gq_filtered_rows = []

        for i in xrange(len(gq_score)):
                if genotype[i] in ['0/1', '1/0'] and float(qual_score[i]) >= input_qual:
                        if int(gq_score[i]) >= input_het:
                                qual_gq_filtered_rows.append(data[i])
                elif genotype[i] == '1/1' and float(qual_score[i]) >= input_qual:
                        if int(gq_score[i]) >= input_hom:
                                qual_gq_filtered_rows.append(data[i])

        print 'Filtering files by Quality Score and GQ'

        return qual_gq_filtered_rows


def filter_snp_indel(qual_gq_filtered_rows):
        """filter rows into snps, and indels"""

        length = []
        snps = []
        indels = []

        for i in xrange(len(qual_gq_filtered_rows)):
                temp = qual_gq_filtered_rows[i][4].split(',')

                for j in xrange(len(temp)):
                        length.append(len(temp[j]))

                if all(i == 1 for i in length) and len(qual_gq_filtered_rows[i][3]) == 1:
                        snps.append(qual_gq_filtered_rows[i])

                else:
                        indels.append(qual_gq_filtered_rows[i])

                length = []
                temp = []

        print 'SNPs and Multiallelic SNPs stored in buffer'

        return snps, indels


def write_processed_variants(file_name, snps, indels, headers, current_path, key_data):
        """defines naming schemes for snps/indels, then writes to file with original vcf format"""

        for i in key_data:
                if i[2] == file_name:
                        snp_name = str(i[0]) + '_' + str(i[1]) + '-' + i[-1] + '_RNA' + '_snp.vcf'
                        indel_name = str(i[0]) + '_' + str(i[1]) + '-' + i[-1] + '_RNA' + '_indel.vcf'

        with open(snp_name, 'w') as w:
                for line in headers:
                        w.write(line + '\n')

                for row in snps:
                        for index, column in enumerate(row):
                                if index != len(row) - 1:
                                        w.write('%s\t' % column)
                                else:
                                        w.write('%s\n' % column)

        with open(indel_name, 'w') as x:
                for line in headers:
                        x.write(line + '\n')
                                                                                                                                                                                             52,1-8        41%
                for row in indels:
                        for index, column in enumerate(row):
                                if index != len(row) - 1:
                                        x.write('%s\t' % column)
                                else:
                                        x.write('%s\n' % column)

        print 'Wrote SNP and Indel Files'


def main():

        if len(sys.argv) != 2:
                print 'Please call script: {0} <keyfile>'.format(sys.argv[0])
                sys.exit(1)

        first_input = str(raw_input("""Do you want default filtering settings? ('y' or 'n') """))
        toggle = first_input.lower()

        if toggle.isalpha():
                if toggle == 'y':
                        input_qual = 20
                        input_het = 20
                        input_hom = 40

                elif toggle == 'n':
                        try:
                                input_qual = int(raw_input("Quality score? "))
                                input_het = int(raw_input("Heterozygous GQ? "))
                                input_hom = int(raw_input("Homozygous GQ? "))
                                if input_qual < 0 or input_het < 0 or input_hom < 0:
                                        print 'Please enter positive integers'
                                        sys.exit(1)

                        except ValueError:
                                print 'Please enter only integers'
                                sys.exit(1)
                else:
                        print "Please enter only y or n"
                        sys.exit(1)
        else:
                print "Please enter only y or n"
                sys.exit(1)

        file_names = glob.glob('*.vcf')
        current_path = os.getcwd()


                                                                                                                                                                                             154,0-1       81%
        key_file = sys.argv[1]
        key_headers, key_data = read_data(key_file)

        family_ids = []
        for i in key_data:
                family_ids.append(i[0])
        family_key = set(family_ids)

        RNA_files = [i[2] for i in key_data]

        for index, file_name in enumerate(file_names):                  #parse in series, and write out final
                if file_name in RNA_files:
                        headers, data = read_data(file_name)
                        qual_gq_filtered_rows = filter_qual(data, input_qual, input_het, input_hom)
                        snps, indels = filter_snp_indel(qual_gq_filtered_rows)
                        write_processed_variants(file_name, snps, indels, headers, current_path, key_data)
                        print file_name, str(round((100*(index + 1)) / len(file_names), 3)) + '% complete\n'
                else:
                        continue

if __name__ == "__main__":
        main()
                                                                                                                                                                                             176,1-8       Bot
