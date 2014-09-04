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


def compare(ref_data, query_data, file_name):

        query_out_data = []
        ref_out_data = []
        multiallelic_locations = []

        lquery_data = query_data
        lref_data = ref_data

        ref_alt_key = [i[4] for i in lref_data]
        query_alt_key = [i[4] for i in lquery_data]

        ref_data_key = [tuple((i[0], i[1], i[3])) for i in lref_data]
        query_data_key = [tuple((i[0], i[1], i[3])) for i in lquery_data]

        common = set(ref_data_key).intersection( set(query_data_key) )

        for i in common:
                ref_index = ref_data_key.index(i)
                query_index = query_data_key.index(i)

                if query_alt_key[query_index] == ref_alt_key[ref_index]:
                        query_out_data.append(query_data[query_index])
                        ref_out_data.append(ref_data[ref_index])

                else:
                        query_temp = query_alt_key[query_index].split(',')
                        ref_temp = ref_alt_key[ref_index].split(',')

                        if any(j in query_temp for j in ref_temp):
                                query_out_data.append(query_data[query_index])
                                ref_out_data.append(ref_data[ref_index])
                                multiallelic_locations.append(query_data_key[query_index])
        
        return query_out_data, ref_out_data, multiallelic_locations


def write_to_file(ref_file, query_file, ref_headers, query_headers, query_out_data, ref_out_data, multiallelic_locations):

        new_ref_file = ref_file[0:-4] + '_shared.vcf'
        new_query_file = query_file[0:-4] + '_shared.vcf'

        with open(new_ref_file, 'w') as w:
                for line in ref_headers:
                        w.write(line + '\n')

                for row in ref_out_data:
                        for index, column in enumerate(row):
                                if index != len(row) - 1:
                                        w.write('%s\t' % column)
                                else:
                                        w.write('%s\n' % column)

        with open(new_query_file, 'w') as w:
                for line in query_headers:
                        w.write(line + '\n')

                for row in query_out_data:
                        for index, column in enumerate(row):
                                if index != len(row) - 1:
                                        w.write('%s\t' % column)
                                else:
                                        w.write('%s\n' % column)


        hold = new_ref_file.split('_')[0:2]
        wobblefile = hold[0] + '_' + hold[1]

        with open('%s_incomplete_match_coordinates.txt' % wobblefile, 'w') as w:
                for line in multiallelic_locations:
                        w.write(str(line) + '\n')


def findbcf():

        print 'Locating bcftools installation...'

        find_bcftools_path = commands.getstatusoutput('find ~/ -name bcftools -executable -type f -print 2>/dev/null')

        bcfdir = find_bcftools_path[1]

        print 'Done!'

        return bcfdir


def bcfmerge(rna, wes, bcfdir):
        """bgzip and tabix output files from gsearch, and use bcftools to merge matching files"""

        bgzipcommand = "bgzip %s; bgzip %s" % (rna, wes)

        rnatemp = rna[0:-4] + '_sorted.vcf'
        westemp = wes[0:-4] + '_sorted.vcf'

        rnatabix = rna + '.gz'
        westabix = wes + '.gz'
        
        bgzipcommand = "bgzip %s; bgzip %s" % (rnatemp, westemp)

        tabixcommand = "tabix -p vcf %s; tabix -p vcf %s" % (rnatabix, westabix)

        sortcommand = 'vcf-sort %s > %s; vcf-sort %s > %s' % (rna, rnatemp, wes, westemp)
        
        os.system(sortcommand)

        #print 'bgzipping gsearch output files...'
        os.system(bgzipcommand)
        #print 'Done!'
        #print '-----------------------------------------------'
        #print 'Indexing with tabix...'
        os.system(tabixcommand)
        #print 'Done!'
        #print '-----------------------------------------------'

        mergename = rna.split('_')[0] + '_' + rna.split('_')[1] + '_' + rna.split('_')[2] + '_WES_' + rna.split('_')[3] + '_' + rna.split('_')[4][0:-4]+ '_sorted_merged.vcf'

        mergecommand = '%s merge -O v %s %s > %s' % (bcfdir, rnatabix, westabix, mergename)

        #print 'Merging gsearch output files...'
        os.system(mergecommand)
        #print 'Done!'

        #print '-----------------------------------------------'
        
        return

def main():

        RNA_files = glob.glob('*_RNA_snp.vcf')
        WES_files = glob.glob('*_WES_snp.vcf')

        if len(RNA_files) != len(WES_files):
                print 'Missing files'
                sys.exit(1)
        counter = 0
        print 'Processing file comparisons...'
        sys.stdout.write('Progress: 0% \r')
        sys.stdout.flush()
        for rna_name in RNA_files:
                for wes_name in WES_files:
                        if rna_name.split('-')[0] == wes_name.split('-')[0]:
                                ref_file = rna_name
                                query_file = wes_name

                                ref_headers, ref_data = read_data(ref_file)
                                query_headers, query_data = read_data(query_file)

                                query_out_data, ref_out_data, multiallelic_locations = compare(ref_data, query_data, ref_file)
                                write_to_file(ref_file, query_file, ref_headers, query_headers, query_out_data, ref_out_data, multiallelic_locations)

                counter += 1
                progress = counter / len(RNA_files) * 100

                if progress != 100:
                        sys.stdout.write('Progress: %d%%   \r' % progress)
                        sys.stdout.flush()
                else:
                        sys.stdout.write('Progress: 100%\n')
                        print 'Done!'

        print '----------------------------------------------------'

        bcfdir = findbcf()

        rnapool = glob.glob('*_RNA_snp_shared.vcf')
        wespool = glob.glob('*_WES_snp_shared.vcf')

        counter = 0
        print 'Merging matching files...'
        sys.stdout.write('Progress: 0% \r')
        sys.stdout.flush()
        for rna in rnapool:
                for wes in wespool:
                        if rna.split('-')[0] == wes.split('-')[0]:
                                rna_match = rna
                                wes_match = wes

                                bcfmerge(rna, wes, bcfdir)

                counter += 1
                progress = counter / len(rnapool) * 100

                if progress != 100:
                        sys.stdout.write('Progress: %d%%   \r' % progress)
                        sys.stdout.flush()
                else:
                        sys.stdout.write('Progress: 100%\n')
                        print 'Done!'

        #outdir = makefolder()

        #print 'Organizing files...'

        #os.system('cd ' + outdir + '; mkdir matching; mv output1.vcf matching; mv output2.vcf matching; mkdir nonmatching; mv unmatched.output1.vcf nonmatching; mv unmatched.output2.vcf nonmatching; mkdir merged')

        print 'Complete!'



if __name__ == "__main__":
        main()



