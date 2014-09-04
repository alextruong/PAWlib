#!/usr/bin/env python


import os
import sys
import commands

def vepcaller(filepath, family_id):
        
        print 'Locating VEP installation...'

        veploc = commands.getstatusoutput('find ~/ -name variant_effect_predictor.pl -print 2>/dev/null')


        print 'Done!'

        print 'Running VEP...'
        noextension = file[0:-4]

        command = veploc + "/ perl variant_effect_predictor.pl --port 3337 --cache --refseq -i --symbol " + family_id + " -o " + filepath + "/" + noextension + "_annotated.vcf")

        print 'Done!'

        return


def main():
        filepath = str(raw_input("Please specify directory where file to annotate is located (full path if necessary): "))
        
        parent_files = glob.glob('*_WES_snp.vcf')
        monoallelic_files = glob.glob('*_RNA-hom_WES-het.vcf')
        full_merge = glob.glob('*_merged.vcf')

        
        family_list = []
        for family in parent_files:
                family_id = family.split('_')[0]
                if family_id not in family_list:
                        family_list.append(family_id)
                else:
                        continue


        for family_id in family_list:
                if family_id.split('_')[1][0:-2] in ['fa', 'mo']:
                        vepcaller(filepath, family_id)

                else:
                        continue
        
        for monoallelic in monoallelic_files:
                        vepcaller(filepath, monoallelic)

        for merge in full_merge:
                if merge.split('_')[0][0:-2] in ['p1', 's1']:
                        vepcaller(filepath, merge)

if __name__ == "__main__":
        main()
