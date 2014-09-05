#!/usr/bin/env python


import os
import sys
import commands
import glob


def find_vep():
        print 'Locating VEP installation...'
        
        find_veploc = commands.getstatusoutput('find ~/ -name variant_effect_predictor.pl -print 2>/dev/null')
        
        veploc = find_veploc[1]
        
        print 'Done!'
        
        return veploc


def vepcaller(filename, veploc):
        
        print 'Running VEP on %s...' % filename
        noextension = filename[0:-4]

        command = veploc + "/ perl variant_effect_predictor.pl --port 3337 --cache --refseq --symbol -i " + filename + " -o " + noextension + "_annotated.vcf"

        os.system(command)

        print 'Done!'

        return


def main():

        veploc = find_vep()

        monoallelic_files = glob.glob('*_RNA-hom_WES-het.vcf')
        full_merge = glob.glob('*_merged.vcf')

        for monoallelic in monoallelic_files:
                        vepcaller(monoallelic, veploc)

        for merge in full_merge:
                if merge.split('_')[0][0:-2] in ['p1', 's1']:
                        vepcaller(merge, veploc)

if __name__ == "__main__":
        main()
