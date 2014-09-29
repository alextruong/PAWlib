#!/usr/bin/env python


import os
import sys
import commands
import glob


def find_vep():
        #locates the VEP executable
        print 'Locating VEP installation...'
        
        find_veploc = commands.getstatusoutput('find ~/ -name variant_effect_predictor.pl -print 2>/dev/null')
        
        veploc = find_veploc[1]
        
        print 'Done!'
        
        return veploc


def vepcaller(filename, veploc):
        
        #the port 3377 is for hg19
        print 'Running VEP on %s...' % filename
        noextension = filename[0:-4]

        command = veploc + " --port 3337 --cache --refseq --symbol -i " + filename + " -o " + noextension + "_annotated.vcf"

        os.system(command)

        print 'Done!'

        return


def main():

        veploc = find_vep()

        #grabs all of the monoallelic potential variant files, and potential rna-editing files
        monoallelic_files = glob.glob('*_RNA-hom_WES-het.vcf')
        RNA_editing_files = glob.glob('*_RNA-het_WES-hom.vcf')
        full_merge = glob.glob('*_merged.vcf')

        #annotates each file in the lists above
        for monoallelic in monoallelic_files:
                        vepcaller(monoallelic, veploc)
                        
        for RNA_editing in RNA_editing_files:
                        vepcaller(RNA_editing, veploc)

if __name__ == "__main__":
        main()
