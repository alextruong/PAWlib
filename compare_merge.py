#!/usr/bin/env python

from __future__ import division
import os
import sys
import commands
import glob


def makefolder():
        """Make working directory for script outputs"""

        newdir = True;

        while newdir == True:
                outdir = str(raw_input('Please name output file directory (to be created here; please ensure it is empty for best results): '))

                if os.path.isdir(outdir) == True:
                        print 'It looks like this directory already exists.'

                        ifdelete = None

                        while ifdelete not in ('y', 'n'):
                                ifdelete = str(raw_input('Do you want to delete the existing directory and write an empty instance (y/n)? '))
                                if ifdelete == 'y':
                                        print 'Deleting existing directory...'
                                        os.system('rm -rf %s' % outdir)
                                        print 'Making local outfile directory in same location as script execution called %s...' % outdir
                                        os.system('mkdir %s' % outdir)

                                        newdir = False

                                elif ifdelete == 'n':
                                        print 'Please type a different directory name.'
                                else:
                                        print 'please choose a valid option.'

                else:
                        print 'Making local outfile directory in same location as script execution called %s...' % outdir
                        os.system('mkdir %s' % outdir)

                        newdir = False

        print 'Done!'
        print '-----------------------------------------------'

        return outdir
