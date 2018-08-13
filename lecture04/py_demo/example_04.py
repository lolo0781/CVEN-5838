#! /usr/bin/env python

import sys
import getopt

def myExample(argv):
    
    print 'argv = ', argv
    InputFileName = 'NULL'

    # command line options and arguments 
    try:
        opts, args = getopt.getopt(argv, 'h i:')

    except:
        print 'Error in user input options'
        exit(0)

    for opt, arg in opts:
        print 'opt = ', opt
        print 'arg = ', arg

        if opt == '-h':
            print 'No help available'
            exit(0)

        if opt == '-i':
            InputFileName = arg
            print 'Setting input file name to ', InputFileName


    # Opening and reading and writing files
    try:
        f = open(InputFileName, 'r')
        g = open(InputFileName + '.out', 'w'

        for line in f:
            print line

        f.close()
        g.close()

    except:
        print 'There was a problem opening the input file', InputFileName
    

if __name__=="__main__":

    print '*** Executing Main Routine ***'
    
    # store all entries in the list from 1 on
    args = sys.argv[1:]

    myExample(args)
    
