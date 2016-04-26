#!/usr/bin/python2

import sys
import argparse
from itertools import islice
from os import system
from os.path import splitext, basename, isfile

def get_read_lengths(files):    
    '''Make a dictionary of file names, read lengths, and total reads'''
    file_di = {}
    len_di = {}
    tot = {}
    for file in files:
        base = basename(file)
        tot[base] = 0
        file_di[base] = {}
        with open(file, 'r') as f:
            for l in islice(f, 1, None, 4):
                seq = l.rstrip()
                read_len = len(seq)
                len_di[read_len] = 1
                try:
                    file_di[base][read_len] += 1
                except KeyError:
                    file_di[base][read_len] = 1
                tot[base] += 1
    return file_di, len_di, tot


def write_lenDist_output(file_di, len_di, tot):
    '''Calculate read length ratio and write to output fore each sample'''
    with open("lenDist.tsv", "w") as f:
        for sample in sorted(file_di):
            f.write('\t{0}'.format(sample))
        f.write('\n')
        for some_len in sorted(len_di):
            f.write(str(some_len))
            for sample in sorted(file_di):
                re=0
                if some_len in file_di[sample]:
                    re=file_di[sample][some_len]/float(tot[sample])
                try:
                    f.write('\t{}'.format(str(re)))
                except:
                    print "WARNING: Writing failed for {}".format(sample)
            f.write('\n')


def create_length_dist_image(arg):
    '''Using R, creates a length distribution image for each sample'''
    if arg.image:
        assert isfile("lenDistGraph.R"), \
                "The image producing script (lenDistGraph.R) not present"
        system('Rscript --vanilla lenDistGraph.R lenDist.tsv')
        

def main(arg):
    file_di, len_di, tot = get_read_lengths(arg.inputs)
    write_lenDist_output(file_di, len_di, tot)
    create_length_dist_image(arg)


if __name__ == '__main__':
     parser = argparse.ArgumentParser(
             description='Analyzed the combined blast results')
     parser.add_argument(
             'inputs', 
             action='store', 
             nargs='+',
             help='FILE.fq, trimmed fastq file for each sample ')
     parser.add_argument(
             '--image', 
             action='store_true', 
             help='Indicates a length distribution image be produced')
     arg = parser.parse_args()
     main(arg)
