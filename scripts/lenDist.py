#!/usr/bin/python2

import sys
from itertools import islice
from os import system
from os.path import splitext, basename

file_di = {}
len_di = {}
tot = {}
for file in sys.argv[1:]:
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
                print "why"
                print re
                print some_len
                print sample
                print
        f.write('\n')

system('Rscript --vanilla lenDistGraph.R lenDist.tsv')
        



