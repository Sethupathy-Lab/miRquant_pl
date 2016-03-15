#!/usr/bin/python2

import sys
import os
import re

adaptor_begin = 'TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC'
adaptor_end = 'ATCTCGTATGCCGTCTTCTGCTTG'

if sys.argv[1].isdigit():
    adaptor_pos = int(sys.argv[1])
    fastqs = sys.argv[2:]
else:
    fastqs = sys.argv[1:]

for file in fastqs:
    path = os.path.dirname(file) + '/'
    base = os.path.basename(file)
    basename = os.path.splitext(base)[0]
    if 'adaptor_pos' not in locals():
        c = 0
        for bit in basename.split('_'):
            if re.match("^[C,T,G,A]*$", bit):
                    adapt_bit = bit
                    c += 1
        if c != 1:
            print "Could not determine identify adaptor sequence for:", base
            continue
        else:
            adaptor = adaptor_begin + adapt_bit + adaptor_end
            with open('{}{}.adaptor'.format(path, basename), 'w') as f:
                print "Output saved as:", f.name
                f.write('{}\n'.format(adaptor))
    else:
        adapt_bit = basename.split('_')[adaptor_pos]
        if not re.match("^[C,T,G,A]*$", adapt_bit):
            print "Specified adaptor", adapt_bit, "incorrect format for:", base
            continue
        adaptor = adaptor_begin + adapt_bit + adaptor_end
        with open('{}{}.adaptor'.format(path, basename), 'w') as f:
            print "Output saved as:", f.name
            f.write('{}\n'.format(adaptor))

