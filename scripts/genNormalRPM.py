#!/usr/bin/python2

###############################################################################
#
# This script calculates normalizes the data across samples by calculating the 
# reads assigned per million mapped reads.
#
# Usage: python script.py species path_to_files > path_to_output
#
#   species: species abbreviation (mmu, hsa, ect.) to get only miRs
#            set 'all' to get RPMM for all reads
#
###############################################################################

import sys

if sys.argv[1] in ['mmu', 'hsa', 'rno', 'cast', 'all']:
    species = sys.argv[1]
else:
    print 'species not recognized, check script to see species and append list'
    sys.exit()

total_counts, datout, list = {}, {}, {}
for file in sys.argv[2:]:
    with open(file, 'r') as f:
        header = f.readline().split('\t')
        count_i = header.index('Count')
        tot = f.readline().split('\t')[count_i]
        total_counts[file] = 0
        datout[file] = {}
        for l in f:
            l = l.split('\t')
            total_counts[file] += float(l[count_i])
            if species == 'all' or species in l[0]:
                datout[file][l[0]] = l[count_i]
                list[l[0]] = 1

thresh = 100
miR_over_thresh = {}        
for file in datout:
    miR_over_thresh[file] = []
    print "\t{}".format(file),
print
for miR in list:
    print miR,
    for file in datout:
        re = 0
        if miR in datout[file]:
            re = 1000000*(float(datout[file][miR])/total_counts[file])
            if re >= thresh and 'miR' in miR:
                miR_over_thresh[file].append([miR, re])
        print "\t{}".format(re),
    print


with open('miRs_over_RMM_{}.tsv'.format(str(thresh)), 'w') as f:
    max_len = max([len(v) for v in miR_over_thresh.values()])
    high_RRM_miRs = [[[k, '']] + 
            [['miR', 'RRM']] + 
            v + [['', '']] * (max_len - len(v)) 
            for k, v in miR_over_thresh.iteritems()]
    high_RRM_miRs_out = zip(*high_RRM_miRs)
    for line in high_RRM_miRs_out:
        tmp = []
        for item in line:
            tmp += item
        f.write('{}\n'.format('\t'.join(map(str, tmp))))
