#!/usr/bin/python2

###############################################################################
#
# This script calculates normalizes the data across samples by calculating the 
# reads assigned per million mapped reads.
#
# Usage: python script.py species path_to_files
#
#   species: species abbreviation (mmu, hsa, rno, cast)
#   path_to_files: full path to TAB_lenDist_summary.txt
#
# Outputs:
#   RPMM_all.tsv: RPMM for all contigs
#   RPMM_miRs_only.tsv: RPMM for only the miRs
#   RPMM_miRs_over_100.tsv: RPMM for miRs in which RPMM was over 100 for at
#                           least one sample
#
###############################################################################

import sys

def check_species_input():
    '''This checks that a the species input makes sense.  If species input
    known to be correct, but not in list below, add species abbreviation'''
    if sys.argv[1] in ['mmu', 'hsa', 'rno', 'cast']:
        sys.stdout = open('RPM_all.tsv', 'w')
        return sys.argv[1]
    else:
        print '\nSpecies not recognized, known species: mmu, hsa, rno, cast'
        print 'If species abbreviation known to be correct, edit script.\n'
        print 'Usage: python genNormalRPMM.py species path_to_files > output_name\n'
        sys.exit()


def get_data_from_file():
    '''Brings in data and counts from TAB_lenDist file'''
    total_counts, datout, li = {}, {}, {}
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
                datout[file][l[0]] = l[count_i]
                li[l[0]] = 1
    return datout, li, total_counts
    

def calculate_RPMM_and_seperate_miRs(datout, li, total_counts, species):
    '''Calculates the RPMM for each line.  Separates the miRs from all parts,
    and makes a list of miRs that are over a threshold (default 100) for at
    least one of the samples'''
    thresh = 100
    sample_miRs = {}        
    miRs_over_thresh = []
    sorted_files = sorted(datout)
    for file in sorted_files:
        sample_miRs[file] = {}
        print "\t{}".format(file),
    print
    c = 0
    for miR in li:
        print miR,
        for file in sorted_files:
            re = 0
            if miR in datout[file]:
                re = 1000000*(float(datout[file][miR])/total_counts[file])
                if species in miR:
                    c += 1
                    sample_miRs[file][miR] = re
                    if re >= thresh:
                        miRs_over_thresh.append(miR)
            print "\t{}".format(re),
        print
    miRs_over_thresh = set(miRs_over_thresh)
    return sample_miRs, miRs_over_thresh, thresh


def write_miRs_only_RPMM(sample_miRs):
    '''Writes RPMM file for only the miRs to an output file, called
    called RPMM_miRs_only.tsv'''
    all_miRs = set([miR for key in sample_miRs for miR in sample_miRs[key]])
    sample_list = sorted(sample_miRs.keys())
    with open('RPM_miRs_only.tsv', 'w') as f:
        f.write('\t{}\n'.format('\t'.join(sample_list)))
        for miR in all_miRs:
            f.write(miR)
            for sample in sample_list:
                try:
                    f.write('\t{}'.format(sample_miRs[sample][miR]))
                except KeyError:
                    f.write('\t{}'.format('0'))
            f.write('\n')
    

def write_miRs_over_100(sample_miRs, miRs_over_thresh, thresh):
    '''Writes miRs where RPMM is greater than threshold for at least one sample
    to an output file (called RPMM_miRs_over_(threshold).tsv'''
    sample_list = sorted(sample_miRs.keys())
    with open('RPM_miRs_over_{}.tsv'.format(str(thresh)), 'w') as f:
        f.write('\t{}\n'.format('\t'.join(sample_list)))
        for miR in miRs_over_thresh:
            f.write(miR)
            for sample in sample_list:
                try:
                    f.write('\t{}'.format(sample_miRs[sample][miR]))
                except KeyError:
                    f.write('\t{}'.format('0'))
            f.write('\n')
    

def main():
    species = check_species_input()
    datout, li, total_counts = get_data_from_file()
    sample_miRs, miRs_over_thresh, thresh = calculate_RPMM_and_seperate_miRs(
            datout, li, total_counts, species)
    write_miRs_only_RPMM(sample_miRs)
    write_miRs_over_100(sample_miRs, miRs_over_thresh, thresh)

if __name__ == '__main__':
    main()
