#!/usr/bin/python2

###############################################################################
#
# Usage: python process_all_summary2tab.py SCRIPTDIR SPEC FULL_PATH2FILES
#
#  Where: SCRIPTDIR = Directory where the perl script summary2Tab_clust.pl is
#         SPEC = The species this is for (hsa or mmu)
#         FULL_PATH2FILES = Full path to any unique file in g1Results directory
#          -eg: /proj/seth_lab/users/USER/smrnapipeline/SAMPLE/(cont.)
#                  IntermediateFiles/g1Results/shift_summary.txt
#
###############################################################################

import sys
import os
import subprocess as sp


def run_summary2Tab_clust():
    sum2tab = '{}/summary2Tab_clust.pl'.format(sys.argv[1])
    spec = sys.argv[2]
    for file in sys.argv[3:]:
        dir = os.path.dirname(file)
        print dir
        os.chdir(dir)
        os.system('perl {} lenDist_summary.txt 0 {}'.format(sum2tab, spec))
        os.system('perl {} 3p_summary.txt 0 {}'.format(sum2tab, spec))
        os.system('perl {} ed_summary.txt 0 {}'.format(sum2tab, spec))
        os.system('(head -n 1 TAB_3p_summary.txt && grep {} TAB_3p_summary.txt \
                | sort -k 6,6nr) > TAB_3p_summary_miR.txt'.format(spec))
        proc = sp.Popen([('head -n 2 TAB_lenDist_summary.txt | tail -n 1')],
                stdout = sp.PIPE, shell=True)
        countLine = proc.communicate()[0]
        countParts = countLine.split('\t')
        count = float(countParts[5])
        miRcount = 0
        with open('TAB_3p_summary_miR.txt', 'r') as f:
            for l in f:
                l = l.strip('\n')
                try:
                    miRcount += float(l.split('\t')[5])
                except:
                    pass
        os.system('mv TAB_3p_summary_miR.txt TAB_lenDist_summary.txt \
                TAB_3p_summary.txt TAB_ed_summary.txt Shrimp_results.bed ../../')
        dirParts = dir.split('/')
        dirParts.pop()
        dirParts.pop()
        topDir='/'.join(dirParts)
        Name = dirParts.pop()[:-1]
        os.chdir(topDir)
        os.system('tar -zcvf {}.tgz IntermediateFiles'.format(Name))
        os.system('echo "Mapped: {}" >> {}.stats'.format(count, Name))
        os.system('echo "miRMapped: {}" >> {}.stats'.format(miRcount, Name))
        

def write_summary_table():
    Table, TKeys = {}, []
    for f in sys.argv[3:]:
        dir = os.path.dirname(f)
        dirParts = dir.strip('\n').split('/')
        dirParts.pop()
        dirParts.pop()
        topDir = '/'.join(dirParts)
        Name = dirParts.pop()[:-1]
        projDir = '/'.join(dirParts)
        fi = '{}/{}.stats'.format(topDir, Name)
        Table[Name] = {}
        with open(fi, 'r') as f:
            for l in f:
                A, B = l.strip('\n').split(':')
                Table[Name][A] = B
                TKeys.append(A)
    for n in Table:
        print '\t{}'.format(n)
    for k in TKeys:
        print k, 
        for n in Table:
            print '\t{}'.format(Table[n][k])


def main():
    run_summary2Tab_clust()
    write_summary_table()


if __name__ == '__main__':
    main()
