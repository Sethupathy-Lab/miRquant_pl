#!/usr/bin/python2

###############################################################################
#
# Usage: python generate_mapping_info_table.py PATH_TO_FILE.stats
#    eg: python generate_mapping_info_table.py .../smallRNA/PROJECT/*/*.stats
#
# Output saved as MappingInfoTable.txt
#
###############################################################################

import sys
import os
from os.path import basename, dirname, abspath


def make_tRNA_summary_file():
    '''Gets the # of tRNA mapped reads and writes the output to SAMPLE.stats'''
    for file in sys.argv[1:]:
        dir = dirname(file)
        os.chdir(dir)
        os.system('(head -n 1 TAB_3p_summary.txt && grep -P "chr.*tRNA" TAB_3p_summary.txt | sort -k 6,6nr) > TAB_3p_summary_tRNA.txt')
        tRNAcount = 0
        with open("TAB_3p_summary_tRNA.txt", 'r') as f:
            f.readline()
            for l in f:
                l = l.rstrip().split('\t')
                tRNAcount += float(l[5])
        name = dir.split('/')[-1]
        with open("{}stats".format(name), 'a') as f:
            f.write('tRNAMapped: {}\n'.format(str(int(tRNAcount))))


def mapping_stats_dict():
    '''For each sample input, opens the SAMPLE./SAMPLE.stats file and imports
    the contents into an output dictionary.'''
    out_di = {}
    for file in sys.argv[1:]:
        dir = dirname(file)
        name = dir.split('/')[-1]
        out_di[name] = {}
        out_di[name]["File"] = file
        with open(file, 'r') as fi:
            for l in fi:
                a, b = l.rstrip().split(":")
                if "." in b:
                    b = b.split(".")[0]
                out_di[name][a] = b
    return out_di


def calculate_additional_stats(out_di):
    '''For each sample input, calculates trimming and mapping efficencies.
    Loads the calculated efficency into the output dictionary'''
    for name in out_di.keys():
        tot_r = float(out_di[name]["TotReads"])
        trim_r = float(out_di[name]["TrimmReads"])
        short_r = float(out_di[name]["ShortReads"])
        EM_r = float(out_di[name]["EMhits"])
        Emiss_r = float(out_di[name]["EMmiss"])
        map_r = float(out_di[name]["Mapped"])
        miRmap_r = float(out_di[name]["miRMapped"])
        tRNAmap_r = float(out_di[name]["tRNAMapped"])
        out_di[name]["TrimReadPer"] = "{0:.2f}".format(trim_r / tot_r * 100)
        out_di[name]["ShortReadPer"] = "{0:.2f}".format(short_r / tot_r * 100)
        out_di[name]["EMReadPer"] = "{0:.2f}".format(EM_r / trim_r * 100)
        out_di[name]["EMissReadPer"] = "{0:.2f}".format(Emiss_r / trim_r * 100)
        out_di[name]["MapReadPer"] = "{0:.2f}".format(map_r / trim_r * 100)
        out_di[name]["miRMapPer"] = "{0:.2f}".format(miRmap_r / map_r * 100)
        out_di[name]["tRNAMapPer"] = "{0:.2f}".format(tRNAmap_r / map_r * 100)
    return out_di


def output_line_headers():
    '''Make list of list of out_di keys and line header values'''
    return [["File","File"], 
            ["TotReads","Total Reads"], 
            ["TrimmReads","Trimmed Reads"],
            ["TrimReadPer","Percent Trimmed Reads"],
            ["ShortReads","Too Short Reads"],
            ["ShortReadPer","Percent Too Short"],
            ["EMhits","Exact Match Reads"],
            ["EMReadPer","Percent Exact Matches"],
            ["EMmiss","Mismatch Reads"],
            ["EMissReadPer","Percent Mismatched"],
            ["Mapped","Mapped Reads"],
            ["MapReadPer","Percent Mapped"],
            ["miRMapped","miR Mapped Reads"],
            ["miRMapPer","Percent miR Mapped"],
            ["tRNAMapped","tRNA Mapped Reads"],
            ["tRNAMapPer","Percent tRNA Mapped"]]
        
        
def write_mapping_file(out_di, out_dir, line_head_li):
    '''Writes the output dictionary to a file, called MappingInfoTable.csv'''
    with open('{}/MappingInfoTable.csv'.format(out_dir), 'w') as f:
        f.write('Sample_name,{}\n'.format(','.join(sorted(out_di))))
        for item in line_head_li:
            f.write('{}'.format(item[1]))
            for sample in sorted(out_di):
                f.write(',{}'.format(out_di[sample][item[0]]))
            f.write('\n')


def main():
    out_dir = os.getcwd()
    make_tRNA_summary_file()
    out_di = mapping_stats_dict()
    out_di = calculate_additional_stats(out_di)
    line_head_li = output_line_headers()
    write_mapping_file(out_di, out_dir, line_head_li)


if __name__ == '__main__':
    main()
