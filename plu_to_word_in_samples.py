#!/usr/bin/python
import string
import sys
import os
import operator
from optparse import OptionParser

def loadAlign(fname):
    align_start = list()
    align_end = list()
    align_phn = list()
    falign = open(fname)
    line = falign.readline()
    while line:
        if '-' not in line:
            s, e, p = [int(x) for x in line.strip().split()]
            align_start.append(s)
            align_end.append(e)
            align_phn.append(p)
        line = falign.readline()
    falign.close()
    return align_start, align_end, align_phn

def loadTimeStamps(fname):
    start = list()
    end = list()
    fin = open(fname)
    line = fin.readline()
    while line:
        tokens = line.strip().split()
        s = int(tokens[0])
        e = int(tokens[1])
        start.append(s)
        end.append(e)
        line = fin.readline()
    fin.close()
    return start, end

def main():
    ########## parser #######
    parser = OptionParser()
    parser.add_option("-l", dest = "plu", help = "plu file, e.g., /data/sls/scratch/chiaying/7lectures/18.06-1999/L02/prs_plu_800/0000.prs")
    parser.add_option("-a", dest = "align", help = "plu_btm file, e.g., /usr/users/chiaying/AG/sample_rhs/ag_base/exps/18.06.1999/L02/Syl/800/0000")  
    parser.add_option("-b", dest = "bound_phn", help = "bound in frame num file, e.g. /data/sls/scratch/chiaying/7lectures/18.06-1999/L02/bounds_in_phn/0000.phn")
    parser.add_option("-s", dest = "bound_sample", help = "bound in sample file, e.g. /data/sls/scratch/chiaying/7lectures/18.06-1999/L02/bounds_in_samples/0000.lab")
    parser.add_option("-o", dest = "output", help = "words in sample file, e.g., /data/sls/scratch/chiaying/7lectures/18.06-1999/L02/ag_word_trans_800/000.wrd") 

    (options, args) = parser.parse_args()
    if len(args) != 0:
       parser.error("incorrect number of arguments, \
            please type -h for more help.")
    
    if not (os.path.isdir(os.path.dirname(options.output))):
       os.makedirs(os.path.dirname(options.output))

    align_start, align_end, align_phn = loadAlign(options.align) # in frame
    start_frame, end_frame = loadTimeStamps(options.bound_phn)
    start_sample, end_sample = loadTimeStamps(options.bound_sample)

    if not os.path.isdir(os.path.dirname(options.output)):
        os.makedirs(os.path.dirnames(options.output))

    fout = open(options.output, 'w')

    fplu = open(options.plu)
    line = fplu.readline()
    offset = 0
    while line:
        tokens = line.strip().split(':')
        word = '|'.join(tokens[0].strip().split())
        plus_btm = tokens[1].replace("[\\'\\']", '').strip().split()
        if len(plus_btm):
            s = align_start[offset]
            e = align_end[offset + len(plus_btm) - 1]
            s_index = start_frame.index(s)
            e_index = end_frame.index(e)
            fout.write(str(start_sample[s_index]) + ' ' + str(end_sample[e_index]) + ' ' + word + '\n')
            offset += len(plus_btm)
        line = fplu.readline()

    fplu.close()
    fout.close()

if __name__ == "__main__":
    main()
