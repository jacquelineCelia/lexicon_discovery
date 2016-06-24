#!/usr/bin/python
import string
import sys
import os
import operator
from optparse import OptionParser


def ParseFile(parse, fn_tags, dirname):
    if not os.path.isdir(dirname):
        os.makedirs(dirname)

    ftags = open(fn_tags)
    fin = open(parse)
    line = fin.readline()
    while line:
        entry = ''
        WordTokens = line.strip().split('(Words')
        for S in WordTokens:
            if S != '':
                SilenceSeprated = S.split('(SIL')
                for W in SilenceSeprated:
                    if W != '':
                        top = '' 
                        btm = '' 
                        tokens = W.split()
                        for idx, t in enumerate(tokens):
                            if "T_" in t:
                                top += t + ' '
                                btm += tokens[idx + 1] + ' '
                                if idx + 2 < len(tokens):
                                    if '[' in tokens[idx + 2]:
                                        btm += tokens[idx + 2] + ' '
                        if top != '' and btm != '':
                            entry += top + ' : ' + btm.strip() + '\n' 

        tag = ftags.readline().strip().split()[0]
        fout = open(dirname + '/' + tag.zfill(4) + '.prs', 'w')
        fout.write(entry)
        fout.close()
        line = fin.readline()
    fin.close()
    ftags.close()

def main():
    ########## parser #######
    parser = OptionParser()
    parser.add_option("-p", dest = "parse", help = "parse file, e.g., ../debug/ag_base/exps/test/wsj.prs.1000")
    parser.add_option("-t", dest = "tags", help = "tag file, e.g., tag.list")
    parser.add_option("-d", dest = "write_dir", help = "write_directory, e.g., /data/sls/scratch/chiaying/AG/prs_to_plu/debug_test/")

    (options, args) = parser.parse_args()
    if len(args) != 0:
       parser.error("incorrect number of arguments, \
            please type -h for more help.")
    
    if not (os.path.isdir(options.write_dir)):
       os.makedirs(options.write_dir)

    ParseFile(options.parse, options.tags, options.write_dir)

if __name__ == "__main__":
    main()
