#!/usr/users/chiaying/install/epd-7.2-2-rh5-x86_64/bin//python2.7
import string
import sys
import os
import operator
import collections

from optparse import OptionParser
from operator import itemgetter

def loadAlign(fname):
    start = list()
    end = list()
    word = list()
    fin = open(fname)
    line = fin.readline()
    while line:
        tokens = line.strip().split()
        s = float(tokens[0])
        e = float(tokens[1])
        w = ''
        if len(tokens) > 2:
            w = ' '.join(tokens[2:len(tokens)])
        start.append(s)
        end.append(e)
        word.append(w)
        line = fin.readline()
    fin.close()
    return start, end, word

def findClosest(s, ref_start):
    diff = [abs(x - s) for x in ref_start]
    return min(enumerate(diff), key=itemgetter(1))[0] 


def main():
    ########## parser #######
    parser = OptionParser()
    parser.add_option("-r", dest = "reference", help = "e.g., 18.06-1999/L02/global.sec.wrd")
    parser.add_option("-p", dest = "plu", help = "e.g.,  18.06-1999/L02/ag_global_200.fussy_match_Syl_dphmm.wrd")  
    parser.add_option("-t", dest = "threshold", help = "e.g, threshold on the cluster szie for pringint out", default = 0, type = int)

    (options, args) = parser.parse_args()
    if len(args) != 0:
       parser.error("incorrect number of arguments, \
            please type -h for more help.")

    ref_start, ref_end, ref_word = loadAlign(options.reference)
    hyp_start, hyp_end, hyp_word = loadAlign(options.plu)

    plu_dict = dict()

    for index, w in enumerate(hyp_word):
        if w:
            s = hyp_start[index]
            e = hyp_end[index]
            nearest_index_in_ref = findClosest(s, ref_start)
            word = ref_word[nearest_index_in_ref] + ' '
            if e > ref_end[nearest_index_in_ref]:
                nearest_index_in_ref += 1 
                while nearest_index_in_ref < len(ref_end) and e > ref_end[nearest_index_in_ref]:
                    word += ref_word[nearest_index_in_ref] + ' '
                    nearest_index_in_ref += 1
                if nearest_index_in_ref < len(ref_start):
                    if (e - ref_start[nearest_index_in_ref]) / (ref_end[nearest_index_in_ref] - ref_start[nearest_index_in_ref]) >= 0.5:
                        word += ref_word[nearest_index_in_ref]
            if w not in plu_dict:
                plu_dict[w] = list()
            plu_dict[w].append(word.strip())
            
    summary = dict()
    for w in plu_dict:
        if len(plu_dict[w]) >= options.threshold:
            summary[w] = collections.Counter(plu_dict[w]) 
            print len(plu_dict[w]), w, summary[w].most_common(len(summary[w]))
            # print len(plu_dict[w]), w, summary[w]

if __name__ == "__main__":
    main()

