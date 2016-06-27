#!/usr/bin/python
import string
import sys
import os
import operator
from optparse import OptionParser

def loadList(wavelist):
    tags = list()
    waves = list()
    fin = open(wavelist)
    line = fin.readline()
    while line:
        tag, wav = line.strip().split()
        tags.append(tag.zfill(4))
        waves.append(wav)
        line = fin.readline()
    fin.close()
    return tags, waves

def parseWaveTime(waves):
    start_sample = list()
    end_sample = list()
    for w in waves:
        tokens = w.split('::')
        s, e = [int(t) * 16 for t in tokens[-1].strip().split(':')]
        start_sample.append(s)
        end_sample.append(e)
    return start_sample, end_sample

def main():
    ########## parser #######
    parser = OptionParser()
    parser.add_option("-w", dest = "wavelist", help = "e.g., /data/sls/scratch/chiaying/7lectures/18.06-1999/L02/speech_detect.list")
    parser.add_option("-d", dest = "root_dir", help = "root_dir dir, e.g., /data/sls/scratch/chiaying/7lectures/18.06-1999/L02/")  
    parser.add_option("-b", dest = "sub_dir", help = "sub_dir, e.g, ag_word_trans or word_trans")
    parser.add_option("-o", dest = "output", help = "output file that contains the global word align in samples, e.g. /data/sls/scratch/chiaying/7lectures/18.06-1999/L02/global.wrd")
    parser.add_option("-e", dest = "ext", help = "extension, e.g., wrd", default = "wrd")

    (options, args) = parser.parse_args()
    if len(args) != 0:
       parser.error("incorrect number of arguments, \
            please type -h for more help.")

    tags, waves = loadList(options.wavelist)
    start_sample, end_sample = parseWaveTime(waves)

    fout = open(options.output, 'w')

    last_sample = 0
    for index, t in enumerate(tags):
        if last_sample < start_sample[index]:
            fout.write(str(last_sample / 16000.0) + ' ' + str(start_sample[index] / 16000) + '\n')
        last_sample = end_sample[index]
        fwrd = open(options.root_dir + '/' + options.sub_dir + '/' + tags[index].zfill(4) + '.' + options.ext)
        line = fwrd.readline()
        while line:
            tokens = line.strip().split()
            s = int(tokens[0]) 
            e = int(tokens[1])
            w = ''
            if len(tokens) > 2 :
                w = tokens[2]
            fout.write(str((s + start_sample[index]) / 16000.0) + ' ' + str((e + start_sample[index]) / 16000.0) + ' ' + w + '\n')
            line = fwrd.readline()
        fwrd.close()
    fout.close()

if __name__ == "__main__":
    main()

