#!/usr/bin/python
from optparse import OptionParser

def loadList(wavelist):
    tags = list()
    waves = list()
    fin = open(wavelist)
    line = fin.readline()
    while line:
        tag, wav = line.strip().split()
        tags.append(tag)
        waves.append(wav)
        line = fin.readline()
    fin.close()
    return tags, waves

def main():
    ########## parser #######
    parser = OptionParser()
    parser.add_option("-t", dest = "wavelist", help = "e.g., ../18.06-1999/L02/speech_detect.list")
    parser.add_option("-n", dest = "iteration", help = "number of iterations of the ag learning")
    parser.add_option("-d", dest = "root_dir", help = "output dir, e.g., /data/sls/scratch/chiaying/7lectures/18.06-1999/L02/")
    parser.add_option("-o", dest = "output", help = "output file, e.g., ../18.06-1999/L02/plu_to_wrd_in_sample.ctl")
    parser.add_option("-y", dest = "type", help = "syl or word, e.g., word", default = "word")
    parser.add_option("-b", dest = "subdir", help = "sub dir, e.g., sample_rhs", default = "sample_rhs")
    parser.add_option("-a", dest = "acoustic", help = "acoustic dir, e.g, /usr/users/chiaying/AG/18.06-1999/L02/Syl/")

    (options, args) = parser.parse_args()
    if len(args) != 0:
       parser.error("incorrect number of arguments, \
            please type -h for more help.")

    tokens = options.wavelist.split('/')
    lecture = tokens[-2]
    course = tokens[-3]

    print lecture, course

    print options.type
    print options.acoustic
    print options.iteration


    tags, waves = loadList(options.wavelist)
    fout = open(options.output, 'w')
    for idx, tag in enumerate(tags):
        entry = '/usr/users/chiaying/AG/scripts/plu_to_word_in_samples.py -l ' + options.root_dir + '/prs_plu_' + options.type + '_' + options.iteration + '/' + tag.zfill(4) + '.prs '\
                + ' -a ' + options.acoustic + '/' + options.iteration + '/' + tag.zfill(4) + ' -b ' + options.root_dir + '/bounds_in_phn/' + tag.zfill(4) + '.phn -s ' + options.root_dir + '/bounds_in_samples/' + tag.zfill(4) + '.lab -o ' + options.root_dir + '/ag_' + options.type + '_trans_' + options.iteration + '/' + tag.zfill(4) + '.wrd\n' 
        entry += '-------------------\n'
        fout.write(entry)
    fout.close()

if __name__ == "__main__":
    main()
