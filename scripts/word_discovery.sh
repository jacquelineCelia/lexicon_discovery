#!/bin/bash

if [ "$#" -ne 5 ]; then
    echo "./word_discovery.sh _parse_file_ _tag_list_ _iter_num_ _root_dir_" 
    echo "_parse_file_ the parse file"
    echo "_tag_list tags of the sentences" 
    echo "_iter_num_ the iteration number of the parse file"
    echo "_root_dir_ root dir for the output files"
    exit;
fi

prs_file=$1
tag_list=$2
iter=$3
root_dir=$4
type="syl_dphmm"

# convert prs to plu
mkdir -p $root_dir/prs_plu_word_$iter
./prs_to_plu.py -p $prs_file -t $tag_list -d $root_dir/prs_plu_word_$iter

# build ctl for getting plu word in samples
./build_plu_word_in_samples.py -a `dirname $prs_file` -t $tag_list -n $iter -d $root_dir -o $root_dir/plu_to_wrd_in_samples.ctl

mkdir $root_dir/ag_word_trans_$iter 
ctl_exec -control $root_dir/plu_to_wrd_in_samples.ctl -overwrite -distributed

# get global word trans 
./get_global_word_trans.py -w $tag_list -d $root_dir -b ag_word_trans_$iter -o $root_dir/ag_global_$iter.$type.wrd  

# map plu to words
./map_plu_to_wrds.py -t 2 -r $root_dir/global.sec.wrd -p $root_dir/ag_global_$iter.$type.wrd > $root_dir/ag_global_$iter.$type.vocab
