#!/bin/bash

scriptdir=$PWD
prs_file=$1
tag_list=$2
iter=$3
root_dir=$4
type=$5


# convert prs to plu
mkdir -p $root_dir/prs_plu_word_$iter
/usr/users/chiaying/AG/scripts/prs_to_plu.py -p $prs_file -t $tag_list -d $root_dir/prs_plu_word_$iter

# build ctl for getting plu word in samples
/data/sls/scratch/chiaying/7lectures/ctl_scripts/build_plu_word_in_samples.py -a `dirname $prs_file` -t $tag_list -n $iter -d $root_dir -o $root_dir/plu_to_wrd_in_samples.ctl

#mkdir $root_dir/ag_word_trans_$iter 
ctl_exec -control $root_dir/plu_to_wrd_in_samples.ctl -overwrite -distributed

# get global word trans 
/data/sls/scratch/chiaying/7lectures/scripts/get_global_word_trans.py -w $tag_list -d $root_dir -b ag_word_trans_$iter -o $root_dir/ag_global_$iter.$type.wrd  

# map plu to words
/data/sls/scratch/chiaying/7lectures/scripts/map_plu_to_wrds.py -t 2 -r $root_dir/global.sec.wrd -p $root_dir/ag_global_$iter.$type.wrd > $root_dir/ag_global_$iter.$type.vocab
