#!/bin/bash
mkdir Syl_dphmm
/usr/users/chiaying/AG/fussy_match/ag_base/adaptor \
    -A Syl_dphmm/MIT-World-2005-OCW-Friedman-Flat-16May2005.prs \
    -C -D -E -F Syl_dphmm/MIT-World-2005-OCW-Friedman-Flat-16May2005.trace \
    -G Syl_dphmm/MIT-World-2005-OCW-Friedman-Flat-16May2005.wlt \
    -R -1 \
    -n 2000 \
    -r 0 \
    -a 1e-4 \
    -b 1e4 \
    -e 1 \
    -f 1 \
    -g 10 \
    -h 0.1 \
    -w 1 \
    -T 1 \
    -m 0 \
    -d 0 \
    -x 10 \
    -k MIT-World-2005-OCW-Friedman-Flat-16May2005.dphmm.input \
    -l configuration_dphmm \
    -i snapshot_dphmm \
    -j dphmm \
    -o model_id_dphmm \
    -u gaussian_seed_dphmm \
    -v Syl_dphmm \
Syl_dphmm.lt < phn.txt
