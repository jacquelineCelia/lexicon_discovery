/* -*- C++ -*-
 *
 * Copyright (c) 2014
 * Spoken Language Systems Group
 * MIT Computer Science and Artificial Intelligence Laboratory
 * Massachusetts Institute of Technology
 *
 * FILE: mixture.cc
 *                    
 * Chia-ying (Jackie) Lee <chiaying@csail.mit.edu>
 * Jan 2016 
*********************************************************************/

#include <cmath>
#include "mixture.h"

Mixture::Mixture(Config* config) {
    _config = config;
    _dim = _config -> dim();
    _mean.resize(_dim, 0);
    _pre.resize(_dim, 0);
}

void Mixture::set_det() {
    _det = 0;
    for (int i = 0; i < _config -> dim(); ++i) {
        _det += log(_pre[i]);
    }
    _det *= 0.5;
    _det -= 0.5 * (_config -> dim()) * 1.83787622175;
    // log(2*3.1415926) = 1.83787622175
}

void Mixture::Minus(float* data) {
    for (int i = 0; i < _config -> dim(); ++i) {
        _mean[i] -= data[i];
        _pre[i] -= data[i] * data[i];
    }
}

void Mixture::Plus(float* data) {
    for(int i = 0 ; i < _config -> dim(); ++i) {
        _mean[i] += data[i];
        _pre[i] += data[i] * data[i];
    }
}

void Mixture::print() {
    for(int i = 0 ; i < _config -> dim(); ++i) {
        cout << "mean[" << i << "]: " << _mean[i] << " "
            << "pre[" << i << "]: " << _pre[i] << endl;
    }
}

float Mixture::likelihood(float* data) {
    float likelihood = 0;
    for (int i = 0; i < _config -> dim(); ++i) {
        likelihood += (data[i] - _mean[i]) * (data[i] - _mean[i]) * _pre[i];
    }
    likelihood *= -0.5;
    return _det + likelihood;
}
