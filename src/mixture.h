/* -*- C++ -*-
 *
 * Copyright (c) 2014
 * Spoken Language Systems Group
 * MIT Computer Science and Artificial Intelligence Laboratory
 * Massachusetts Institute of Technology
 *
 * FILE: mixture.h 
 *                    
 * Chia-ying (Jackie) Lee <chiaying@csail.mit.edu>
 * Jan 2016 
*********************************************************************/
#ifndef MIXTURE_H
#define MIXTURE_H

#include <iostream>
#include <fstream>
#include <vector>

#include "config.h"

using namespace std;

class Mixture {
 public:
     Mixture(Config*);
     void set_det(float det) {_det = det;}
     void set_det();
     void set_mean(vector<float> mean) {_mean = mean;}
     void set_pre(vector<float> pre) {_pre = pre;}
     void Minus(float*);
     void Plus(float* data);
     void print(); 
     vector<float> mean() {return _mean;}
     vector<float> pre() {return _pre;}
     float det() {return _det;}
     float likelihood(float*);
     ~Mixture() {};
 private:
     Config* _config;
     float _det;
     int _dim;
     vector<float> _mean;
     vector<float> _pre;
};

#endif
