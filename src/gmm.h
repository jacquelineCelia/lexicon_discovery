/* -*- C++ -*-
 *
 * Copyright (c) 2014
 * Spoken Language Systems Group
 * MIT Computer Science and Artificial Intelligence Laboratory
 * Massachusetts Institute of Technology
 *
 * FILE: gmm.h 
 *                    
 * Chia-ying (Jackie) Lee <chiaying@csail.mit.edu>
 * Jan 2016 
*********************************************************************/
#ifndef GMM_H
#define GMM_H

#include "config.h"
#include "mixture.h"

#include <fstream>
#include <iostream>
#include <vector>

using namespace std;

class GMM {
 public:
     GMM(Config*);
     GMM(Config*, int);
     void load(ifstream&);
     void set_mixture_det(int, float);
     void set_mixture_mean(int, vector<float>&);
     void set_mixture_pre(int, vector<float>&);
     void set_weight(vector<float>& weight) {_weight = weight;}
     void ComputeLikehood(vector<float*>, float*);
     void print();
     float ComputeLikehood(float* data);
     Mixture& mixture(int index) {return _mixtures[index];}
     vector<float> weight() {return _weight;}
     vector<float> ComponentLikelihood(float*);
     void Minus(float*, int);
     void Plus(float* data, int index);
     void Save(ofstream&);
     int mix_num() {return _mix_num;}
     ~GMM();
 private:
     Config* _config;
     vector<float> _weight;
     vector<Mixture> _mixtures;
     int _mix_num;
};

#endif
