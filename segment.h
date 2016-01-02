/* -*- C++ -*-
 *
 * Copyright (c) 2014
 * Spoken Language Systems Group
 * MIT Computer Science and Artificial Intelligence Laboratory
 * Massachusetts Institute of Technology
 *
 * All Rights Reserved

 *	FILE: segment.h 
 *										                            *
 *   				      				                            * 
 *   Chia-ying (Jackie) Lee <chiaying@csail.mit.edu>				*
 *   Feb 2014							                            *
*********************************************************************/

#ifndef SEGMENT_H
#define SEGMENT_H

#include "bound.h"
#include "sym.h"

#include <vector>

using namespace std;

class Segment {
 public:
     Segment(symbol symbol) {_id_symbol = symbol;}
     Segment(int, int);
     int id() {return _id;}
     symbol id_symbol() {return _id_symbol;}
     void push_back(Bound*);
     vector<float*> data() {return _data;}
     int frame_num() {return _frame_num;}
     float* frame(int index) {return _data[index];}
     void set_state_seq(vector<int>& state_seq) {_state_seq = state_seq;}
     void set_mix_seq(vector<int>& mix_seq) {_mix_seq = mix_seq;}
     vector<int> state_seq() {return _state_seq;}
     vector<int> mix_seq() {return _mix_seq;} 
     ~Segment();
 private:
     symbol _id_symbol;
     int _id;
     int _frame_num;
     vector<float*> _data;
     vector<int> _state_seq;
     vector<int> _mix_seq;
};

#endif
