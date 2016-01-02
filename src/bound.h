/* -*- C++ -*-
 *
 * Copyright (c) 2016
 * Spoken Language Systems Group
 * MIT Computer Science and Artificial Intelligence Laboratory
 * Massachusetts Institute of Technology
 *
 * FILE: bound.h
 *   				      				                            
 * Chia-ying (Jackie) Lee <chiaying@csail.mit.edu>				    
 * Jan 2016 
*********************************************************************/

#ifndef BOUND_H
#define BOUND_H

#include <vector>
#include <cstring>

using namespace std;

class Bound {
 public:
     Bound(int, int, int);
     void set_data(float**, int);
     int get_frame_num() {return _end_frame - _start_frame + 1;}
     int start_frame() {return _start_frame;}
     int end_frame() {return _end_frame;}
     vector<float*> data() {return _data_ptr;}
     ~Bound();
 private:
     float** _data;
     vector<float*> _data_ptr;
     int _dim;
     int _start_frame;
     int _end_frame;
};

#endif 
