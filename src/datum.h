/* -*- C++ -*-
 *
 * Copyright (c) 2014
 * Spoken Language Systems Group
 * MIT Computer Science and Artificial Intelligence Laboratory
 * Massachusetts Institute of Technology
 *
 * FILE: datum.h 
 *                    
 * Chia-ying (Jackie) Lee <chiaying@csail.mit.edu>
 * Jan 2016 
*********************************************************************/

#ifndef DATUM_H
#define DATUM_H

#include <vector>

#include "bound.h"
#include "config.h"
#include "segment.h"
#include "sym.h"

class Datum {
 public:
     Datum(Config*);
     Datum();
     std::vector<symbol> retrieve_symbols();
     void set_tag(const std::string& tag) {_tag = tag;}
     void set_bounds(std::vector<Bound*> bounds) {_bounds = bounds;}
     void set_frame_num(int n) {_total_frame_num = n;}
     void set_corrupted(bool corrupted) {_corrupted = corrupted;}
     void set_corrupted_type(int type) {_corruption_type = type;}
     void set_segments(vector<Segment*>& segs) {_segments = segs;}
     void print_terminals();
     void ClearSegs(); 
     void Save(const std::string&);
     int get_frame_num() {return _total_frame_num;}
     vector<Bound*> bounds() {return _bounds;}
     vector<Segment*>& segments() {return _segments;} 
     std::string tag() {return _tag;}
     ~Datum();
 private:
     std::vector<Segment*> _segments;
     std::vector<Bound*> _bounds;
     std::string _tag;
     int _total_frame_num;
     bool _corrupted;
     int _corruption_type;
};

#endif
