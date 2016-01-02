/* -*- C++ -*-
 *
 * Copyright (c) 2016
 * Spoken Language Systems Group
 * MIT Computer Science and Artificial Intelligence Laboratory
 * Massachusetts Institute of Technology
 *
 * FILE: datum.h 
 *   				      				                             
 * Chia-ying (Jackie) Lee <chiaying@csail.mit.edu>				    
 * Jan 2016 
*********************************************************************/

#include "datum.h"

#include <fstream>
#include <iostream>
#include <sstream>

using namespace std;

Datum::Datum() {
}

Datum::Datum(Config* config) {
}

vector<symbol> Datum::retrieve_symbols() {
    vector<symbol> terminals;
    vector<Segment*>::iterator iter = _segments.begin(); 
    for (; iter != _segments.end(); ++iter) {
        terminals.push_back((*iter) -> id_symbol());
    }
    return terminals;
}

void Datum::print_terminals() {
    for (size_t i = 0; i < _segments.size(); ++i) {
        cout << _segments[i] -> id() << " ";
    }
    cout << endl;
}

void Datum::ClearSegs() {
    vector<Segment*>::iterator s_iter = _segments.begin();
    for (; s_iter != _segments.end(); ++s_iter) {
        delete *s_iter;
    }
    _segments.clear();
}

void Datum::Save(const string& dir) {
    string filename = dir + _tag;
    // cout << filename << endl;
    ofstream fout(filename.c_str());
    vector<Segment*>::iterator s_iter = _segments.begin();
    int total_frame_num = 0;
    for (; s_iter != _segments.end(); ++s_iter) {
        if (((*s_iter) -> id_symbol()).string_reference() != "\'\'") {
            fout << total_frame_num << " " << \
                total_frame_num + (*s_iter) -> frame_num() - 1 << " " << (*s_iter) -> id_symbol() << endl; 
            total_frame_num += (*s_iter) -> frame_num();
        }
        else {
            fout << "- - " << (*s_iter) -> id_symbol() << endl;
        }
    }
    fout.close();
}

Datum::~Datum() {
    vector<Bound*>::iterator iter = _bounds.begin();
    for (; iter != _bounds.end(); ++iter) {
        delete *iter;
    }
    vector<Segment*>::iterator s_iter = _segments.begin();
    for (; s_iter != _segments.end() ; ++s_iter) {
        delete *s_iter;
    }
}

