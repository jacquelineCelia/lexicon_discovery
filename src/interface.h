/* -*- C++ -*-
 *
 * Copyright (c) 2014
 * Spoken Language Systems Group
 * MIT Computer Science and Artificial Intelligence Laboratory
 * Massachusetts Institute of Technology
 *
 * FILE: interface.h 
 *   				      				                          
 * Chia-ying (Jackie) Lee <chiaying@csail.mit.edu>				
 * Jan 2016 
*********************************************************************/

#ifndef INTERFACE_H
#define INTERFACE_H

#include <map>
#include <string>
#include <utility>
#include <vector>
#include <tr1/unordered_map>

#include "cluster.h"
#include "config.h"
#include "datum.h"
#include "sampler.h"
#include "sym.h"
#include "trie.h"

typedef symbol S;
typedef pair<symbol, std::vector<symbol> > SSs;
typedef map<S, float> S_FL;
typedef tr1::unordered_map<S, S_FL> S_S_FL;

typedef trie<S, S_FL> St_S_FL;
typedef St_S_FL::const_iterator Stit;

class Interface {

 public:
     Interface(std::string&);
     ~Interface();
     bool load_speech_data(const std::string&);
     bool load_config(const std::string&);
     bool load_gaussian(const std::string&);
     bool load_bounds(Datum*, std::ifstream&, std::ifstream&);
     bool load_segs(Datum*, std::ifstream&);
     bool load_speech_model(const std::string&);
     bool load_clusters(const std::string&, \
             const std::string&, const std::string&);
     bool load_hdphmm_clusters(const std::string&);
     bool load_dphmm_clusters(const std::string&, const std::string&);

     bool retrieve_symbols(std::vector<std::vector<symbol> >&);
     bool load_model_id(const std::string&);
     std::vector<symbol> retrieve_symbols(const int);
     std::vector<std::vector<symbol> > resegment(const int, \
             trie<S, map<S, float> >&, trie<S, map<S, float> >&, \
             map<S, float>&, const std::vector<symbol>&, \
             const float pi0, const float r0);
     void initialize_state_mixture_seq();
     void insert_empty_strings();
     void resample_cluster_parameters();
     void show_clusters();
     void show_counters();
     void save_model(unsigned int);
     void save_data(unsigned int);
     bool initialize_counter();
 private:
     std::string get_tag(std::string);

     std::vector<Datum*> _data;
     std::map<symbol, Cluster*> _clusters;
     std::map<symbol, Cluster*> _counters;
     Config _config;
     Sampler _sampler;
     std::string _rootdir;
};

#endif
