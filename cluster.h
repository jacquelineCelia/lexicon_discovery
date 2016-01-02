#ifndef CLUSTER_H
#define CLUSTER_H

#include "bound.h"
#include "config.h"
#include "gmm.h"
#include "prob_list.h"
#include "sym.h"
#include "segment.h"

#include <fstream>
#include <iostream>
#include <vector>

using namespace std;

class Cluster {
 public:
     Cluster(Config*);
     Cluster(Config*, int, int);
     void load(ifstream&); 
     void set_trans(vector<vector<float> > t) {_trans_prob = t;}
     void set_state_mixture_det(int, int, float);
     void set_state_mixture_mean(int, int, vector<float>&);
     void set_state_mixture_pre(int, int, vector<float>&);
     void set_state_weights(int, vector<float>&);
     void set_id(int id); 
     void set_id_symbol(symbol id_symbol) {_id_symbol = id_symbol;}
     void Plus(Segment*);
     void Minus(Segment*);
     void print();
     void Save(ofstream& fout);

     symbol id_symbol() {return _id_symbol;}
     int id() {return _id;}

     int state_num() {return _state_num;}
     int mix_num() {return _emissions.size() ? _emissions[0].mix_num() : 0;}
     vector<vector<float> > transition_probs() {return _trans_prob;}
     vector<vector<float> > ConstructSegProbTable(const vector<Bound*>&);
     vector<vector<ProbList<int> > > MessageBackwardForASegment(Segment*);
     GMM& emission(int index) {return _emissions[index];}
     ~Cluster();
 private:
     Config* _config;
     vector<vector<float> > _trans_prob;
     vector<GMM> _emissions;
     int _state_num;
     int _id;
     symbol _id_symbol;
};

#endif
