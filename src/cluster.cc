/* -*- C++ -*-
 *
 * Copyright (c) 2016
 * Spoken Language Systems Group
 * MIT Computer Science and Artificial Intelligence Laboratory
 * Massachusetts Institute of Technology
 *
 * FILE: cluster.cc
 *   				      				                            
 * Chia-ying (Jackie) Lee <chiaying@csail.mit.edu>				    
 * Jan 2016 
*********************************************************************/
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <omp.h>

#include "cluster.h"
#include "toolkit.h"

#define DEBUG false 
#define MIN_PROB_VALUE -70000000

Cluster::Cluster(Config* config) {
    _config = config;
}

Cluster::Cluster(Config* config, int state_num, int mixture_num) {
    _config = config;
    _state_num = state_num;
    _trans_prob.resize(state_num);
    for (int i = 0 ; i < state_num; ++i) {
        _trans_prob[i].resize(state_num + 1, 0);
    }
    for (int i = 0; i < state_num; ++i) {
        GMM emission(_config, mixture_num);
        _emissions.push_back(emission);
    }
}

void Cluster::set_state_weights(int i, vector<float>& weight) {
    _emissions[i].set_weight(weight);
}

void Cluster::set_state_mixture_det(int i, int j, float det) {
    _emissions[i].set_mixture_det(j, det);
}

void Cluster::set_state_mixture_mean(int i, int j, vector<float>& mean) {
    _emissions[i].set_mixture_mean(j, mean);
}

void Cluster::set_state_mixture_pre(int i, int j, vector<float>& pre) {
    _emissions[i].set_mixture_pre(j, pre);
}

void Cluster::load(ifstream& fin) {
    fin.read(reinterpret_cast<char*> (&_id), sizeof(int));
    set_id(_id);
    int fixed;
    fin.read(reinterpret_cast<char*> (&fixed), sizeof(int));
    fin.read(reinterpret_cast<char*> (&_state_num), sizeof(int));
    // Initialize space for _transition_probs and _emissions
    for (int i = 0; i < _state_num; ++i) {
        GMM emission(_config);
        _emissions.push_back(emission);
        vector<float> trans_prob(_state_num + 1, 0);
        _trans_prob.push_back(trans_prob);
    }
    for (int i = 0; i < _state_num; ++i) {
        fin.read(reinterpret_cast<char*> (&_trans_prob[i][0]), \
                sizeof(float) * (_state_num + 1));
    }
    for (int i = 0; i < _state_num; ++i) {
        _emissions[i].load(fin);
    }
}

void Cluster::print() {
    cout << "Cluster: " << _id_symbol << endl;
    for (int i = 0; i < _state_num; ++i) {
        for (int j = i; j <= _state_num; ++j) {
            cout << "trans[" << i << "][" << j << "]: " << _trans_prob[i][j] << " ";
        }
        cout << endl;
    }
    for (size_t i = 0; i < _emissions.size(); ++i) {
        _emissions[i].print();
    }
}

void Cluster::set_id(int id) {
    _id = id;
    stringstream s;
    s << _id;
    symbol sym(s.str());
    _id_symbol = sym;
}

vector<vector<float> > Cluster::ConstructSegProbTable(const vector<Bound*>& bounds) {
    int b = bounds.size();
    int total_frame_num = 0;
    vector<int> accumulated_frame_nums(b, 0);
    int start_frame = bounds[0] -> start_frame();
    int end_frame = bounds[b - 1] -> end_frame();
    vector<float*> frames(end_frame - start_frame + 1, NULL); 
    for (int i = 0 ; i < b; ++i) {
        vector<float*> bound_data = bounds[i] -> data();
        memcpy(&frames[total_frame_num], &bound_data[0], \
                (bounds[i] -> get_frame_num()) * sizeof(float*));
        total_frame_num += bounds[i] -> get_frame_num();
        accumulated_frame_nums[i] = total_frame_num;
    }
    if (end_frame - start_frame + 1 != total_frame_num) {
        cerr << "Miss matched total frame number!" << endl;
        exit(1);
    }
    float** frame_prob_for_each_state;
    frame_prob_for_each_state = new float* [_state_num];

    for (int i = 0 ; i < _state_num; ++i) {
        frame_prob_for_each_state[i] = new float[total_frame_num];
    }
    // parallel
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
        for (int i = 0 ; i < _state_num; ++i) {
            _emissions[i].ComputeLikehood(frames, \
                    frame_prob_for_each_state[i]);
        }

    if (DEBUG) {
        for (int i = 0; i < _state_num; ++i) {
            for (int j = 0; j < total_frame_num; ++j) {
                cout << frame_prob_for_each_state[i][j] << " ";
            }
            cout << endl;
        }
    }
    vector<float> vec_min_value(b, MIN_PROB_VALUE);
    vector<vector<float> > prob_table(b, vec_min_value);
    int max_duration = _config -> max_duration();

    // parallel
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
        for (int i = 0 ; i < b; ++i) {
            int j = i;
            int start_frame = i == 0 ? 0 : accumulated_frame_nums[i - 1];
            int duration = accumulated_frame_nums[i] - start_frame; 
            int ptr = start_frame;
            vector<float> cur_prob(_state_num, MIN_PROB_VALUE);
            while (ptr < start_frame + duration && j < b \
                    && (duration <= max_duration || (int) i == j)) {
                if (ptr == start_frame) {
                    cur_prob[0] = frame_prob_for_each_state[0][ptr];
                }
                else {
                    vector<float> next_prob(_state_num, 0);
                    for (int k = 0; k < _state_num; ++k) {
                        float summand = MIN_PROB_VALUE;
                        for (int l = 0; l <= k; ++l) {
                            summand = ToolKit::SumLogs(cur_prob[l] + _trans_prob[l][k], summand);
                        }
                        next_prob[k] = summand + frame_prob_for_each_state[k][ptr];
                    }
                    cur_prob = next_prob;
                }
                if (ptr == accumulated_frame_nums[j] - 1) {
                    float next_prob = MIN_PROB_VALUE;
                    for (int k = 0; k < _state_num; ++k) {
                        next_prob = ToolKit::SumLogs(cur_prob[k] + _trans_prob[k][_state_num], next_prob);
                    }
                    prob_table[i][j] = next_prob; 
                    if (++j < b) {
                        duration = accumulated_frame_nums[j] - start_frame;
                    }
                }
                ++ptr;
            }
        }
    // }
    for (int i = 0; i < b; ++i) {
        for (int j = 0; j < b; ++j) {
            if (isnan(prob_table[i][j])) {
               cerr << "Found NaN in  ConstructSegProbTable. Cluster id: " << _id << endl;
               exit(-1);
            }
        }
    }
    for (int i = 0 ; i < _state_num; ++i) {
        delete[] frame_prob_for_each_state[i];
    }
    delete [] frame_prob_for_each_state;
    return prob_table;
}

vector<vector<ProbList<int> > > Cluster::MessageBackwardForASegment(Segment* segment) {
    int frame_num = segment -> frame_num();
    vector<vector<ProbList<int> > > B;
    B.resize(_state_num + 1);
    for (int i = 0 ; i <= _state_num; ++i) {
        B[i].resize(frame_num + 1);
    }
    // Initialization [need to check what the initial value should be!]
    for (int i = 1; i <= _state_num; ++i) {
        B[i][frame_num].push_back(_trans_prob[i - 1][_state_num], -1);
    }
    // Message Backward
    for (int j = frame_num - 1; j > 0; --j) {
       float* data = segment -> frame(j);
       vector<float> emit_probs(_state_num);
       for (int k = 0; k < _state_num; ++k) {
           float emit_prob = _emissions[k].ComputeLikehood(data);
           emit_probs[k] = emit_prob;
       }
       for (int i = 1; i <= _state_num; ++i) {
           for (int k = i; k <= _state_num; ++k) {
               B[i][j].push_back(_trans_prob[i - 1][k - 1] + \
                   emit_probs[k - 1] + B[k][j + 1].value(), k);
           }
       } 
    }
    float emit_prob = _emissions[0].ComputeLikehood(segment -> frame(0));
    B[0][0].push_back(emit_prob + B[1][1].value(), 1); 
    return B;
}

void Cluster::Minus(Segment* segment) {
    vector<int> state_seq = segment -> state_seq();
    vector<int> mix_seq = segment -> mix_seq();
    vector<float*> data = segment -> data();
    if (state_seq.size() != data.size()) {
        cout << "In ClusterCounter::Minus, state_seq and data have different sizes." << endl;
        exit(2);
    }
    else if (mix_seq.size() != data.size()) {
        cout << "In ClusterCounter::Minus, mix_seq and data have different sizes." << endl;
        exit(2);
    }
    else {
        for (size_t i = 0; i < state_seq.size(); ++i) {
            int cur_state = state_seq[i];
            int next_state = i == state_seq.size() - 1 ? \
                             _state_num : state_seq[i + 1];
            --_trans_prob[cur_state][next_state];
            _emissions[cur_state].Minus(data[i], mix_seq[i]);
        }
    }
}

void Cluster::Plus(Segment* segment) {
    const vector<int> state_seq = segment -> state_seq(); 
    const vector<int> mix_seq = segment -> mix_seq();
    const vector<float*> data = segment -> data();
    if (state_seq.size() != data.size()) {
        cout << "In ClusterCounter::Plus, state_seq and data have different sizes." << endl;
        exit(2);
    }
    else if (mix_seq.size() != data.size()) {
        cout << "In ClusterCounter::Plus, mix_seq and data have different sizes." << endl;
        exit(2);
    }
    else {
        for (size_t i= 0 ; i < state_seq.size(); ++i) {
            int cur_state = state_seq[i];
            int next_state = i == state_seq.size() - 1 ? \
                             _state_num : state_seq[i + 1];
            ++_trans_prob[cur_state][next_state];
            _emissions[cur_state].Plus(data[i], mix_seq[i]);
        }
    }
}

void Cluster::Save(ofstream& fout) {
    fout.write(reinterpret_cast<char*> (&_id), sizeof(int));
    fout.write(reinterpret_cast<char*> (&_state_num), sizeof(int));
    for (int i = 0; i < _state_num; ++i) {
        fout.write(reinterpret_cast<char*> (&_trans_prob[i][0]), sizeof(float) * (_state_num + 1));
    }
    for (int i = 0 ; i < _state_num; ++i) {
        _emissions[i].Save(fout);
    }
}

Cluster::~Cluster() {
}

