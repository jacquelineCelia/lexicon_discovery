/* -*- C++ -*-
 *
 * Copyright (c) 2014
 * Spoken Language Systems Group
 * MIT Computer Science and Artificial Intelligence Laboratory
 * Massachusetts Institute of Technology
 *
 * FILE: interface.cc
 *   				      				                          
 * Chia-ying (Jackie) Lee <chiaying@csail.mit.edu>				
 * Jan 2016 
*********************************************************************/
#include <sstream>
#include <iostream>
#include <fstream>
#include <sys/stat.h>

#include "cluster.h"
#include "interface.h"

using namespace std;

Interface::Interface(std::string& rootdir) {
    _rootdir = rootdir;
}

string Interface::get_tag(string s) {
   size_t found_last_slash, found_last_period;
   found_last_slash = s.find_last_of("/");
   found_last_period = s.find_last_of(".");
   return s.substr(found_last_slash + 1, \
     found_last_period - 1 - found_last_slash);
}

void Interface::insert_empty_strings() {
    for (size_t i = 0; i < _data.size(); ++i) {
        _sampler.InsertEmptyStrings(_data[i]);
    }
}

bool Interface::load_speech_data(const string& filename) {
    ifstream flist(filename.c_str(), ifstream::in);
    if (!flist.is_open()) {
        return false;
    }
    while (flist.good()) {
        string fn_index, fn_data, fn_segs;
        flist >> fn_index >> fn_data >> fn_segs;
        if (fn_index != "" && fn_data != "" && fn_segs != "") {
            ifstream findex(fn_index.c_str(), ifstream::in);
            ifstream fdata(fn_data.c_str(), ifstream::binary);
            ifstream fsegs(fn_segs.c_str(), ifstream::in);
            if (!findex.is_open()) {
                cerr << "Cannot open " << findex << endl;
                return false;
            }
            if (!fdata.is_open()) {
                cerr << "Cannot open " << fdata << endl;
                return false;
            }
            if (!fsegs.is_open()) {
                cerr << "Cannot open " << fsegs << endl;
                return false;
            }
            cout << "Loading " << fn_data << "..." << endl;
            Datum* datum = new Datum(&_config);
            string tag = get_tag(fn_data);
            datum -> set_tag(tag);
            if (!load_bounds(datum, findex, fdata)) {
                cerr << "Cannot load speech data" << fn_data << endl;
                return false;
            }
            if (!load_segs(datum, fsegs)) {
                cerr << "Cannot load segs" << fn_segs << endl;
                return false;
            }
            _data.push_back(datum);
            findex.close();
            fdata.close();
            fsegs.close();
        }
   }
   flist.close();
   return true;
}

bool Interface::load_bounds(Datum* datum, ifstream& findex, ifstream& fdata) {
    // Goal : get the data for bounds
    vector<Bound*> bounds;
    int total_frame;
    findex >> total_frame;
    int start_frame = 0;
    int end_frame = 0;
    int label;
    int dim = _config.dim();
    while (end_frame != total_frame - 1) {
        findex >> start_frame >> end_frame >> label;
        int frame_num = end_frame - start_frame + 1;
        Bound* bound = new Bound(dim, start_frame, end_frame);
        // float** data is deleted in bound.cc
        float** data = new float* [frame_num];
        for (int i = 0; i < frame_num; ++i) {
            data[i] = new float [dim];
            fdata.read(reinterpret_cast<char*> (data[i]), \
                    sizeof(float) * dim);
         }
         bound -> set_data(data, frame_num);
         bounds.push_back(bound);
    }
    datum -> set_bounds(bounds);
    datum -> set_frame_num(total_frame);
    return true;
}

bool Interface::load_segs(Datum* datum, ifstream& fsegs) {
    vector<Segment*> segments;
    vector<Bound*> bounds = datum -> bounds();
    int seg_s = 0, seg_e = 0, seg_label;
    size_t ptr = 0;
    while (fsegs >> seg_s >> seg_e >> seg_label) {
        Segment* seg = new Segment(seg_label, seg_e - seg_s + 1);
        while (ptr < bounds.size() && bounds[ptr] -> start_frame() >= seg_s && \
          bounds[ptr] -> end_frame() <= seg_e) {
            seg -> push_back(bounds[ptr]);
            ++ptr;
        }
        segments.push_back(seg);
    }
    datum -> set_segments(segments);
    return true;
}

bool Interface::load_config(const string& fn_config) {
    bool result = _config.load(fn_config);
    _sampler.set_config(&_config); 
    return result; 
}

bool Interface::load_gaussian(const string& fn_gaussian) {
    return _config.load_gaussian(fn_gaussian);
}

bool Interface::load_clusters(const string& fn_model, \
                    const string& type, const string& model_id) {
    cout << type << endl;
    if (type == "hdphmm") {
        return load_hdphmm_clusters(fn_model);
    }
    else if (type == "dphmm") {
        return load_dphmm_clusters(fn_model, model_id);
    }
    else {
        cerr << "Undefined model type: " << type 
            << ". Only \"hdphmm\" or \"dphmm\" is allowed" << endl; 
        return false;
    }
}

bool Interface::load_hdphmm_clusters(const string& fn_snapshot) {
    _clusters.clear();
    // Read the new model from snapshot
    ifstream fsnapshot(fn_snapshot.c_str(), ios::binary);
    int num_clusters;
    fsnapshot.read(reinterpret_cast<char*> (&num_clusters), sizeof(int));
    cout << "Number of clusters: " << num_clusters << endl;
    float pi[num_clusters];
    fsnapshot.read(reinterpret_cast<char*> (pi), sizeof(float) * num_clusters);
    float beta[num_clusters];
    fsnapshot.read(reinterpret_cast<char*> (beta), sizeof(float) * num_clusters);
    float A[num_clusters][num_clusters];
    for (int i = 0; i < num_clusters; ++i) {
        fsnapshot.read(reinterpret_cast<char*>(A[i]), sizeof(float) * num_clusters);
    }
    for (int i = 0; i < num_clusters; ++i) {
        Cluster* c = new Cluster(&_config);
        c -> load(fsnapshot);
        _clusters[c -> id_symbol()] = c;
    }
    fsnapshot.close();
    if (!initialize_counter()) {
        cerr << "Can't initialize_counter" << endl;
        return false;
    }
    return true;
}

void Interface::initialize_state_mixture_seq() {
    vector<Datum*>::iterator d_iter = _data.begin();
    for (; d_iter != _data.end(); ++d_iter) {
        vector<Segment*> segs = (*d_iter) -> segments();
        for (size_t i = 0; i < segs.size(); ++i) {
            if ((segs[i] -> id_symbol()).string_reference() != "\'\'") {
                _sampler.SampleStateMixtureSeq(segs[i], \
                    _clusters[segs[i] -> id_symbol()]);
            }
        }
        _sampler.AddClusterAssignment(*d_iter, _counters);
    }
}

bool Interface::load_dphmm_clusters(const string& fn_model, \
                                    const string& model_id) {
    ifstream fid(model_id.c_str());
    if (!fid.is_open()) {
        return false;
    }
    int id;
    ifstream fin(fn_model.c_str(), ios::binary);
    int data_num, cluster_num;
    if (!fin.good()) {
        cout << fn_model << " cannot be opened." << endl;
        return false;
    }
    fin.read(reinterpret_cast<char*> (&data_num), sizeof(int));
    fin.read(reinterpret_cast<char*> (&cluster_num), sizeof(int));
    cout << "number of clusters " << cluster_num << endl;
    for (int i = 0; i < cluster_num; ++i) {
        fid >> id;
        int member_num, state_num, mixture_num, vector_dim;
        fin.read(reinterpret_cast<char*> (&member_num), sizeof(int));
        fin.read(reinterpret_cast<char*> (&state_num), sizeof(int));
        fin.read(reinterpret_cast<char*> (&mixture_num), sizeof(int));
        fin.read(reinterpret_cast<char*> (&vector_dim), sizeof(int));
        Cluster* new_cluster = new Cluster(&_config, state_num, mixture_num);
        vector<vector<float> > trans;
        for (int j = 0; j < state_num; ++j) {
            vector<float> prob(state_num + 1, 0);
            fin.read(reinterpret_cast<char*> (&prob[0]), \
                    sizeof(float) * (state_num + 1));
            trans.push_back(prob);
        }
        new_cluster -> set_trans(trans);
        for (int j = 0; j < state_num; ++j) {
            vector<float> weights(mixture_num, 0);
            for (int k = 0 ; k < mixture_num; ++k) {
                float w, det;
                vector<float> mean(vector_dim, 0);
                vector<float> pre(vector_dim, 0);
                fin.read(reinterpret_cast<char*> (&w), \
                    sizeof(float));
                fin.read(reinterpret_cast<char*> (&det), \
                    sizeof(float));
                fin.read(reinterpret_cast<char*> (&mean[0]), \
                    sizeof(float) * vector_dim);
                fin.read(reinterpret_cast<char*> (&pre[0]), \
                    sizeof(float) * vector_dim);
                new_cluster -> set_state_mixture_det(j, k, det);
                new_cluster -> set_state_mixture_mean(j, k, mean);
                new_cluster -> set_state_mixture_pre(j, k, pre);
                weights.push_back(w);
            }
            new_cluster -> set_state_weights(j, weights);
        }
        new_cluster -> set_id(id);
        _clusters[new_cluster -> id_symbol()] = new_cluster;
    }
    fin.close();
    if (!initialize_counter()) {
        cerr << "Can't initialize counters" << endl;
        return false;
    }
    return true;
}

bool Interface::initialize_counter() {
    map<symbol, Cluster*>::iterator iter = _clusters.begin();
    for (; iter != _clusters.end(); ++iter) {
        int state_num = (iter -> second) -> state_num();
        int mix_num = (iter -> second) -> mix_num();
        if (mix_num <= 0) {
            cerr << "mix num: " << mix_num << endl;
            return false;
        }
        if (state_num <= 0) {
            cerr << "state num: " << state_num << endl;
            return false;
        }
        Cluster* new_counter = new Cluster(&_config, state_num, mix_num);
        new_counter -> set_id_symbol(iter -> first);
        _counters[iter -> first] = new_counter;
    } 
    return true;
}

/*
bool Interface::load_model_id(const string& model_id) {
    ifstream fid(model_id.c_str());
    if (!fid.is_open()) {
        return false;
    }
    size_t ptr = 0;
    int id;
    while(fid >> id) {
        _clusters[ptr++] -> set_id(id);
    }
    if (ptr != _clusters.size()) {
        cerr << ptr << " " << _clusters.size() << endl;
        cerr << "id number doesn't match model number" << endl;
        return false;
    }
    fid.close();
    return true;
}
*/

bool Interface::retrieve_symbols(vector<vector<symbol> >& trains) {
    trains.clear();
    for (size_t i = 0; i != _data.size(); ++i) {
        vector<symbol> symbols = retrieve_symbols(i);
        if (!symbols.size()) {
            std::cerr << "got empty terminals in retrieve_symbols." << endl;
            return false;
        }
        trains.push_back(symbols);
    }
    return true;
}

vector<symbol> Interface::retrieve_symbols(const int index) {
    return _data[index] -> retrieve_symbols();
}

vector<vector<symbol> > Interface::resegment(const int index, \
  trie<S, S_FL>& rules, trie<S, S_FL>& rule_counts, \
  S_FL& parent_counts,  const vector<symbol>& plu_top, \
  const float pi0, const float r0) {
  return _sampler.Resegment(_data[index], \
          rules, rule_counts, parent_counts, \
          plu_top, _clusters, _counters, pi0, r0);
}

void Interface::show_clusters() {
    map<symbol, Cluster*>::iterator c_iter = _clusters.begin();
    for (; c_iter != _clusters.end(); ++c_iter) {
        c_iter -> second -> print();
    }
}

void Interface::show_counters() {
    map<symbol, Cluster*>::iterator c_iter = _counters.begin();
    for (; c_iter != _counters.end(); ++c_iter) {
        c_iter -> second -> print();
    }
}

void Interface::resample_cluster_parameters() {
    map<symbol, Cluster*>::iterator iter = _clusters.begin();
    for (; iter != _clusters.end(); ++iter) {
        _sampler.SampleClusterParams(iter -> second, _counters[iter -> first]);
    }
}

void Interface::save_data(unsigned int iteration) {
    stringstream iter;
    iter << iteration;
    string output_dir = _rootdir + "/" + iter.str() + "/";
    mkdir(output_dir.c_str(), S_IRWXU | S_IRWXO | S_IRWXG);

    #pragma omp parallel
    {
        #pragma omp for schedule (dynamic, 1)
        for (int d = 0; d < (int) _data.size(); ++d) {
            _data[d] -> Save(output_dir);
        }
    }
}

void Interface::save_model(unsigned int iteration) {
    stringstream iter;
    iter << iteration;

    string path = _rootdir + "/" + iter.str() + "/" + "snapshot";
    ofstream fout(path.c_str(), ios::binary);

    map<symbol, Cluster*>::iterator c_iter = _clusters.begin();
    for (; c_iter != _clusters.end(); ++c_iter) {
        c_iter -> second -> Save(fout);
    }
    fout.close();
}

Interface::~Interface() {
    vector<Datum*>::iterator iter = _data.begin();
    for (; iter != _data.end(); ++iter) {
        delete *iter;
    }

    map<symbol, Cluster*>::iterator c_iter = _clusters.begin();
    for (; c_iter != _clusters.end(); ++c_iter) {
        delete c_iter -> second;
    }

    map<symbol, Cluster*>::iterator counter_iter = _counters.begin();
    for (; counter_iter != _counters.end(); ++counter_iter) {
        delete counter_iter -> second;
    }
}

