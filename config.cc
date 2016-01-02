#include <fstream>
#include <iostream>
#include "config.h"

Config::Config() {
}

bool Config::load(const string& fn_config) {
    ifstream fconfig(fn_config.c_str(), ifstream::in);
    if (!fconfig.good()) {return false;}
    fconfig >> _dim;
    if (!fconfig.good()) {return false;}
    fconfig >> _max_duration; 
    if (!fconfig.good()) {return false;}
    fconfig >> _max_num_units;
    if (!fconfig.good()) {return false;}
    fconfig >> _max_num_epsilons;
    if (!fconfig.good()) {return false;}
    fconfig >> _trans_alpha;
    if (!fconfig.good()) {return false;}
    fconfig >> _gaussian_a0;
    if (!fconfig.good()) {return false;}
    fconfig >> _gaussian_k0;
    if (!fconfig.good()) {return false;}
    fconfig >> _mix_alpha; 
    if (!fconfig.good()) {return false;}
    fconfig >> _empty_string_prior; 
    fconfig.close();
    return true;
}

bool Config::load_gaussian(const string& fn_gaussian) {
    ifstream fgaussian(fn_gaussian.c_str(), ios::binary);
    if (!fgaussian.good()) {
        cout << "Cannot load Gaussian Prior" << endl;
        return false;
    }
    cout << "Loading Gaussian" << endl;
    float weight;
    fgaussian.read(reinterpret_cast<char*> (&weight), sizeof(float));
    float mean[_dim];
    float pre[_dim];
    fgaussian.read(reinterpret_cast<char*> (mean), sizeof(float) * _dim);
    fgaussian.read(reinterpret_cast<char*> (pre), sizeof(float) * _dim);
    _gaussian_u0.assign(mean, mean + _dim);
    _gaussian_b0.assign(pre, pre + _dim);
    for (int i = 0; i < _dim; ++i) {
        _gaussian_b0[i] = _gaussian_a0 / _gaussian_b0[i];
    }
    fgaussian.close(); 
    return true;
}

Config::~Config() {
}
