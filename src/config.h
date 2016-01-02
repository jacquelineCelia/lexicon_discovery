#ifndef CONFIG_H
#define CONFIG_H

#include <string>
#include <vector>

using namespace std;

class Config {
 public:
    bool load(const string&);
    bool load_gaussian(const string& fn_gaussian);
    int dim() {return _dim;}
    int max_duration() {return _max_duration;}
    int max_num_units() {return _max_num_units;}
    int max_num_epsilons() {return _max_num_epsilons;}
    float transition_alpha() {return _trans_alpha;}
    float gaussian_a0() {return _gaussian_a0;}
    float gaussian_k0() {return _gaussian_k0;}
    float mix_alpha() {return _mix_alpha;}
    float empty_string_prior() {return _empty_string_prior;}
    vector<float> gaussian_b0() {return _gaussian_b0;}
    vector<float> gaussian_u0() {return _gaussian_u0;}
    Config();
    ~Config();
 private:
    int _dim;
    int _max_duration;
    int _max_num_units;
    int _max_num_epsilons;
    float _trans_alpha;
    float _gaussian_a0;
    float _gaussian_k0;
    float _mix_alpha;
    float _empty_string_prior;
    vector<float> _gaussian_b0;
    vector<float> _gaussian_u0;
};

#endif
