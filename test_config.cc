#include <iostream>
#include <string>

#include "config.h"

using namespace std;

int main() {
    string filename("configuration");
    Config config;

    config.load(filename);

    cout << "dimension: " << config.dim() << endl;
    cout << "max duration: " << config.max_duration() << endl;
    cout << "max num units: " << config.max_num_units() << endl;
    cout << "max num epsilons: " << config.max_num_epsilons() << endl;
    cout << "trans alpha: " << config.transition_alpha() << endl; 
    cout << "gaussian a0: " << config.gaussian_a0() << endl;
    cout << "gaussian k0: " << config.gaussian_k0() << endl;
    cout << "mix alpha: " << config.mix_alpha() << endl;
    cout << "empty string prior: " << config.empty_string_prior() << endl;

    config.load_gaussian("gaussian");

    vector<float> gaussian_u0 = config.gaussian_u0();
    vector<float> gaussian_b0 = config.gaussian_b0();

    for (size_t i = 0; i < gaussian_u0.size(); ++i) {
        cout << gaussian_u0[i] << " " << gaussian_b0[i] << endl;
    }

    return 0;
}
