#include "rvg.h"

using namespace std;
using namespace boost;
using namespace boost::random;

RandomVarGen::RandomVarGen() {
    // size_t SEED = time(0);
    size_t SEED = 0;
    generator.seed(SEED);
}

/*
float RandomVarGen::GetGammaSample(float a, float b) {
    gamma_distribution<> dist(a, b);
    return dist(generator);
}
*/

float RandomVarGen::GetUniformSample() {
    uniform_real<> dist(0, 1);
    float sample = dist(generator);
    while (sample == 0 || sample == 1) {
        sample = dist(generator);
    }
    return sample;
}

float RandomVarGen::GetNormalSample(float u, float v) {
    normal_distribution<> dist(u, v);
    return dist(generator);
}

RandomVarGen::~RandomVarGen() {
}
