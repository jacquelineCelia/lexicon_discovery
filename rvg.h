#ifndef RVG_H
#define RVG_H

#include <ctime>

#include <boost/random/mersenne_twister.hpp>  
#include <boost/random/variate_generator.hpp>             
#include <boost/random/linear_congruential.hpp>
#include <boost/random/gamma_distribution.hpp> 
#include <boost/random/normal_distribution.hpp> 
#include <boost/random/uniform_real.hpp> 

typedef boost::mt19937 base_generator_type;

class RandomVarGen{
 public:
     RandomVarGen();
     float GetUniformSample();
     float GetGammaSample(float a, float b);
     float GetNormalSample(float u, float v);
     ~RandomVarGen();
 private:
     base_generator_type generator;
};

#endif
