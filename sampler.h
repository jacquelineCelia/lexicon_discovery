#ifndef SAMPLER_H
#define SAMPLER_H

#include <map>
#include <vector>
#include <mkl.h>
#include <tr1/unordered_map>

#include "config.h"
#include "cluster.h"
#include "config.h"
#include "datum.h"
#include "prob_list.h"
#include "rvg.h"
#include "sym.h"
#include "trie.h"
#include "utility.h"

#define RHS_DEBUG false 
#define SEG_DEBUG false 
#define PLU_DEBUG false 
#define DEBUG false 
#define MSG_DEBUG false 
#define B_DEBUG false 

using namespace std;

typedef vector<symbol> Ss;
typedef symbol S;
typedef pair<symbol, std::vector<symbol> > SSs;

typedef map<S, float> S_FL;
typedef tr1::unordered_map<S, S_FL> S_S_FL;
typedef trie<S, S_FL> St_S_FL;
typedef St_S_FL::const_iterator Stit;

class Sampler {
 public:
  Sampler();
  void set_config(Config* config) {_config = config;}
  vector<vector<symbol> > Resegment(Datum*, \
          trie<S, S_FL>&, trie<S, S_FL>&, S_FL&, \
          const vector<symbol>&, map<symbol, Cluster*>&, \
          map<symbol, Cluster*>&, const float pi0, const float r0);

  void ComputeSegProbGivenCluster(Datum*, const map<symbol, Cluster*>&, \
          trie<S, vector<vector<float> > >&) ;
  
  void MessageBackward(const Ss&, Datum*, \
          map<symbol, vector<vector<float> > >&, \
          trie<S, S_FL>&, \
          vector<vector<ProbList<int> > >&, \
          vector<vector<ProbList<int> > >&);

  void SampleForward(const Ss&, Datum* datum, \
          const Ss&, \
          map<symbol, vector<vector<float> > >&, \
          trie<S, vector<vector<float> > >&, \
          trie<S, S_FL>&, \
          vector<vector<ProbList<int> > >&,
          vector<vector<ProbList<int> > >&,
          vector<Ss>&, vector<Segment*>&);

  Ss SampleRHS(symbol, const float, \
          const Ss&, \
          trie<S, vector<vector<float> > >& prob_clusters, \
          trie<S, S_FL>&, \
          vector<Bound*>&, const int, const int);

  void ComputeSegProbGivenMultiCs(Datum*, const Ss&, \
          trie<S, S_FL>&, trie<S, vector<vector<float> > >&); 

  int SampleIndexFromLogDistribution(vector<float>);
  int SampleIndexFromDistribution(vector<float>);

  void SampleStateMixtureSeq(Segment*, Cluster*);

  void RemoveClusterAssignment(Datum*, map<symbol, Cluster*>&);
  void AddClusterAssignment(Datum*, map<symbol, Cluster*>&); 

  void SampleClusterParams(Cluster*, Cluster* = NULL);

  vector<float> SampleDirFromGamma(int, float*, float = -70000000);

  vector<float> SampleGaussianPre(vector<float>, vector<float>, float);

  float UpdateGammaRate(float, float, float, float, float);

  vector<float> SampleGaussianMean(vector<float>, vector<float>, float);
  
  void InsertEmptyStrings(Datum*);

  void FastComputeSegProbGivenPluTop(Datum*, \
          const vector<symbol>&, \
          const vector<symbol>&, \
          trie<S, S_FL>&, \
          trie<S, vector<vector<float> > >&, \
          map<symbol, vector<vector<float> > >&);

  float PLUT2BProb(const Ss& plu_top, const vector<Ss>& plu_btm, \
          trie<S, S_FL>& ,S_FL&);
  float IncreasePLUT2B(const Ss& plu_top, const vector<Ss>& plu_btm, \
          trie<S, S_FL>&, S_FL&);

  ~Sampler();
 private:
  vector<symbol> GetPluCandidates(map<symbol, Cluster* >&);
  vector<symbol> GetPluCandidates(map<symbol, vector<vector<float> > >&);

  RandomVarGen rvg;
  Config* _config;
  VSLStreamStatePtr stream; 
};

#endif
