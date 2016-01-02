// py-cky.h
//
// (c) Mark Johnson, 27th January 2006, last modified 2nd May, 2013

#ifndef PY_CKY_H
#define PY_CKY_H

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
// #include <ext/hash_map>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <utility>
#include <vector>
#include <tr1/unordered_map>
#include <omp.h>

#include "earley.h"
#include "gammadist.h"
#include "mt19937ar.h"
#include "slice-sampler.h"
#include "sym.h"
#include "xtree.h"
#include "trie.h"
#include "utility.h"

extern int debug;

//! Suppose there are n samples occupying m tables.
//! Then the probability that the n+1 sample occupies
//! table 1 <= k <= m is:
//!
//!  P(x_{n+1} = k) = (n_k - a)/(n + b)
//!
//! and the probability that the n+1 sample occupies
//! the new table m+1 is:
//!
//!  P(x_{n+1} = m+1) = (m*a + b)/(n + b)
//!
//! The probability of a configuration in which a 
//! restaurant contains n customers at m tables,
//! with n_k customers at table k is:
//!
//!
//!  a^{-m} G(m+b/a)  G(b)                 G(n_k-a)
//!         -------- ------  \prod_{k=1}^m --------
//!          G(b/a)  G(n+b)                 G(1-a)
//!
//! where G is the Gamma function.

inline float power(float x, float y) { return y == 1 ? x : powf(x, y); }
inline double power(double x, double y) { return y == 1 ? x : pow(x, y); }

#ifndef QUADPREC
  typedef double F;
#else
  #include "quadmath.h"
  typedef __float128 F;
  inline __float128 power(__float128 x, __float128 y) { return y == 1 ? x : pow(double(x), double(y)); }
  inline __float128 log(__float128 x) {return log(double(x));}
#endif

#define MIN_PROB_VALUE -70000000

// inline long double power(long double x, long double y) { return powl(x, y); }

typedef symbol S;
typedef std::vector<S> Ss;
typedef std::vector<Ss> Sss;

typedef std::map<S,F> S_F;
// typedef tr1::unordered_map<S,F> S_F;

typedef std::pair<S,Ss> SSs;
typedef std::map<SSs,F> SSs_F;

//! readline_symbols() reads all of the symbols on the current
//! line into syms
//
inline 
std::istream& readline_symbols(std::istream& is, Ss& syms) {
  syms.clear();
  std::string line;
  if (std::getline(is, line)) {
    std::istringstream iss(line);
    std::string s;
    while (iss >> s)
      syms.push_back(s);
  }
  return is;
}  // readline_symbols()


//! A default_value_type{} object is used to read an object from a stream,
//! assigning a default value if the read fails.  Users should not need to
//! construct such objects, but should use the default_value() function instead.
//
template <typename object_type, typename default_type>
struct default_value_type {
  object_type& object;
  const default_type defaultvalue;
  default_value_type(object_type& object, const default_type defaultvalue)
    : object(object), defaultvalue(defaultvalue) { }
};

//! default_value() is used to read an object from a stream, assigning a
//! default value if the read fails.  It returns a default_value_type{}
//! object, which does the actual reading.
//
template <typename object_type, typename default_type>
default_value_type<object_type,default_type>
default_value(object_type& object, const default_type defaultvalue=default_type()) {
  return default_value_type<object_type,default_type>(object, defaultvalue);
}

//! This operator>>() reads default_value_type{} from an input stream.
//
template <typename object_type, typename default_type>
std::istream& operator>> (std::istream& is, 
			  default_value_type<object_type, default_type> dv) {
  if (is) {
    if (is >> dv.object)
      ;
    else {
      is.clear(is.rdstate() & ~std::ios::failbit);  // clear failbit
      dv.object = dv.defaultvalue;
    }
  }
  return is;
}

// inline F random1() { return rand()/(RAND_MAX+1.0); }
inline F random1() { return mt_genrand_res53(); }


//! A pycfg_type is a CKY parser for a py-cfg
//
struct pycfg_type {

  pycfg_type() 
    : estimate_theta_flag(false), predictive_parse_filter(false), 
      default_weight(1), default_pya(1e-1), default_pyb(1e3),
      pya_beta_a(0), pya_beta_b(0), pyb_gamma_s(0), pyb_gamma_c(0) { }

  typedef unsigned int U;
  typedef std::pair<U,U> UU;

  typedef std::map<S,U> S_U;

  typedef std::map<S,UU> S_UU;

  typedef tr1::unordered_map<S,S_F> S_S_F;

  typedef trie<S, S_F> St_S_F;
  typedef St_S_F::const_iterator Stit;

  typedef catcounttree_type tree;

  typedef std::set<tree*> sT;

  typedef trie<S,sT> St_sT;

  typedef std::vector<tree*> Ts;

  typedef std::map<S,Ts> S_Ts;

  //! If estimate_theta_flag is true, then we estimate the generator 
  //! rule weights using a Dirichlet prior
  //
  bool estimate_theta_flag;

  //! If predictive_parse_filter is true, then first do a deterministic
  //! Earley parse of each sentence and use this to filter the nondeterministic
  //! CKY parses
  //
  bool predictive_parse_filter;

  //! predictive_parse_filter_grammar is the grammar used by the Earley parser
  //
  earley::grammar predictive_parse_filter_grammar;

  //! start is the start symbol of the grammar
  //
  S start;

  //! rhs_parent_weight maps the right-hand sides of rules
  //! to rule parent and rule weight 
  //
  St_S_F rhs_parent_weight;
  St_S_F filtered_rhs_parent_weight;

  //! unarychild_parent_weight maps unary children to a vector
  //! of parent-weight pairs
  //
  S_S_F unarychild_parent_weight;
  S_S_F filtered_unarychild_parent_weight;

  //! parent_weight maps parents to the sum of their rule weights
  //
  S_F parent_weight;

  //! default_weight is the default weight for rules with no explicit
  //! weight.  Used when grammar is read in.
  //
  F default_weight;

  //! rule_priorweight is the prior weight of rule
  //
  SSs_F rule_priorweight;

  //! parent_priorweight is the prior weight the parent
  //
  S_F parent_priorweight;

  //! terms_pytrees maps terminal strings to their PY trees
  //
  St_sT terms_pytrees;
  St_sT filtered_terms_pytrees;

  //! parent_pyn maps parents to the number of times they have been expanded
  //
  S_U parent_pyn;

  //! parent_pym maps parents to the number of distinct PY tables for parent
  //
  S_U parent_pym;

  //! rule table
  std::map<S, std::vector<Ss> > rules_look_up;


  F default_pya;   //!< default value for pya
  F default_pyb;   //!< default value for pyb

  F pya_beta_a;    //!< alpha parameter of Beta prior on pya
  F pya_beta_b;    //!< beta parameter of Beta prior on pya

  F pyb_gamma_s;   //!< s parameter of Gamma prior on pyb
  F pyb_gamma_c;   //!< c parameter of Gamma prior on pyb

  S_F parent_pya;  //!< pya value for parent
  S_F parent_pyb;  //!< pyb value for parent

  //! get_pya() returns the value of pya for this parent
  //
  F get_pya(S parent) const { 
    S_F::const_iterator it = parent_pya.find(parent);
    return (it == parent_pya.end()) ? default_pya : it->second;
  }  // pycfg_type::get_pya()

  //! set_pya() sets the value of pya for this parent, returning
  //! the old value for pya
  //
  F set_pya(S parent, F pya) {
    F old_pya = default_pya;
    S_F::iterator it = parent_pya.find(parent);
    if (it != parent_pya.end())
      old_pya = it->second;
    if (pya != default_pya)
      parent_pya[parent] = pya;
    else // pya == default_pya
      if (it != parent_pya.end())
	parent_pya.erase(it);
    return old_pya;
  }  // pycfg_type::set_pya()

  //! get_pyb() returns the value of pyb for this parent
  //
  F get_pyb(S parent) const { 
    S_F::const_iterator it = parent_pyb.find(parent);
    return (it == parent_pyb.end()) ? default_pyb : it->second;
  }  // pycfg_type::get_pyb()

  //! sum_pym() returns the sum of the pym for all parents
  //
  U sum_pym() const {
    U sum = 0;
    cforeach (S_U, it, parent_pym)
      sum += it->second;
    return sum;
  }  // pycfg_type::sum_pym()

  //! terms_pytrees_size() returns the number of trees in terms_pytrees.
  //
  U terms_pytrees_size() const {
    U size = 0;
    terms_pytrees.for_each(terms_pytrees_size_helper(size));
    return size;
  }  // pycfg_type::terms_pytrees_size()

  struct terms_pytrees_size_helper {
    U& size;
    terms_pytrees_size_helper(U& size) : size(size) { }

    template <typename Words, typename TreePtrs>
    void operator() (const Words& words, const TreePtrs& tps) {
      size += tps.size();
       // Jackie: TreePtrs -- set of trees
    }  // pycfg_type::terms_pytrees_size_helper::operator()    

  };  // pycfg_type::terms_pytrees_size_helper{}

  //! rule_weight() returns the weight of rule parent --> rhs
  //
  template <typename rhs_type>
  F rule_weight(S parent, const rhs_type& rhs) const {
      assert(!rhs.empty());
      if (rhs.size() == 1) {
          S_S_F::const_iterator it = unarychild_parent_weight.find(rhs[0]);
          if (it == unarychild_parent_weight.end())
              return 0;
          else
              return dfind(it -> second, parent);
    }
    else {  // rhs.size() > 1
        Stit it = rhs_parent_weight.find(rhs);
        if (it == rhs_parent_weight.end())
            return 0;
        else
            return dfind(it->data, parent);
    }
  }  // pycfg_type::rule_weight()

  //! rule_prob() returns the probability of rule parent --> rhs
  //
  template <typename rhs_type>
  F rule_prob(S parent, const rhs_type& rhs) const {
      assert(!rhs.empty());
      F parentweight = afind(parent_weight, parent);
      F ruleweight = rule_weight(parent, rhs);
      assert(ruleweight > 0);
      assert(parentweight > 0);
      return ruleweight/parentweight;
  }  // pycfg_type::rule_prob()

  F PLUT2BProb(const Ss& plu_top, const Sss& plu_btm) {
      F prob = 1;
      assert(plu_top.size() == plu_btm.size());
      for (size_t i = 0; i < plu_top.size(); ++i) {
          F ruleprob = rule_prob(plu_top[i], plu_btm[i]);
          if (debug == -1) {
              std::cerr << "tree_prob " 
                  << plu_top[i] << " --> " << plu_btm[i] 
                  << ", prob = " << ruleprob << std::endl;
          }
          prob *= ruleprob;
      }
      return prob;
  }

  //! tree_prob() returns the probability of the tree under the current
  //! model
  //
  F tree_prob(const tree* tp) const {
      if (tp -> children.empty()) 
          return 1;
      F pya = get_pya(tp->cat);
      if (pya == 1) { // no cache
          F prob = 1;
          Ss children;
          cforeach(tree::ptrs_type, it, tp->children) {
              children.push_back((*it)->cat);
              prob *= tree_prob(*it);
          }
          prob *= rule_prob(tp->cat, children);
          return prob;
      }
      F pyb = get_pyb(tp->cat);
      U pym = dfind(parent_pym, tp->cat);
      U pyn = dfind(parent_pyn, tp->cat);
      if (tp->count > 0) { // existing node
          assert(tp->count <= pyn);
          assert(pym > 0);
          F prob = (tp->count - pya)/(pyn + pyb);
          assert(finite(prob)); assert(prob > 0); assert(prob <= 1);
          return prob;
      }
      // new node
      F prob = (pym * pya + pyb)/(pyn + pyb);
      assert(finite(prob)); assert(prob > 0); assert(prob <= 1);
      Ss children;
      cforeach(tree::ptrs_type, it, tp->children) {
          children.push_back((*it)->cat);
          prob *= tree_prob(*it);
      }
      prob *= rule_prob(tp->cat, children);
      if (prob < 0) {
          std::cerr << "## pycfg_type::tree_prob(" << *tp << ") = " 
	  	  << prob << std::endl;
      }
      assert(finite(prob)); assert(prob <= 1); assert(prob >= 0);
      // assert(prob > 0); 
      return prob;
  }  // pycfg_type::tree_prob()

  //! incrrule() increments the weight of the rule parent --> rhs,
  //! returning the probability of this rule under the old grammar.
  //
  template <typename rhs_type>
  F incrrule(S parent, const rhs_type& rhs, F weight = 1) {
      assert(!rhs.empty());
      assert(weight >= 0);
      F& parentweight = parent_weight[parent];
      F parentweight0 = parentweight;
      F rhsweight0;
      parentweight += weight;
      if (rhs.size() == 1) {
          F& rhsweight = unarychild_parent_weight[rhs[0]][parent];
          rhsweight0 = rhsweight;
          rhsweight += weight;
      }
      else {  // rhs.size() > 1
          F& rhsweight = rhs_parent_weight[rhs][parent];
          rhsweight0 = rhsweight;
          rhsweight += weight;
      }
      assert(parentweight0 >= 0);
      assert(rhsweight0 >= 0);
      return rhsweight0/parentweight0;
  }  // incrrule()

  //! decrrule() decrements the weight of rule parent --> rhs,
  //! returning the probability of this rule under the new grammar,
  //! and deletes the rule if it has weight 0.
  //
  template <typename rhs_type>
  F decrrule(S parent, const rhs_type& rhs, F weight = 1) {
      assert(weight >= 0);
      assert(!rhs.empty());
      F rhsweight;
      F parentweight = (parent_weight[parent] -= weight);
      assert(parentweight >= 0);
      if (parentweight == 0) {
          parent_weight.erase(parent);
      }
      if (rhs.size() == 1) {
          S_F& parent1_weight = unarychild_parent_weight[rhs[0]];
          rhsweight = (parent1_weight[parent] -= weight);
          assert(rhsweight >= 0);
          if (rhsweight == 0) {
              parent1_weight.erase(parent);
              if (parent1_weight.empty()) {
                  unarychild_parent_weight.erase(rhs[0]);
              }
          }
      }
      else {  // non-unary rule
          S_F& parent1_weight = rhs_parent_weight[rhs];
          rhsweight = (parent1_weight[parent] -= weight);
          if (rhsweight == 0) {
              parent1_weight.erase(parent);
              if (parent1_weight.empty()) {
                  rhs_parent_weight.erase(rhs);
              }
          }
      }
      return rhsweight/parentweight;
  }  // pycfg_type::decrrule()

  //! IncreasePLUT2B()
  //! increment the part of rules from PLU_Top to PLU_Btm
  F IncreasePLUT2B(const Ss& TOPs, const Sss& BTMs) {
      F prob = 1;
      F weight = 1;
      assert(TOPs.size() == BTMs.size());
      for (size_t i = 0; i < TOPs.size(); ++i) {
          assert(BTMs[i].size()> 0);
          F ruleprob = incrrule(TOPs[i], BTMs[i], \
                  estimate_theta_flag * weight);
          if (debug == -1) {
              std::cerr << "Increase " << TOPs[i] 
                  << " --> " << BTMs[i] 
                  << ", prob = " << ruleprob << std::endl;
          }
          prob *= ruleprob;
      }
      return prob;
  }

  //! DecreasePLUT2B()
  //! increment the part of rules from PLU_Top to PLU_Btm
  F DecreasePLUT2B(const Ss& TOPs, const Sss& BTMs) {
      F prob = 1;
      F weight = 1;
      assert(TOPs.size() == BTMs.size());
      for (size_t i = 0; i < TOPs.size(); ++i) {
          assert(BTMs[i].size()> 0);
          F ruleprob = decrrule(TOPs[i], BTMs[i], \
                  estimate_theta_flag * weight);
          if (debug == -1) {
              std::cerr << "Decrease " << TOPs[i] 
                  << " --> " << BTMs[i] 
                  << ", prob = " << ruleprob << std::endl;
          }
          prob *= ruleprob;
      }
      return prob;
  }

  //! incrtree() increments the cache for tp, increments
  //! the rules if the cache count is appropriate, and returns
  //! the probability of this tree under the original model.
  //
  F incrtree(tree* tp) {
      U weight = 1; // added by Jackie
      if (tp -> children.empty())  { 
          return 1;  // terminal node
      }
      assert(weight >= 0);
      F pya = get_pya(tp->cat);    // PY cache statistics
      F pyb = get_pyb(tp->cat);
      if (pya == 1) { // don't table this category
          F prob = 1;
            {
                Ss children;
                cforeach (tree::ptrs_type, it, tp -> children)
                    children.push_back((*it)->cat);
                prob *= incrrule(tp->cat, children, estimate_theta_flag*weight);
            }
            cforeach (tree::ptrs_type, it, tp->children)
                prob *= incrtree(*it);
          return prob;
      }
      else if (tp->count > 0) {  // old PY table entry
          U& pyn = parent_pyn[tp->cat];
          F prob = (tp->count - pya)/(pyn + pyb);
          assert(finite(prob)); assert(prob > 0); assert(prob <= 1);
          tp->count += weight;              // increment entry count
          pyn += weight;                    // increment PY count
          return prob;
      } 
      else { // new PY table entry
          {
              Ss terms;
              tp->terminals(terms);
              bool inserted ATTRIBUTE_UNUSED = terms_pytrees[terms].insert(tp).second;
              assert(inserted);
          }
          U& pym = parent_pym[tp->cat];
          U& pyn = parent_pyn[tp->cat];
          F prob = (pym*pya + pyb)/(pyn + pyb);  // select new table
          assert(finite(prob)); assert(prob > 0); assert(prob <= 1);
          tp->count += weight;              // increment count
          pym += 1;                         // one more PY table entry
          pyn += weight;                    // increment PY count
          {
              Ss children;
              cforeach (tree::ptrs_type, it, tp->children)
                  children.push_back((*it)->cat);
              prob *= incrrule(tp->cat, children, estimate_theta_flag*weight);
          }
          cforeach (tree::ptrs_type, it, tp->children)
              prob *= incrtree(*it);
          return prob;
      }
  }  // pycfg_type::incrtree()

  // added by Jackie
  void retrieve_rule_counts(Ss& plu_top, trie<S, std::map<S, float> >& current_rule_counts, std::map<S, float>& current_parent_counts) {
      std::set<symbol> unique_plu(plu_top.begin(), plu_top.end());
      std::set<symbol>::iterator plu_iter = unique_plu.begin();
      for (; plu_iter != unique_plu.end(); ++plu_iter) {
          current_parent_counts[*plu_iter] = parent_weight[*plu_iter];
          Sss& rhss = rules_look_up[*plu_iter];
          for (size_t i = 0; i < rhss.size(); ++i) {
              if (rhss[i].size() == 1) {
                  S_F& parent_counts = unarychild_parent_weight[rhss[i][0]];
                  if (parent_counts.find(*plu_iter) != parent_counts.end()) {
                      float count = float(parent_counts[*plu_iter]);
                      current_rule_counts[rhss[i]][*plu_iter] = count;
                  }
              }
              else {
                  S_F& parent_counts = rhs_parent_weight.find(rhss[i]) -> data;
                  if (parent_counts.find(*plu_iter) != parent_counts.end()) {
                      float count = float(parent_counts[*plu_iter]); 
                      current_rule_counts[rhss[i]][*plu_iter] = count;
                  }
              }
          }
      }
  }

  // added by Jackie
  int retrieve_rule_weights(Ss& plu_top, trie<S, std::map<S, float> >& rule_weights_snapshot, float min = 0) {
      std::set<symbol> unique_plu(plu_top.begin(), plu_top.end());
      std::set<symbol>::iterator plu_iter = unique_plu.begin();
      int counter = 0;
      int eliminated_rule = 0;
      for (; plu_iter != unique_plu.end(); ++plu_iter) {
          Sss& rhss = rules_look_up[*plu_iter];
          for (size_t i = 0; i < rhss.size(); ++i) {
              if (rhss[i].size() == 1) {
                  S_F& parent_counts = unarychild_parent_weight[rhss[i][0]];
                  if (parent_counts.find(*plu_iter) != parent_counts.end()) {
                      float prob = float(parent_counts[*plu_iter]) / float(afind(parent_weight, *plu_iter));
                      if (prob >= min) {
                          ++counter;
                          prob = prob == 0 ? MIN_PROB_VALUE : log(prob);
                          rule_weights_snapshot[rhss[i]][*plu_iter] = prob;
                          // std::cout << "rule(" << *plu_iter << " --> " << rhss[i] << ") = " << prob << std::endl;
                      }
                      else {
                          ++eliminated_rule;
                      }
                  }
              }
              else {
                  S_F& parent_counts = rhs_parent_weight.find(rhss[i]) -> data;
                  if (parent_counts.find(*plu_iter) != parent_counts.end()) {
                      float prob = float(parent_counts[*plu_iter]) / float(afind(parent_weight, *plu_iter));
                      if (prob >= min) {
                          ++counter;
                          prob = prob == 0 ? MIN_PROB_VALUE : log(prob);
                          rule_weights_snapshot[rhss[i]][*plu_iter] = prob;
                          // std::cout << "rule(" << *plu_iter << " --> " << rhss[i] << ") = " << prob << std::endl;
                      }
                      else {
                          ++eliminated_rule;
                      }
                  }
              }
          }
      }
      // std::cerr << "# eliminated rules = " << eliminated_rule << endl;
      return counter;
  }

  float filter_rules(F thres) {
      unsigned int total_rule_num = 0;
      unsigned int eliminated_rule_num = 0;
      // unsigned int check_total_num = 0;
      // filter out unary rules
      filtered_unarychild_parent_weight.clear();
      cforeach(S_S_F, it, unarychild_parent_weight) {
          S child = it -> first;
          const S_F& rules = it -> second;
          cforeach(S_F, ptr, rules) {
              S parent = ptr -> first;
              F count = ptr -> second;
              F rule_prob = count / afind(parent_weight, parent);
              if (rule_prob >= thres) {
                  filtered_unarychild_parent_weight[child][parent] = count;
              }
              else {
                  // std::cerr << "eliminating " << parent << " --> " << child << " = " << rule_prob << std::endl; 
                  ++eliminated_rule_num;
              }
              ++total_rule_num;
          }
      }
      filtered_rhs_parent_weight.clear();
      rhs_parent_weight.for_each(EliminateBinaryRules((*this), total_rule_num, \
                                                      eliminated_rule_num, \
                                                      thres, \
                                                      filtered_rhs_parent_weight)); 
      filtered_terms_pytrees.clear();
      terms_pytrees.for_each(EliminateCachedTables((*this), total_rule_num, \
                                                   eliminated_rule_num, \
                                                   thres, filtered_terms_pytrees));
      return (float) eliminated_rule_num / total_rule_num;
  } 

  struct CountRuleNums {
      unsigned int& total_rule_num; 

      CountRuleNums(unsigned int& t)
         : total_rule_num(t) {}
      
      template <typename KeyType, typename DataType>
      void operator() (const KeyType& children, const DataType& rules) {
          cforeach (typename DataType, it, rules) {
              ++total_rule_num;
          }
      }  
  };  

  struct EliminateCachedTables {
      const pycfg_type& g;
      unsigned int& total_rule_num; 
      unsigned int& eliminated_rule_num;
      F thres;
      St_sT& filtered_terms_pytrees;

      EliminateCachedTables(pycfg_type& grammar, unsigned int& t, unsigned int& e, \
                            F p, St_sT& trees)
        : g(grammar), total_rule_num(t), eliminated_rule_num(e), \
          thres(p), filtered_terms_pytrees(trees) {
      }

      template <typename KeyType, typename Trees>
      void operator() (const KeyType& children, const Trees& trees) {
          cforeach (typename Trees, it, trees) {
              F pya = g.get_pya((*it) -> cat);
              F pyb = g.get_pyb((*it) -> cat);
              U pym = dfind(g.parent_pym, (*it) -> cat);
              U pyn = dfind(g.parent_pyn, (*it) -> cat);
              F rule_prob;

              if ((*it) -> count > 0) {
                  rule_prob = ((*it) -> count - pya) / (pyn + pyb);
              }
              else {
                 // std::cerr << "getting things from 0 " << std::endl;
                 rule_prob = (pym * pya + pyb) / (pyn + pyb);
              }
              if (rule_prob >= thres) {
                  filtered_terms_pytrees[children].insert((*it));
                  // std::cerr << "accepting " << (*it) -> cat 
                  //          << " --> " << children << " = " << rule_prob << std::endl; 
              }
              else {
                  // std::cerr << "eliminating " << (*it) -> cat 
                  //          << " --> " << children << " = " << rule_prob << std::endl; 
                  ++eliminated_rule_num;
              }
              ++total_rule_num;
          }
      }  
  };

  struct EliminateBinaryRules {
      const pycfg_type& g;
      unsigned int& total_rule_num; 
      unsigned int& eliminated_rule_num;
      F thres;
      St_S_F& filtered_rhs_parent_weight;

      EliminateBinaryRules(pycfg_type& grammar, \
                           unsigned int& t, \
                           unsigned int& e, \
                           F p, \
                           St_S_F& binary_rules)
         : g(grammar), total_rule_num(t), \
           eliminated_rule_num(e), \
           thres(p), \
           filtered_rhs_parent_weight(binary_rules) {
      }
      
      template <typename KeyType, typename DataType>
      void operator() (const KeyType& children, const DataType& rules) {
          cforeach (typename DataType, it, rules) {
              S parent = it -> first;
              F count = it -> second;
              F rule_prob = count / afind(g.parent_weight, parent);
              if (rule_prob >= thres) {
                  filtered_rhs_parent_weight[children][parent] = count;
                  // std::cerr << "accepting " << parent << " --> " 
                  //          << children << " = " << rule_prob << std::endl; 
              }
              else {
                  // std::cerr << "eliminating " << parent << " --> " 
                  //          << children << " = " << rule_prob << std::endl; 
                  ++eliminated_rule_num;
              }
              ++total_rule_num;
          }
      }  
  };  


  // added by Jackie
  void find_rule_prob_thres(tree* tp, Sss& plu_btms, \
          S& parent, Ss& rhs, F& min_prob) {
      // Check whether the non-terminal is cached
      F pya = get_pya(tp -> cat);
      if (pya == 1) {
          if (tp -> is_plu_top()) {
              Ss children = plu_btms[0]; 
              assert(children.size() > 0);
              plu_btms.erase(plu_btms.begin());
              F branch_prob = rule_prob(tp -> cat, children);
              // std::cerr << "branch_prob (" << tp -> cat << " --> " << children << " = " << branch_prob << std::endl;
              if (branch_prob < min_prob) {
                  min_prob = branch_prob;
                  parent = tp -> cat;
                  rhs = children;
              }
          }
          else {
              Ss children;
              cforeach(tree::ptrs_type, it, tp -> children) {
                  children.push_back((*it) -> cat);
              }
              F branch_prob = rule_prob(tp -> cat, children);
              if (branch_prob < min_prob) {
                  min_prob = branch_prob;
                  parent = tp -> cat;
                  rhs = children;
              }
              // std::cerr << "branch_prob (" << tp -> cat << " --> " << children << " = " << branch_prob << std::endl;
              cforeach(tree::ptrs_type, it, tp -> children) {
                  find_rule_prob_thres(*it, plu_btms, parent, rhs, min_prob);
              }
          }
      }
      else {
          F pyb = get_pyb(tp -> cat);
          U pym = dfind(parent_pym, tp -> cat);
          U pyn = dfind(parent_pyn, tp -> cat);
          if (tp -> count > 0) { // existing node
              assert(tp -> count <= pyn);
              assert(pym > 0);
              F branch_prob = (tp -> count - pya) / (pyn + pyb);
              Ss children;
              tp -> terminals(children);
              if (branch_prob < min_prob) {
                  min_prob = branch_prob;
                  parent = tp -> cat;
                  rhs = children;
              }
              // std::cerr << "branch_prob (" << tp -> cat << " --> " << children << " = " << branch_prob << std::endl;
              std::vector<catcounttree_type*> leaf_nodes;
              tp -> terminal_ptrs(leaf_nodes);
              for (size_t i = 0; i < leaf_nodes.size(); ++i) {
                  find_rule_prob_thres(leaf_nodes[i], plu_btms, parent, rhs, min_prob);
              }
          }
          else {
              // std::cerr << "0 scenario" << std::endl;
              F branch_prob = (pym * pya + pyb) / (pyn + pyb);
              Ss children;
              cforeach(tree::ptrs_type, it, tp -> children) {
                  children.push_back((*it) -> cat);
              }
              branch_prob *= rule_prob(tp -> cat, children);
              if (branch_prob < min_prob) {
                  min_prob = branch_prob;
                  parent = tp -> cat;
                  rhs = children;
              }
              cforeach(tree::ptrs_type, it, tp -> children) {
                  find_rule_prob_thres((*it), plu_btms, parent, rhs, min_prob);
              }
          }
      }
  } 

  // added by Jackie
  // find the lowest probable rule for plu_top -> plu_btm
  float find_rule_prob_thres(Ss& plu_top, Sss& plu_btm, S& parent, Ss& rhs) {
      if (plu_top.size() != plu_btm.size()) {
          std::cerr << "plu top size != plu btm size in find_min_prob_rule" << std::endl;
      }
      float min_prob = 100000;
      for (size_t i = 0; i < plu_top.size(); ++i) {
          if (plu_btm[i].size() == 1) {
              float rule_prob = \
                float(unarychild_parent_weight[plu_btm[i][0]][plu_top[i]]) / float(afind(parent_weight, plu_top[i])) ;
              if (rule_prob < min_prob) {
                  min_prob = rule_prob;
                  parent = plu_top[i];
                  rhs = plu_btm[i];
              }
          }
          else {
              float rule_prob = \
                float((rhs_parent_weight.find(plu_btm[i]) -> data)[plu_top[i]]) / float(afind(parent_weight, plu_top[i]));
              if (rule_prob < min_prob) {
                  min_prob = rule_prob;
                  parent = plu_top[i];
                  rhs = plu_btm[i];
              }
          }
      }
      // min_prob *= float(random1());
      assert(min_prob >= 0);
      assert(min_prob <= 1);
      return min_prob;
  }

  //! decrtree() decrements the cache for tp, decrements
  //! the rules if the cache count is appropriate, and returns
  //! the probability of this tree under the new model.
  //
  F decrtree(tree* tp) {
      U weight = 1; // Modified by Jackie
      if (tp -> children.empty())  {
          return 1;  // terminal node
      }
      F pya = get_pya(tp->cat);    // PY cache statistics
      if (pya == 1) {  // don't table this category
          F prob = 1;
          {
              Ss children;
              cforeach (tree::ptrs_type, it, tp->children) {
                  children.push_back((*it)->cat);
              }
              F ruleprob = decrrule(tp->cat, children, estimate_theta_flag*weight);
              assert(ruleprob > 0);
              prob *= ruleprob;
          }
          cforeach (tree::ptrs_type, it, tp->children)  {
              prob *= decrtree(*it);
          }
          return prob;
      }
      assert(weight <= tp->count);
      tp -> count -= weight;
      assert(afind(parent_pyn, tp->cat) >= weight);
      const U pyn = (parent_pyn[tp->cat] -= weight);
      F pyb = get_pyb(tp->cat);
      if (tp->count > 0) {  // old PY table entry
          assert(pyn > 0);
          F prob = (tp->count - pya)/(pyn + pyb);
          assert(finite(prob)); assert(prob > 0); assert(prob <= 1);
          return prob;
      } 
      else { // tp->count == 0, remove PY table entry
          {
              Ss terms;
              tp->terminals(terms);
              sT& pytrees = terms_pytrees[terms];
              sT::size_type nerased ATTRIBUTE_UNUSED = pytrees.erase(tp);
              assert(nerased == 1);
              if (pytrees.empty()) {
                  terms_pytrees.erase(terms);
              }
          }
          // Bug: when pym or pyn goes to zero and the parent is erased, 
          // and then the reference to pym or pyn becomes a dangling reference
          // U& pym = parent_pym[tp->cat];
          // pym -= 1;                         // reduce cache count
          assert(parent_pym.count(tp->cat) > 0);
          const U pym = --parent_pym[tp->cat];
          if (pym == 0) {
              parent_pym.erase(tp->cat);
          }
          if (pyn == 0) {
              parent_pyn.erase(tp->cat);
          }
          F prob = (pym*pya + pyb)/(pyn + pyb);  // select new table
          assert(finite(prob)); assert(prob > 0); assert(prob <= 1);
          {
              Ss children;
              cforeach (tree::ptrs_type, it, tp->children)
                  children.push_back((*it)->cat);
              prob *= decrrule(tp->cat, children, estimate_theta_flag*weight);
          }
          assert(prob > 0);
          cforeach (tree::ptrs_type, it, tp->children) {
              prob *= decrtree(*it);
          }
          // assert(prob > 0);
          return prob;
      }
  }  // pycfg_type::decrtree()

  //! read() reads a grammar from an input stream (implements >> )
  //
  std::istream& read(std::istream& is) {
    start = symbol::undefined();
    F weight = default_weight;
    F pya = default_pya;
    F pyb = default_pyb;
    S parent;
    while (is >> default_value(weight, default_weight) 
	      >> default_value(pya, default_pya)
	      >> default_value(pyb, default_pyb)
	      >> parent >> " -->") {
      if (weight<=0)
          weight=default_weight;
      if (start.is_undefined())
          start = parent;
      Ss rhs;
      readline_symbols(is, rhs);
      rules_look_up[parent].push_back(rhs);
      if (debug >= 100000)
          std::cerr << "# " << weight << '\t' << parent << " --> " << rhs << std::endl;

      incrrule(parent, rhs, weight);
      if (pya != default_pya)
          parent_pya[parent] = pya;
      if (pyb != default_pyb)
          parent_pyb[parent] = pyb;
      rule_priorweight[SSs(parent,rhs)] += weight;
      parent_priorweight[parent] += weight;
    }
    return is;
  }  // pycfg_type::read()

  //! write() writes a grammar (implements << )
  //
  std::ostream& write(std::ostream& os) const {
    assert(start.is_defined());
    write_rules(os, start);
    cforeach (S_F, it, parent_weight)
      if (it->first != start) 
	write_rules(os, it->first);
    return os;
  }  // pycfg_type::write()

  std::ostream& write_rules(std::ostream& os, S parent) const {
    rhs_parent_weight.for_each(write_rule(os, parent));
    cforeach (S_S_F, it0, unarychild_parent_weight) {
      S child = it0->first;
      const S_F& parent_weight = it0->second;
      cforeach (S_F, it1, parent_weight)
	if (it1->first == parent)
	  os << it1->second << '\t' << parent 
	     << " --> " << child << std::endl;
    }
    bool old_compact_trees_flag = compact_trees;  // save old flag
    compact_trees = false;  // turn off compact_trees
    terms_pytrees.for_each(write_pycache(os, parent));
    compact_trees = old_compact_trees_flag;
    return os;
  }  // pycfg_type::write_rules()

  //! write_rule{} writes a single rule
  //
  struct write_rule {
    std::ostream& os;
    S parent;

    write_rule(std::ostream& os, symbol parent) : os(os), parent(parent) { }

    template <typename Keys, typename Value>
    void operator() (const Keys& rhs, const Value& parentweights) {
      cforeach (typename Value, pwit, parentweights) 
	if (pwit->first == parent) {
	  os << pwit->second << '\t' << parent << " -->";
	  cforeach (typename Keys, rhsit, rhs)
	    os << ' ' << *rhsit;
	  os << std::endl;
	}
    }  // pycfg_type::write_rule::operator()

  };  // pycfg_type::write_rule{}
  
  //! write_pycache{} writes the cache entries for a category
  //
  struct write_pycache {
    std::ostream& os;
    S parent;
    
    write_pycache(std::ostream& os, S parent) : os(os), parent(parent) { }

    template <typename Words, typename TreePtrs>
    void operator() (const Words& words, const TreePtrs& tps) {
      cforeach (typename TreePtrs, tpit, tps) 
	if ((*tpit)->cat == parent)
	  os << (*tpit) << std::endl;
    }  // pycfg_type::write_pycache::operator()
  };  // pycfg_type::write_pycache{}

  //! logPcorpus() returns the log probability of the corpus trees
  //
  F logPcorpus() const {
      F logP = 0;
      // grammar part
      cforeach (SSs_F, it, rule_priorweight) {
          S parent = it->first.first;
          const Ss& rhs = it->first.second;
          F priorweight = it->second;
          F weight = rule_weight(parent, rhs);
          logP += lgamma(weight) - lgamma(priorweight);
      }
      if (debug >= 5000)
          TRACE1(logP);
      cforeach (S_F, it, parent_priorweight) {
          S parent = it->first;
          F priorweight = it->second;
          F weight =dfind(parent_weight, parent);
          logP += lgamma(priorweight) - lgamma(weight);
      }
      if (debug >= 5000)
          TRACE1(logP);
      assert(logP <= 0);
      // PY adaptor part
      cforeach (S_U, it, parent_pyn) {
          S parent = it->first;
          U pyn = it->second;
          U pym = afind(parent_pym, parent);
          F pya = get_pya(parent);
          F pyb = get_pyb(parent);
          logP += lgamma(pyb) - lgamma(pyn+pyb);
          for (U i = 0; i < pym; ++i)
              logP += log(i*pya + pyb);
      }
      if (debug >= 5000)
          TRACE1(logP);
      terms_pytrees.for_each(logPcache(*this, logP));
      if (debug >= 5000)
          TRACE1(logP);
      assert(logP <= 0);
      return logP;
  }  // pycfg_type::logPcorpus()

  struct logPcache {
      const pycfg_type& g;
      F& logP;

      logPcache(const pycfg_type& g, F& logP) : g(g), logP(logP) { }
      
      template <typename Words, typename TreePtrs>
          void operator() (const Words& words, const TreePtrs& tps) {
              cforeach (typename TreePtrs, it, tps) {
                  S parent = (*it)->cat;
                  U count = (*it)->count;
                  F pya = g.get_pya(parent);
                  logP += lgamma(count-pya) - lgamma(1-pya);
              }
          }  // pycfg_type::logPcache::operator()
  };  // pycfg_type::logPcache{}

  //! logPrior() returns the prior probability of the PY a and b values
  //
  F logPrior() const {
    F sumLogP = 0;
    if (pyb_gamma_s > 0 && pyb_gamma_c > 0)
      cforeach (S_U, it, parent_pyn) {
	S parent = it->first;
	F pya = get_pya(parent);
	assert(pya >= 0);
	assert(pya <= 1);
	F pyb = get_pyb(parent);
	assert(pyb >= 0);
	if (pya_beta_a > 0 && pya_beta_b > 0 && pya > 0) {
	  F logP = pya_logPrior(pya, pya_beta_a, pya_beta_b);
	  if (debug >= 2000)
	    TRACE5(parent, logP, pya, pya_beta_a, pya_beta_b);
	  sumLogP += logP;
	}
	F logP = pyb_logPrior(pyb, pyb_gamma_c, pyb_gamma_s);
	if (debug >= 2000)
	  TRACE5(parent, logP, pyb, pyb_gamma_c, pyb_gamma_s);
	sumLogP += logP;
      }
    return sumLogP;
  }  // pycfg_type::logPrior()

  //! pya_logPrior() calculates the Beta prior on pya.
  //
  static F pya_logPrior(F pya, F pya_beta_a, F pya_beta_b) {
    F prior = lbetadist(pya, pya_beta_a, pya_beta_b);     //!< prior for pya
    return prior;
  }  // pycfg_type::pya_logPrior()

  //! pyb_logPrior() calculates the prior probability of pyb 
  //! wrt the Gamma prior on pyb.
  //
  static F pyb_logPrior(F pyb, F pyb_gamma_c, F pyb_gamma_s) {
    F prior = lgammadist(pyb, pyb_gamma_c, pyb_gamma_s);  // prior for pyb
    return prior;
  }  // pcfg_type::pyb_logPrior()

  //////////////////////////////////////////////////////////////////////
  //                                                                  //
  //                      Resample pyb                                //
  //                                                                  //
  //////////////////////////////////////////////////////////////////////
  
  //! resample_pyb_type{} is a function object that returns the part of log prob that depends on pyb.
  //! This includes the Gamma prior on pyb, but doesn't include e.g. the rule probabilities
  //! (as these are a constant factor)
  //
  struct resample_pyb_type {
    typedef double F;

    U pyn, pym;
    F pya, pyb_gamma_c, pyb_gamma_s;
    resample_pyb_type(U pyn, U pym, F pya, F pyb_gamma_c, F pyb_gamma_s) 
      : pyn(pyn), pym(pym), pya(pya), pyb_gamma_c(pyb_gamma_c), pyb_gamma_s(pyb_gamma_s)
    { }

    //! operator() returns the part of the log posterior probability that depends on pyb
    //
    F operator() (F pyb) const {
      assert(pyb > 0);
      F logPrior = pyb_logPrior(pyb, pyb_gamma_c, pyb_gamma_s);  //!< prior for pyb
      F logProb = 0;
      logProb += (pya == 0 ? pym*log(pyb) : pym*log(pya) + lgamma(pym + pyb/pya) - lgamma(pyb/pya));
      logProb += lgamma(pyb) - lgamma(pyn+pyb);
      return logProb+logPrior;
    }
  };  // pcfg_type::resample_pyb_type{}

  //! resample_pyb() samples new values for pyb for each adapted nonterminal.
  //! It returns the log prior prob of new b values.
  //
  void resample_pyb() {
    U niterations = 20;   //!< number of resampling iterations
    // std::cerr << "\n## resample_pyb(), initial parent_pya = " << parent_pya << ", parent_pyb = " << parent_pyb << std::endl;
    cforeach (S_U, it, parent_pyn) {
      S parent = it->first;
      U pyn = it->second;
      U pym = afind(parent_pym, parent);
      F pya = get_pya(parent);
      F pyb = get_pyb(parent);
      resample_pyb_type pyb_logP(pyn, pym, pya, pyb_gamma_c, pyb_gamma_s);
      // pyb = slice_sampler1d(pyb_logP, pyb, random1, 0.0, std::numeric_limits<double>::infinity(), 0.0, niterations, 100*niterations);
      pyb = slice_sampler1dp(pyb_logP, pyb, random1, 1, niterations);
      parent_pyb[parent] = pyb;
      // parent_bap[parent].first += naccepted;
      // parent_bap[parent].second += nproposed;
    }
  }  // pcfg_type::resample_pyb()

  //////////////////////////////////////////////////////////////////////
  //                                                                  //
  //                   Resample pya and pyb                           //
  //                                                                  //
  //////////////////////////////////////////////////////////////////////

  //! resample_pya_type{} calculates the part of the log prob that depends on pya.
  //! This includes the Beta prior on pya, but doesn't include e.g. the rule probabilities
  //! (as these are a constant factor)
  //
  struct resample_pya_type {
    U pyn, pym;
    F pyb, pya_beta_a, pya_beta_b;
    const Ts& trees;
    
    resample_pya_type(U pyn, U pym, F pyb, F pya_beta_a, F pya_beta_b, const Ts& trees) 
      : pyn(pyn), pym(pym), pyb(pyb), pya_beta_a(pya_beta_a), pya_beta_b(pya_beta_b), trees(trees)
    { }

    //! operator() returns the part of the log posterior probability that depends on pya
    //
    F operator() (F pya) const {
      F logPrior = pya_logPrior(pya, pya_beta_a, pya_beta_b);     //!< prior for pya
      F logProb = 0;
      F lgamma1a = lgamma(1-pya);
      cforeach (Ts, it, trees) {
	U count = (*it)->count;
	logProb += lgamma(count-pya) - lgamma1a;
      }
      logProb += (pya == 0 ? pym*log(pyb) : pym*log(pya) + lgamma(pym + pyb/pya) - lgamma(pyb/pya));
      return logPrior + logProb;
    }   // pycfg_type::resample_pya_type::operator()

  };  // pycfg_type::resample_pya_type{}
  
  //! resample_pya() samples new values for pya for each adapted nonterminal
  //
  void resample_pya(const S_Ts& parent_trees) {
    U niterations = 20;   //!< number of resampling iterations
    // std::cerr << "\n## Initial parent_pya = " << parent_pya << ", parent_pyb = " << parent_pyb << std::endl;
    cforeach (S_U, it, parent_pyn) {
      S parent = it->first;
      F pya = get_pya(parent);
      if (pya == 0)   // if this nonterminal has pya == 0, then don't resample
	continue;
      F pyb = get_pyb(parent);
      U pyn = it->second;
      U pym = afind(parent_pym, parent);
      const Ts& trees = afind(parent_trees, parent);
      resample_pya_type pya_logP(pyn, pym, pyb, pya_beta_a, pya_beta_b, trees);
      pya = slice_sampler1d(pya_logP, pya, random1, std::numeric_limits<double>::min(), 1.0, 0.0, niterations);
      parent_pya[parent] = pya;
    }
  }  // pycfg_type::resample_pya()

  //! resample_pyab_parent_trees_helper{} constructs parent_trees from terms_pytrees.
  //
  struct resample_pyab_parent_trees_helper {
    S_Ts& parent_trees;
    resample_pyab_parent_trees_helper(S_Ts& parent_trees) : parent_trees(parent_trees) { }

    template <typename Words, typename TreePtrs>
    void operator() (const Words& words, const TreePtrs& tps) {
      cforeach (typename TreePtrs, it, tps) {
	S parent = (*it)->cat;
	parent_trees[parent].push_back(*it);
      }
    }  // pycfg_type::resample_pyab_parent_trees_helper::operator()
  };  // pycfg_type::resample_pyab_parent_trees_helper{}

  //! resample_pyab() resamples both pya and pyb for each adapted nonterminal.
  //
  void resample_pyab() {
    const U niterations = 5;  //!< number of alternating samples of pya and pyb
    S_Ts parent_trees;
    terms_pytrees.for_each(resample_pyab_parent_trees_helper(parent_trees));
    for (U i=0; i<niterations; ++i) {
      resample_pyb();
      resample_pya(parent_trees);
    }
    resample_pyb();
  }  // pycfg_type::resample_pyab()

  //! write_adaptor_parameters() writes out adaptor parameters to a file
  //
  std::ostream& write_adaptor_parameters(std::ostream& os) const {
    cforeach (S_F, it, parent_priorweight) {
      S parent = it->first;
      F pya = get_pya(parent);
      if (pya == 1)
	continue;
      U pym = dfind(parent_pym, parent);
      U pyn = dfind(parent_pyn, parent);
      F pyb = get_pyb(parent);
      os << ' ' << parent << ' ' << pym << ' ' << pyn << ' ' << pya << ' ' << pyb;
    }
    return os;
  }  // pycfg_type::write_adaptor_parameters()

  //! initialize_predictive_parse_filter() initializes the predictive
  //! parse filter by building the grammar that the Earley parser requires
  //
  void initialize_predictive_parse_filter() {
    predictive_parse_filter = true;
    cforeach (SSs_F, it, rule_priorweight) {
      const SSs& rule = it->first;
      const Ss& children = rule.second;
      assert(!children.empty());
      S child1 = children.front();
      predictive_parse_filter_grammar.add_rule(it->first, 
					       children.size() == 1 
					       && !parent_priorweight.count(child1));
    }
  }  // pycfg_type::initialize_predictive_parse_filter();

};  // pycfg_type{}

//! operator>> (pycfg_type&) reads a pycfg_type g, setting g.start
//! to the parent of the first rule read.
//
inline
std::istream& operator>> (std::istream& is, pycfg_type& g) {
  return g.read(is);
}  // operator>> (pycfg_type&)

inline
std::ostream& operator<< (std::ostream& os, const pycfg_type& g) {
  return g.write(os);
}  // operator<< (pycfg_type&)

namespace std { namespace tr1 {
    template <> struct hash<pycfg_type::Stit> 
      : public std::unary_function<pycfg_type::Stit, std::size_t> {
      size_t operator()(const pycfg_type::Stit t) const
      {
	return size_t(&(*t));
      }  // ext::hash<pycfg_type::Stit>::operator()
    };  // ext::hash<pycfg_type::Stit>{}
  }  } // namespace std::tr1
  
static const F unaryclosetolerance = 1e-7;

class pycky {

public:

  const pycfg_type& g;
  F anneal;         // annealing factor (1 = no annealing)
  
  pycky(const pycfg_type& g, F anneal=1) : g(g), anneal(anneal) { }

  typedef pycfg_type::tree tree;
  typedef pycfg_type::U U;
  typedef pycfg_type::S_S_F S_S_F;
  typedef pycfg_type::St_S_F St_S_F;
  typedef pycfg_type::Stit Stit;

  typedef std::vector<S_F> S_Fs;
  typedef tr1::unordered_map<Stit,F> Stit_F;
  typedef std::vector<Stit_F> Stit_Fs;

  typedef pycfg_type::sT sT;

  typedef pycfg_type::St_sT St_sT;
  typedef St_sT::const_iterator StsTit;
  typedef std::map<StsTit, F> StsTit_F; // modified by Jackie
  typedef std::vector<StsTit> StsTits;
  typedef std::set<S> sS;


  //! index() returns the location of cell in cells[]
  //
  static U index(U i, U j) { return j*(j-1)/2+i; }

  //! ncells() returns the number of cells required for sentence of length n
  //
  static U ncells(U n) { return n*(n+1)/2; }
  
  Ss terminals;
  S_Fs inactives;
  Stit_Fs actives;
  std::vector<StsTit_F> tree_leaf_probs; // modified by Jackie
  std::vector<sS> top_terminals; // added by Jackie
  StsTits pytits;

  typedef std::vector<sS> sSs;
  sSs predicteds;

  // added by Jackie
  //! get_top_terminals() return plu_tops for plu_bottoms 
  void get_top_terminals(const S_F& expansion, sS& top_terminals); 

  // added by Jackie
  //! get_top_terminals() return plu_tops for plu_bottoms 
  void get_top_terminals(const S& plu_btm, \
                         sS& top_terminals, float min_prob); 

  // added by Jackie
  //! get_top_terminals() return plu_tops for plu_bottoms 
  void get_top_terminals(const S& plu_btm, sS& top_terminals, S_S_F&); 

  // added by Jackie
  //! get_trees() return plu_tops for plu_bottoms 
  void get_trees(const sS& top_terminals, StsTit_F& trees, \
                 S_F& probs, St_sT&);

  // added by Jackie
  // ! get_expansion_prob() returns prob(plu_top -> plu_btm) for a terminal
  void get_expansion_prob(const S& terminal, S_F& prob, S_S_F&);

  // added by Jackie 
  //! extend_inside()
  /*
  template <typename terminals_type>
  F extend_inside(const terminals_type& terminals, float min_prob) {
      return extend_inside(terminals, g.start, min_prob);
  }
  */

  // added by Jackie 
  //! extend_inside()
  template <typename terminals_type>
  F extend_inside(const terminals_type& terminals, \
          S_S_F& unary_rules, \
          St_S_F& rhs_rules, \
          St_sT& terms) {
      return extend_inside(terminals, g.start, unary_rules, rhs_rules, terms);
  }

  // added by Jackie
  //! extend_inside() constructs the inside table, and returns the probability
  //! of the start symbol rewriting to the terminals.
  template <typename terminals_type>
  F extend_inside(const terminals_type& terminals0, S start, \
          S_S_F& unary_rules, St_S_F& rhs_rules, St_sT& terms) {
      terminals = terminals0;
      if (debug >= 10000) {
          std::cerr << "# cky::extend_inside() terminals = " << terminals << std::endl;
      }
      U n = terminals.size();

      top_terminals.clear();  // to store top terminals and to get tables
      top_terminals.resize(ncells(n));
      inactives.clear(); // for unary rules 
      inactives.resize(ncells(n));
      actives.clear(); // for multinary rules 
      actives.resize(ncells(n));
      tree_leaf_probs.clear(); // to store prob(top -> btm)
      tree_leaf_probs.resize(ncells(n));

#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (U i = 0; i < n; ++i) {   // terminals
          inactives[index(i, i + 1)][terminals[i]] = 1;
          // get prob(plu_top -> plu_btm) for the terminals
          S_F terminal_expansion_prob;
          get_expansion_prob(terminals[i], terminal_expansion_prob, unary_rules);

          // get top_terminals
          get_top_terminals(terminals[i], top_terminals[index(i, i + 1)], unary_rules);

          // get trees that expand to top_terminals
          StsTit_F& tree_leaf_prob = tree_leaf_probs[index(i, i + 1)];
          get_trees(top_terminals[index(i, i + 1)], tree_leaf_prob, terminal_expansion_prob, terms);
           
          StsTit_F::const_iterator tree_iter = tree_leaf_prob.begin(); 
          for (; tree_iter != tree_leaf_prob.end(); ++tree_iter) {
              if ((tree_iter -> first) != terms.end()) {
                  add_pycache((tree_iter -> first) -> data, \
                          inactives[index(i, i + 1)], tree_iter -> second);
              }  
          }
        
          inside_unaryclose(inactives[index(i, i + 1)], actives[index(i, i + 1)], NULL, unary_rules, rhs_rules);

          if (debug >= 20000) {
              std::cerr << "# cky::extend_inside() inactives[" << i << "," << i+1 << "] = " 
                        << inactives[index(i,i+1)] << std::endl;
          }
          if (debug >= 20100) {
              std::cerr << "# cky::extend_inside() actives[" << i << "," << i+1 << "] = " 
                << actives[index(i,i+1)] << std::endl;
          }
      }
      /*
      for (U i = 0; i < n ; ++i) {
          std::cerr << "termianl[" << i << "] = " << terminals[i] << std::endl;
          std::cerr << "# top_terminals[" << index(i, i + 1) << "] = " << top_terminals[index(i, i + 1)].size() << std::endl;
      }
      */

      for (U gap = 2; gap <= n; ++gap) { // non-terminals 
#ifdef _OPENMP
#pragma omp parallel for 
#endif
          for (U left = 0; left <= n - gap; ++left) {
              U right = left + gap;
              S_F& parentinactives = inactives[index(left, right)];
              Stit_F& parentactives = actives[index(left, right)];

              for (U mid = left + 1; mid < right; ++mid) {
                  const S_F& rightinactives = inactives[index(mid, right)];
                  if (rightinactives.empty()) {
                      continue;
                  }
                  Stit_F& leftactives = actives[index(left, mid)];
                  cforeach (Stit_F, itleft, leftactives) {
                      const Stit leftactive = itleft -> first;
                      const F leftprob = itleft -> second;
                      cforeach (S_F, itright, rightinactives) {
                          S rightinactive = itright -> first;
                          const F rightprob = itright -> second;
                          const Stit parentactive = leftactive -> find1(rightinactive);
                          if (parentactive != leftactive -> end()) {
                              F leftrightprob = leftprob * rightprob;
                              cforeach (S_F, itparent, parentactive -> data) {
                                  S parent = itparent -> first;
                                  if (parent.string_reference().find("T_") != std::string::npos) {
                                      F rule_prob = (itparent -> second) / afind(g.parent_weight, parent);
                                      parentinactives[parent] += leftrightprob * power(rule_prob, anneal);
                                  }
                                  else {
                                      parentinactives[parent] += leftrightprob * \
                                        power(itparent -> second / afind(g.parent_weight, parent), anneal);
                                  }
                              }
                              if (!parentactive -> key_trie.empty())
                                  parentactives[parentactive] += leftrightprob;
                          }
                      }
                  }
              } 
              // PY correction
              foreach (S_F, it, parentinactives) {
                  F pya = g.get_pya(it -> first);    // PY cache statistics
                  if (pya == 1.0) {
                      continue;
                  }
                  F pyb = g.get_pyb(it -> first);
                  U pym = dfind(g.parent_pym, it -> first);
                  U pyn = dfind(g.parent_pyn, it -> first);
                  it -> second *= power((pym * pya + pyb) / (pyn + pyb), anneal);
              }

              StsTit_F& tree_leaf_prob = tree_leaf_probs[index(left, right)];
              if (gap <= 2) {
                  get_top_terminals(parentinactives, top_terminals[index(left, right)]);
                  get_trees(top_terminals[index(left, right)], tree_leaf_prob, parentinactives, terms);
              }

              // need to get all the py trees for this cell
              // get top_terminals that map to one terminal and two terminals
              
              StsTit_F temp;
              for (U mid = right - 1 ; mid >= right - 2 && mid > left ; --mid) {
                  const StsTit_F& tables = tree_leaf_probs[index(left, mid)];
                  const S_F& plu_btm_probs = inactives[index(mid, right)];
                  StsTit_F::const_iterator t_iter = tables.begin();
                  for (; t_iter != tables.end(); ++t_iter) {
                      if ((t_iter -> first) != terms.end()) {
                          const std::map<S, trie<S, sT> >& next_tables = (t_iter -> first) -> key_trie;
                          std::map<S, trie<S, sT> >::const_iterator next_ptr = next_tables.begin();
                          for (; next_ptr != next_tables.end(); ++next_ptr) {
                              if (&(next_ptr -> second) != terms.end()) {
                                  F tree_expand_prob = (t_iter -> second) * dfind(plu_btm_probs, next_ptr -> first);
                                  if (tree_expand_prob > 0) {
                                      tree_leaf_prob[&(next_ptr -> second)] +=  tree_expand_prob;
                                  }
                              }
                          }
                      }
                  }
              }

              StsTit_F::const_iterator iter = tree_leaf_prob.begin();

              for (; iter != tree_leaf_prob.end(); ++iter) {
                  if ((iter -> first) != terms.end()) {
                      add_pycache((iter -> first) -> data, parentinactives, iter -> second);
                  }
              }
              inside_unaryclose(parentinactives, parentactives, NULL, unary_rules, rhs_rules);

              if (debug >= 20000) {
                  std::cerr << "# cky::extend_inside() inactives[" << left << "," << right 
                      << "] = " << parentinactives << std::endl;
              }
              if (debug >= 20100)
                  std::cerr << "# cky::extend_inside() actives[" << left << "," << right << "] = " 
                            << parentactives << std::endl;
          }
      }
      if (inactives[index(0, n)].find(start) == inactives[index(0, n)].end()) {
          std::cerr << "Can't find anything" << std::endl;
      }
      return dfind(inactives[index(0, n)], start);
  }  // pycky::extend_inside()

  // added by Jackie
  //! add_pycache() 
  void add_pycache(const sT& tps, S_F& inactives, const F prob) const {
      cforeach (sT, it, tps) {
          symbol cat = (*it) -> cat;
          F pya = g.get_pya(cat);    // PY cache statistics
          if (pya == 1.0)
              continue;
          F pyb = g.get_pyb(cat);
          U pyn = dfind(g.parent_pyn, cat);
          inactives[cat] += power( ((*it)->count - pya)/(pyn + pyb), anneal) * prob;
      }
  }  // pycky::add_cache()

  void add_pycache(const sT& tps, S_F& inactives) const {
      cforeach (sT, it, tps) {
          symbol cat = (*it) -> cat;
          F pya = g.get_pya(cat);    // PY cache statistics
          if (pya == 1.0)
              continue;
          F pyb = g.get_pyb(cat);
          U pyn = dfind(g.parent_pyn, cat);
          inactives[cat] += power(((*it) -> count - pya) / (pyn + pyb), anneal);
      }
  }  // pycky::add_cache()

  void inside_unaryclose(S_F& inactives, \
                         Stit_F& actives, \
                         const sS* predictedparents, \
                         const S_S_F& unary_rules, \
                         const St_S_F& rhs_rules) const {
      F delta = 1;
      S_F delta_prob1 = inactives;
      S_F delta_prob0;
      while (delta > unaryclosetolerance) {
          delta = 0;
          delta_prob0 = delta_prob1;
          delta_prob1.clear();
          cforeach (S_F, it0, delta_prob0) {
              S child = it0 -> first;
              S_S_F::const_iterator it = unary_rules.find(child);
              if (it != unary_rules.end()) {
                  const S_F& parent_count = it -> second;
                  cforeach (S_F, it1, parent_count) {
                      S parent = it1 -> first;
                      if (parent.string_reference().find("T_") != std::string::npos) {
                          F prob = it0 -> second;
                          F pya = g.get_pya(parent);
                          F expansion_prob;
                          if (pya == 1) {
                              expansion_prob = it1 -> second / afind(g.parent_weight, parent);
                          }
                          else {
                              F pyb = g.get_pyb(parent);
                              U pym = dfind(g.parent_pym, parent);
                              U pyn = dfind(g.parent_pyn, parent);
                              expansion_prob = (it1 -> second) / afind(g.parent_weight, parent) \
                                               * (pym * pya + pyb) / (pyn + pyb);
                          }
                          prob *= power(expansion_prob, anneal); 
                          delta_prob1[parent] += prob;
                          delta = std::max(delta, prob/(inactives[parent] += prob));
                      }
                      else {
                          F prob = it0 -> second; // child's prob so far
                          F pya = g.get_pya(parent);
                          if (pya == 1)
                              prob *= power(it1 -> second/afind(g.parent_weight, parent), anneal); // prob(parent -> child)
                          else {
                              F pyb = g.get_pyb(parent);
                              U pym = dfind(g.parent_pym, parent);
                              U pyn = dfind(g.parent_pyn, parent);
                              prob *= power(it1 -> second / afind(g.parent_weight, parent) * (pym * pya + pyb) / (pyn + pyb), anneal);
                          }
                          delta_prob1[parent] += prob;
                          delta = std::max(delta, prob/(inactives[parent] += prob));
                      }
                  }
              }
          }
      }
      cforeach (S_F, it0, inactives) {
          Stit it1 = rhs_rules.find1(it0 -> first);
          if (it1 != rhs_rules.end())
              actives[it1] += it0 -> second;
      }
  } // pycky::inside_unaryclose()

  // added by Jackie
  void random_plu(Sss&, StsTit, const F prob, U left, U right, St_sT& terms);

  // added by Jackie
  tree* extend_random_inactive(const S parent, F parentprob, 
	 const U left, const U right, Sss& new_terminals, \
     S_S_F& unary_rules, St_S_F& rhs_rules, St_sT& terms);

  // added by Jackie
  void extend_random_active(const Stit parent, F parentprob, const U left, const U right, \
          tree::ptrs_type& siblings, Sss& new_terminals, \
          S_S_F& unary_rules, St_S_F& rhs_rules, St_sT& terms);

  // added by Jackie
  tree* extend_random_tree(S s, Sss& new_terminals, S_S_F& unary_rules, St_S_F& rhs_rules, St_sT& terms) {
      U n = terminals.size();
      return extend_random_inactive(s, afind(inactives[index(0, n)], s), 0, n, new_terminals, unary_rules, rhs_rules, terms);
  } 

  // added by Jackie
  tree* extend_random_tree(Sss& new_plu_btm, \
          S_S_F& unary_rules, \
          St_S_F& rhs_rules, \
          St_sT& terms) { 
      return extend_random_tree(g.start, new_plu_btm, unary_rules, rhs_rules, terms); }

  //! inside() constructs the inside table, and returns the probability
  //! of the start symbol rewriting to the terminals.
  //
  template <typename terminals_type>
  F inside(const terminals_type& terminals0, S start) {
      terminals = terminals0;
      if (debug >= 10000)
          std::cerr << "# cky::inside() terminals = " << terminals << std::endl;
      
      U n = terminals.size();
      if (g.predictive_parse_filter) {
          earley(g.predictive_parse_filter_grammar, start, terminals, predicteds);
          if (!predicteds[index(0,n)].count(start)) 
              std::cerr << "## " << HERE << " Error: earley parse failed, terminals = " 
                  << terminals << std::endl << exit_failure;
      }
      
      inactives.clear();
      inactives.resize(ncells(n));
      actives.clear();
      actives.resize(ncells(n));
      pytits.clear();
      pytits.resize(ncells(n));

#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (U i = 0; i < n; ++i) {   // terminals
          pytits[index(i, i+1)] = g.terms_pytrees.find1(terminals[i]);  // PY cache
          inactives[index(i,i+1)][terminals[i]] = 1;
          StsTit& pytit = pytits[index(i,i+1)];
          if (pytit != g.terms_pytrees.end())
              add_pycache(pytit->data, inactives[index(i,i+1)]);
          inside_unaryclose(inactives[index(i,i+1)], actives[index(i,i+1)], NULL, g.unarychild_parent_weight, g.rhs_parent_weight);
          
          if (debug >= 20000)
              std::cerr << "# cky::inside() inactives[" << i << "," << i+1 << "] = " 
                  << inactives[index(i,i+1)] << std::endl;
          if (debug >= 20100)
              std::cerr << "# cky::inside() actives[" << i << "," << i+1 << "] = " 
                  << actives[index(i,i+1)] << std::endl;
          if (debug >= 20100) {
              std::cerr << "# cky::inside() pytits[" << i << "," << i+1 << "] = ";
              if (pytits[index(i, i+1)] == g.terms_pytrees.end())
                  std::cerr << "()" << std::endl;
              else
                  std::cerr << pytits[index(i, i+1)]->data << std::endl;
          }
      }

    for (U gap = 2; gap <= n; ++gap) // non-terminals
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (U left = 0; left <= n-gap; ++left) {
            U right = left + gap;
            sS* predictedparents = g.predictive_parse_filter ? &predicteds[index(left,right)] : NULL;
            const StsTit& pytit0 = pytits[index(left, right-1)];
            StsTit& pytit = pytits[index(left, right)];
            if (pytit0 == g.terms_pytrees.end())
                pytit = g.terms_pytrees.end();
            else
                pytit = pytit0->find1(terminals[right-1]);
            S_F& parentinactives = inactives[index(left,right)];
            Stit_F& parentactives = actives[index(left,right)];
            for (U mid = left+1; mid < right; ++mid) {
                const S_F& rightinactives = inactives[index(mid,right)];
                if (rightinactives.empty())
                    continue;
                Stit_F& leftactives = actives[index(left,mid)];
                cforeach (Stit_F, itleft, leftactives) {
                    const Stit leftactive = itleft->first;
                    const F leftprob = itleft->second;
                    cforeach (S_F, itright, rightinactives) {
                        S rightinactive = itright->first;
                        const F rightprob = itright->second;
                        const Stit parentactive = leftactive->find1(rightinactive);
                        if (parentactive != leftactive->end()) {
                            F leftrightprob = leftprob * rightprob;
                            cforeach (S_F, itparent, parentactive->data) {
                                S parent = itparent->first;
                                if (g.predictive_parse_filter && !predictedparents->count(parent))
                                    continue;
                                parentinactives[parent] += leftrightprob * power(itparent->second/afind(g.parent_weight, parent), anneal);
                            }
                            if (!parentactive->key_trie.empty())
                                parentactives[parentactive] += leftrightprob;
                        }
                    }
                }
            }
            // PY correction
            foreach (S_F, it, parentinactives) {
                F pya = g.get_pya(it->first);    // PY cache statistics
                if (pya == 1.0)
                    continue;
                F pyb = g.get_pyb(it->first);
                U pym = dfind(g.parent_pym, it->first);
                U pyn = dfind(g.parent_pyn, it->first);
                it->second *= power( (pym*pya + pyb)/(pyn + pyb), anneal);
            }
            if (pytit != g.terms_pytrees.end())
                add_pycache(pytit->data, parentinactives);
            inside_unaryclose(parentinactives, parentactives, predictedparents, g.unarychild_parent_weight, g.rhs_parent_weight);
            if (debug >= 20000)
                std::cerr << "# cky::inside() inactives[" << left << "," << right 
                    << "] = " << parentinactives << std::endl;
            if (debug >= 20100)
                std::cerr << "# cky::inside() actives[" << left << "," << right << "] = " 
                    << parentactives << std::endl;
            if (debug >= 20100) {
                std::cerr << "# cky::inside() pytits[" << left << "," << right << "] = ";
                if (pytits[index(left, right)] == g.terms_pytrees.end())
                    std::cerr << "()" << std::endl;
                else
                    std::cerr << pytits[index(left, right)]->data << std::endl;
            }
        }
    return dfind(inactives[index(0,n)], start);
  }  // pycky::inside()


  //! random_tree() returns a random parse tree for terminals
  tree* random_tree(S s) {
      U n = terminals.size();
      return random_inactive(s, afind(inactives[index(0, n)], s), 0, n);
  } 
  // pycky::random_tree

  tree* random_tree() { return random_tree(g.start); }

  //! random_inactive() returns a random expansion for an inactive edge
  //
  tree* random_inactive(const S parent, F parentprob, 
			const U left, const U right) const {
      if (left + 1 == right && parent == terminals[left])
          return new tree(parent);
      
      F probthreshold = parentprob * random1();  
      F probsofar = 0;
      F pya = g.get_pya(parent);
      F rulefactor = 1;
      
      if (pya != 1) { // get tree from cache
          F pyb = g.get_pyb(parent);
          U pyn = dfind(g.parent_pyn, parent);
          const StsTit& pytit = pytits[index(left, right)];
          if (pytit != g.terms_pytrees.end())
              cforeach (sT, it, pytit->data) {
                  if ((*it)->cat != parent)
                      continue;
                  probsofar += power( ((*it)->count - pya)/(pyn + pyb), anneal);
                  if (probsofar >= probthreshold)
                      return *it;
              }
          U pym = dfind(g.parent_pym, parent);
          rulefactor = (pym*pya + pyb)/(pyn + pyb);
      }
      
      // tree won't come from cache, so cons up new node
      tree* tp = new tree(parent);
      rulefactor /=  afind(g.parent_weight, parent);
      const S_F& parentinactives = inactives[index(left, right)];
    
      // try unary rules
      cforeach (S_F, it0, parentinactives) {
          S child = it0->first;
          F childprob = it0->second;
          S_S_F::const_iterator it1 = g.unarychild_parent_weight.find(child);
          if (it1 != g.unarychild_parent_weight.end()) {
              const S_F& parent1_weight = it1->second;
              probsofar += childprob * power(dfind(parent1_weight, parent)*rulefactor, anneal);
              if (probsofar >= probthreshold) {
                  tp->children.push_back(random_inactive(child, childprob, left, right));
                  return tp;
              }
          }
      }
      
      // try binary rules
      for (U mid = left+1; mid < right; ++mid) {
          const Stit_F& leftactives = actives[index(left,mid)];
          const S_F& rightinactives = inactives[index(mid,right)];
          cforeach (Stit_F, itleft, leftactives) {
              const Stit leftactive = itleft->first;
              const F leftprob = itleft->second;
              cforeach (S_F, itright, rightinactives) {
                  S rightinactive = itright->first;
                  const F rightprob = itright->second;
                  const Stit parentactive = leftactive->find1(rightinactive);
                  if (parentactive != leftactive->end()) {
                      S_F::const_iterator it = parentactive->data.find(parent);
                      if (it != parentactive->data.end()) {
                          probsofar += leftprob * rightprob * power(it->second*rulefactor, anneal);
                          if (probsofar >= probthreshold) {
                              random_active(leftactive, leftprob, left, mid, tp->children);
                              tp->children.push_back(random_inactive(rightinactive, rightprob, mid, right));
                              return tp;
                          }
                      }
                  }
              }
          }
      }
      
      std::cerr << "\n## Error in pycky::random_inactive(), parent = " << parent
	      << ", left = " << left << ", right = " << right 
	      << ", probsofar = " << probsofar 
	      << " still below probthreshold = " << probthreshold 
	      << std::endl;
      return tp;
  }  // pycky::random_inactive()
  
  void random_active(const Stit parent, F parentprob, \
          const U left, const U right, 
		  tree::ptrs_type& siblings) const ;
}; // pycky{}

struct resample_pycache_helper {
    typedef catcounttree_type tree;
    pycfg_type& g;
    pycky& p;
    resample_pycache_helper(pycfg_type& g, pycky& p) : g(g), p(p) { }

    template <typename Words, typename TreePtrs>
    void operator() (const Words& words, TreePtrs& tps) {
        foreach (typename TreePtrs, tit, tps) {
            tree* tp0 = *tit;
            Ss words;
            tp0 -> terminals(words); // plu_top
            S start = tp0 -> category();  // word

            F old_pya = g.set_pya(start, 1.0);
            F pi0 = g.decrtree(tp0);

            if (pi0 < 0)
                std::cerr << "## pi0 = " << pi0 << ", tp0 = " << tp0 << std::endl;
            assert(pi0 >= 0);

            F r0 = g.tree_prob(tp0);
            assert(r0 >= 0);
            
            F tprob = p.inside(words, start);   // parse string
            if (tprob <= 0)
                std::cerr << "## Error in resample_pycache(): words = " << words << ", tprob = " << tprob
                            << ", tp0 = " << tp0 << std::endl << "## g = " << g << std::endl;
            assert(tprob >= 0);

            tree* tp1 = p.random_tree(start);

            F r1 = g.tree_prob(tp1);
            assert(r1 >= 0);
            
            if (tp0->generalize() == tp1->generalize()) {  // ignore top count
                g.incrtree(tp0);
                tp1->selective_delete();
            }
            else {  // *tp1 != *tp0, do acceptance rejection
                F pi1 = g.incrtree(tp1);
                F pi1r0 = pi1 * r0;
                F pi0r1 = pi0 * r1;
                F accept = (pi0r1 > 0) ? power(pi1r0/pi0r1, p.anneal) : 2.0; // accept if there has been an underflow
                if (random1() <= accept) {
                    // modified by Jackie
                    std::vector<catcounttree_type* > tp0_leaf_nodes, tp1_leaf_nodes;
                    tp0 -> terminal_ptrs(tp0_leaf_nodes); 
                    tp1 -> terminal_ptrs(tp1_leaf_nodes);
                    if (tp0_leaf_nodes.size() != tp1_leaf_nodes.size()) {
                        std::cerr << "tp1 leaf node number != tp0 leaf node number" << std::endl;
                        exit(-1);
                    }
                    else if (tp0_leaf_nodes.size() == 0) {
                        std::cerr << "tp0 leaf node number == 0!" << std::endl;
                        exit(-1);
                    }
                    else {
                        for (size_t i = 0; i < tp0_leaf_nodes.size(); ++i) {
                            std::swap(*tp0_leaf_nodes[i], *tp1_leaf_nodes[i]);
                        }
                    }
                    tp0->generalize().swap(tp1->generalize());  // don't swap top counts
                    tp1->selective_delete();
                }
                else {  // don't accept
                    g.decrtree(tp1);
                    g.incrtree(tp0);
                    tp1->selective_delete();
                }
            }
            g.set_pya(tp0->category(), old_pya);
        }
    }  // resample_pycache_helper::operator()
};  // resample_pycache_helper{}

//! resample_pycache() resamples the strings associated with each cache
inline void resample_pycache(pycfg_type& g, pycky& p) {
    resample_pycache_helper h(g, p);
    p.g.terms_pytrees.for_each(h);
}  // resample_pycache()


#endif // PY_CKY_H
