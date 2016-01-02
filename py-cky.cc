#include "py-cky.h"

// added by Jackie
//! get_top_terminals() return plu_tops for plu_bottoms 
void pycky::get_top_terminals(const S& plu_btm, sS& top_terminals, S_S_F& unary_rules) {
    S_S_F::const_iterator it = unary_rules.find(plu_btm);
    if (it != unary_rules.end()) {
        const S_F& rules = it -> second; 
        cforeach (S_F, it1, rules) {
            S parent = it1 -> first;
            if (parent.string_reference().find("T_") != std::string::npos) {
                top_terminals.insert(parent);
            }
        }
    }
}

/*
// added by Jackie
//! get_top_terminals() return plu_tops for plu_bottoms 
void pycky::get_top_terminals(const S& plu_btm, Ss& top_terminals, float min_prob) {
    S_S_F::const_iterator it = g.unarychild_parent_weight.find(plu_btm);
    if (it != g.unarychild_parent_weight.end()) {
        const S_F& rules = it -> second; 
        cforeach (S_F, it1, rules) {
            S parent = it1 -> first;
            if (parent.string_reference().find("T_") != std::string::npos) {
                float rule_prob = it1 -> second / afind(g.parent_weight, parent);
                if (rule_prob >= min_prob) {
                    top_terminals.push_back(parent);
                }
                else {
                    std::cerr << "eliminating " << parent << " --> " << plu_btm << " = " << rule_prob << std::endl;
                }
            }
        }
    }
}
*/

// added by Jackie
//! get_top_terminals() return plu_tops for plu_bottoms 
void pycky::get_top_terminals(const S_F& inactives, sS& top_terminals) {
    S_F::const_iterator iter = inactives.begin();
    for (; iter != inactives.end(); ++iter) {
        if ((iter -> first).string_reference().find("T_") != std::string::npos) {
            top_terminals.insert(iter -> first);
        }
    }
}

// added by Jackie
// ! get_expansion_prob() returns prob(plu_top -> plu_btm) for a terminal
void pycky::get_expansion_prob(const S& terminal, S_F& prob, S_S_F& unary_rules) {
    S_S_F::const_iterator it = unary_rules.find(terminal);
    if (it != unary_rules.end()) {
        const S_F& parent_weight = it-> second;
        cforeach (S_F, it1, parent_weight) {
            S plu_top = it1 -> first;
            prob[plu_top] = power(it1 -> second/afind(g.parent_weight, plu_top), anneal);
        }
    }
}

// added by Jackie
//! get_trees() return trees that expand each top_terminal 
void pycky::get_trees(const sS& top_terminals, StsTit_F& trees, S_F& expansion_probs, St_sT& terms) {
    sS::const_iterator iter = top_terminals.begin();
    for (; iter != top_terminals.end(); ++iter) {
        if (terms.find1(*iter) != terms.end()) {
            if (expansion_probs.find(*iter) != expansion_probs.end()) {
                trees[terms.find1(*iter)] += expansion_probs[*iter];
            }
            else {
                std::cerr << "Can't find prob(plu_top -> plu_btm) in get_trees" << std::endl;
                exit(-1);
            }
        }
    }
}

void pycky::random_plu(Sss& plu_btm, StsTit current_trie_node, \
                       const F parentprob, U left, U right, St_sT& terms) {
    F probthreshold = parentprob * random1();  
    F probsofar = 0;
    bool flag = false;
    for (U mid = right - 1; mid >= right - 2 && mid > left; --mid) {
        StsTit_F& left_trees = tree_leaf_probs[index(left, mid)];
        StsTit_F::const_iterator iter = left_trees.begin(); 
        for (; iter != left_trees.end(); ++iter) {
            sS::const_iterator s_iter = \
                            top_terminals[index(mid, right)].begin();
            for (; s_iter != top_terminals[index(mid, right)].end(); ++ s_iter) {
                StsTit ptr = (iter -> first) -> find1(*s_iter);
                if (ptr == current_trie_node) {
                    flag = true;
                    probsofar += (iter -> second) * inactives[index(mid, right)][*s_iter];
                    if (probsofar >= probthreshold) {
                        Ss tmp(terminals.begin() + mid, terminals.begin() + right);
                        plu_btm.insert(plu_btm.begin(), tmp);
                        random_plu(plu_btm, iter -> first, iter -> second, left, mid, terms);
                        return;
                    } 
                }
            }
        }
    }
    if (flag == false) {
        sS::const_iterator s_iter = top_terminals[index(left, right)].begin();
        for (; s_iter != top_terminals[index(left, right)].end(); ++s_iter) {
            StsTit ptr = terms.find1(*s_iter);
            if (ptr == current_trie_node) {
                flag = true;
                Ss tmp(terminals.begin() + left, terminals.begin() + right);
                plu_btm.insert(plu_btm.begin(), tmp);
                probsofar = tree_leaf_probs[index(left, right)][ptr];
                return;
            }
        }
        if (flag == false) {
            std::cerr << "random_plu is wrong." << std::endl;
            exit(-1);
        }
    } 
}

catcounttree_type* pycky::extend_random_inactive(const S parent, \
        F parentprob, const U left, const U right, Sss& new_terminals, \
        S_S_F& unary_rules, St_S_F& rhs_rules, St_sT& terms) {

    // get to here only when creating a new tree.
    if (parent.string_reference().find("T_") != std::string::npos) {
        catcounttree_type* tp = new catcounttree_type(parent);
        Ss plu_btm(terminals.begin() + left, terminals.begin() + right);
        new_terminals.push_back(plu_btm);
        return tp; 
    }
    // no tree leaves for plu_btm are created
    
    F probthreshold = parentprob * random1();  
    F probsofar = 0;
    F pya = g.get_pya(parent);
    F rulefactor = 1;

    if (pya != 1) { // get tree from cache. 
        F pyb = g.get_pyb(parent);
        U pyn = dfind(g.parent_pyn, parent);
        const StsTit_F& tree_leaf_prob = tree_leaf_probs[index(left, right)];
        StsTit_F::const_iterator iter = tree_leaf_prob.begin();
        for (; iter != tree_leaf_prob.end(); ++iter) {
            cforeach (sT, it, (iter -> first) -> data) {
                if ((*it) -> cat != parent)
                    continue;
                probsofar += power( ((*it) -> count - pya) / (pyn + pyb), anneal) * (iter -> second);
                if (probsofar >= probthreshold) {
                    Sss plu_btm;
                    random_plu(plu_btm, iter -> first, iter -> second, left, right, terms);
                    new_terminals.insert(new_terminals.end(), plu_btm.begin(), plu_btm.end());
                    Ss plu_top;
                    (*it) -> terminals(plu_top);
                    // std::cerr << "assigning " << plu_btm << " to " << plu_top << std::endl;
                    return *it;
                }
            }
        }
        U pym = dfind(g.parent_pym, parent);
        rulefactor = (pym * pya + pyb) / (pyn + pyb);
    }

    // tree won't come from cache, so cons up new node

    tree* tp = new tree(parent);
    rulefactor /=  afind(g.parent_weight, parent);
    const S_F& parentinactives = inactives[index(left, right)];
   
    // try unary rules
    cforeach (S_F, it0, parentinactives) {
        S child = it0 -> first;
        F childprob = it0 -> second;
        S_S_F::const_iterator it1 = unary_rules.find(child);
        if (it1 != unary_rules.end()) {
            const S_F& parent1_weight = it1 -> second;
            probsofar += childprob * power(dfind(parent1_weight, parent) * rulefactor, anneal);
            if (probsofar >= probthreshold) {
                tp -> children.push_back(extend_random_inactive(child, childprob, \
                            left, right, new_terminals, unary_rules, rhs_rules, terms));
                return tp;
            }
        }
    }
    // try binary rules

    for (U mid = left + 1; mid < right; ++mid) {
        const Stit_F& leftactives = actives[index(left, mid)];
        const S_F& rightinactives = inactives[index(mid, right)];
        cforeach (Stit_F, itleft, leftactives) {
            const Stit leftactive = itleft->first;
            const F leftprob = itleft -> second;
            cforeach (S_F, itright, rightinactives) {
                S rightinactive = itright -> first;
                const F rightprob = itright -> second;
                const Stit parentactive = leftactive -> find1(rightinactive);
                if (parentactive != leftactive -> end()) {
                    S_F::const_iterator it = parentactive -> data.find(parent);
                    if (it != parentactive -> data.end()) {
                        probsofar += leftprob * rightprob * power(it -> second * rulefactor, anneal);
                        if (probsofar >= probthreshold) {
                            extend_random_active(leftactive, leftprob, left, mid, tp->children, new_terminals, unary_rules, rhs_rules, terms);
                            tp -> children.push_back(extend_random_inactive(rightinactive, \
                                        rightprob, mid, right, new_terminals, unary_rules, rhs_rules, terms));
                            return tp;
                        }
                    }
                }
            }
        }
    }
    std::cerr << "\n## Error in pycky::extend_random_inactive(), parent = " << parent
	      << ", left = " << left << ", right = " << right 
	      << ", probsofar = " << probsofar 
	      << " still below probthreshold = " << probthreshold 
	      << std::endl;
    exit(-1);
    return tp;
}  // pycky::extend_random_inactive()

void pycky::extend_random_active(const Stit parent, F parentprob, \
        const U left, const U right, tree::ptrs_type& siblings, \
        Sss& new_terminals, S_S_F& unary_rules, St_S_F& rhs_rules, St_sT& terms) {
    
    F probthreshold = random1() * parentprob;
    F probsofar = 0;

    // unary rule
    const S_F& parentinactives = inactives[index(left, right)];
    cforeach (S_F, it, parentinactives) {
        if (rhs_rules.find1(it -> first) == parent) {
            probsofar += it -> second;
            if (probsofar >= probthreshold) {
                siblings.push_back(extend_random_inactive(it -> first, it -> second, left, right, new_terminals, unary_rules, rhs_rules, terms));
                return;
            }
            break;  // only one unary child can possibly generate this parent
        }
    }

    // binary rules
    for (U mid = left + 1; mid < right; ++mid) {
        const Stit_F& leftactives = actives[index(left, mid)];
        const S_F& rightinactives = inactives[index(mid, right)];
        cforeach (Stit_F, itleft, leftactives) {
            const Stit leftactive = itleft -> first;
            const F leftprob = itleft -> second;
            cforeach (S_F, itright, rightinactives) {
                S rightinactive = itright -> first;
                const F rightprob = itright->second;
                if (parent == leftactive -> find1(rightinactive)) {
                    probsofar += leftprob * rightprob;
                    if (probsofar >= probthreshold) {
                        extend_random_active(leftactive, leftprob, left, mid, siblings, new_terminals, unary_rules, rhs_rules, terms);
                        siblings.push_back(extend_random_inactive(rightinactive, rightprob, mid, right, new_terminals, unary_rules, rhs_rules, terms));
                        return;
                    }
                }
            }
        }
    }

    std::cerr << "## Error in pycky::extend_random_active(), parent = " << parent
	      << ", left = " << left << ", right = " << right 
	      << ", probsofar = " << probsofar << ", probthreshold = " << probthreshold 
	      << std::endl;
    exit(-1);
}

void pycky::random_active(const Stit parent, F parentprob, const U left, const U right, 
		     tree::ptrs_type& siblings) const {
      F probthreshold = random1() * parentprob;
      F probsofar = 0;
      // unary rule
      
      const S_F& parentinactives = inactives[index(left, right)];
      cforeach (S_F, it, parentinactives)
          if (g.rhs_parent_weight.find1(it->first) == parent) {
              probsofar += it->second;
              if (probsofar >= probthreshold) {
                  siblings.push_back(random_inactive(it->first, it->second, left, right));
                  return;
              }
              break;  // only one unary child can possibly generate this parent
          }

      // binary rules
      for (U mid = left + 1; mid < right; ++mid) {
          const Stit_F& leftactives = actives[index(left, mid)];
          const S_F& rightinactives = inactives[index(mid, right)];
          cforeach (Stit_F, itleft, leftactives) {
              const Stit leftactive = itleft->first;
              const F leftprob = itleft->second;
              cforeach (S_F, itright, rightinactives) {
                  S rightinactive = itright->first;
                  const F rightprob = itright->second;
                  if (parent == leftactive->find1(rightinactive)) {
                      probsofar += leftprob * rightprob;
                      if (probsofar >= probthreshold) {
                          random_active(leftactive, leftprob, left, mid, siblings);
                          siblings.push_back(random_inactive(rightinactive, rightprob, mid, right));
                          return;
                      }
                  }
              }
          }
      }
      
      std::cerr << "## Error in pycky::random_active(), parent = " << parent
	      << ", left = " << left << ", right = " << right 
	      << ", probsofar = " << probsofar << ", probthreshold = " << probthreshold 
	      << std::endl;
      return;
}  // pycky::random_active()


