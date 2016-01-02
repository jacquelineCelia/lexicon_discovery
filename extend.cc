  template <typename terminals_type>
  F extend_inside(const terminals_type& terminals0, S start) {
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
          get_expansion_prob(terminals[i], terminal_expansion_prob);

          // get top_terminals
          get_top_terminals(terminals[i], top_terminals[index(i, i + 1)]);

          // get trees that expand to top_terminals
          StsTit_F& tree_leaf_prob = tree_leaf_probs[index(i, i + 1)];
          get_trees(top_terminals[index(i, i + 1)], tree_leaf_prob, terminal_expansion_prob);
           
          StsTit_F::const_iterator tree_iter = tree_leaf_prob.begin(); 
          for (; tree_iter != tree_leaf_prob.end(); ++tree_iter) {
              if ((tree_iter -> first) != g.filtered_terms_pytrees.end()) {
                  add_pycache((tree_iter -> first) -> data, \
                          inactives[index(i, i + 1)], tree_iter -> second);
              }  
          }
        
          inside_unaryclose(inactives[index(i, i + 1)], actives[index(i, i + 1)], \
                  g.predictive_parse_filter ? &predicteds[index(i,i+1)] : NULL);

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
                                      parentinactives[parent] += leftrightprob * power(itparent -> second / afind(g.parent_weight, parent), anneal);
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
                  get_trees(top_terminals[index(left, right)], tree_leaf_prob, parentinactives);
              }

              // need to get all the py trees for this cell
              // get top_terminals that map to one terminal and two terminals
              for (U mid = right - 1 ; mid >= right - 2 && mid > left ; --mid) {
                  const StsTit_F& tables = tree_leaf_probs[index(left, mid)];
                  StsTit_F::const_iterator iter = tables.begin();
                  for (; iter != tables.end(); ++iter) {
                      if ((iter -> first) != g.filtered_terms_pytrees.end()) {
                          Ss::const_iterator Ss_iter = \
                          top_terminals[index(mid, right)].begin();
                          for (; Ss_iter != top_terminals[index(mid, right)].end(); ++Ss_iter) {
                              StsTit ptr = (iter -> first) -> find1(*Ss_iter);
                              if (ptr != g.filtered_terms_pytrees.end()) {
                                  F plu_btm_prob = inactives[index(mid, right)][*Ss_iter]; // already annealed
                                  tree_leaf_prob[ptr] += (iter -> second) * plu_btm_prob;
                              }
                          }
                      }
                  }
              }

              StsTit_F::const_iterator iter = tree_leaf_prob.begin();

              for (; iter != tree_leaf_prob.end(); ++iter) {
                  if ((iter -> first) != g.filtered_terms_pytrees.end()) {
                      add_pycache((iter -> first) -> data, parentinactives, iter -> second);
                  }
              }
              inside_unaryclose(parentinactives, parentactives, NULL);

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
