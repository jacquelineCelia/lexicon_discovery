#include <vector>
#include <iostream>

#include "sym.h"
#include "xtree.h"
#include "utility.h"

template <typename tree_type>
void test() {

  tree_type* x1 = new tree_type("L1");
  tree_type* x2 = new tree_type("L2_T");
  tree_type* x3 = new tree_type("L3_T");

  tree_type* x4 = new tree_type("L4");
  tree_type* x5 = new tree_type("L5_T");
  tree_type* x6 = new tree_type("L6_T");

  x1->children.push_back(x2);
  x1->children.push_back(x3);

  x4->children.push_back(x5);
  x4->children.push_back(x6);

  vector<symbol> btm;
  btm.push_back("1");
  btm.push_back("2");
  x2 -> insert(1, true, btm);
  x2 -> set_current_plu_btm(1);

  btm.clear();
  btm.push_back("3");
  btm.push_back("4");

  x3 -> insert(1, true, btm);
  x3 -> set_current_plu_btm(1);

  x1 -> write_tree(1);
  cerr << endl;

  vector<tree_type*> leaf_nodes_0, leaf_nodes_1;

  x1 -> terminal_ptrs(leaf_nodes_0);
  x4 -> terminal_ptrs(leaf_nodes_1);
  
  if (leaf_nodes_0.size() != leaf_nodes_1.size()) {
      cerr << "mismatched number of leaf nodes!!" << endl;
      exit(-1);
  }
  if (leaf_nodes_0.size() == 0) {
      cerr << "zero leaf nodes!" << endl;
      exit(-1);
  }
  for (size_t i = 0 ; i < leaf_nodes_0.size(); ++i) {
      swap(*leaf_nodes_0[i], *leaf_nodes_1[i]);
  }

  x1 -> write_tree(1);
  cerr << endl;
  x4 -> write_tree(1);
  cerr << endl;

  x1 -> generalize().swap(x4 -> generalize());

  x1 -> write_tree(1);
  cerr << endl;
  x4 -> write_tree(1);
  cerr << endl;

  delete x1, x2, x3, x4, x5, x6;
}

int main(int argc, char** argv) {
  test<catcounttree_type>();

}  // main()
