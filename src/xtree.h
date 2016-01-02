// xtree.h
//
// Uses the "Barton and Nackman trick" as described in
// http://osl.iu.edu/~tveldhui/papers/techniques/techniques01.html#l14

#ifndef XTREE_H
#define XTREE_H

#include "custom-allocator.h"  // must be first

#include <algorithm>
#include <iostream>
#include <set>
#include <vector>
#include <utility>
#include <map>

#include "sym.h"

extern bool compact_trees;

//! xtree_type{} is a generic tree node type.  It consists of a generic label
//! together with a sequence of children pointers.
//
template <typename special_type> 
struct xtree_type {
  typedef xtree_type<special_type> general_type;

  typedef special_type* ptr_type;
  typedef std::vector<ptr_type> ptrs_type;
  typedef std::vector<symbol> Ss;

  ptrs_type children;
  symbol    cat;

  xtree_type<special_type>(symbol cat=symbol::undefined()) : cat(cat) { }

  //! delete_tree() deletes this node and all of its children
  //
  void delete_tree() {
    for (typename ptrs_type::iterator it = children.begin(); it != children.end(); ++it)
      (*it)->specialize().delete_tree();
    delete this;
  }  // xtree_type::delete_tree()
  
  //! category() returns the category of this node
  //
  symbol category() const { return cat; }

  //! specialize() converts an xtree pointer to the more specific treeptr type
  //
  special_type& specialize() { return static_cast<special_type&>(*this); }
  const special_type& specialize() const { return static_cast<const special_type&>(*this); }

  //! generalize() converts an xtree pointer to the more general treeptr type
  //
  general_type& generalize() { return static_cast<general_type&>(*this); }
  const general_type& generalize() const { return static_cast<general_type&>(*this); }

  //! swap() swaps the contents of two tree nodes
  //
  void swap(general_type& t) {
    std::swap(children, t.children);
    std::swap(cat, t.cat);
  }

  /* added by Jackie */
  //! is_plu_top() returns whether this node is for a plu_top sym
  bool is_plu_top() {
      return cat.string_reference().find("T_") == std::string::npos ? false : true;
  } 

  //! terminals() returns the terminal yield of the tree
  //
  template <typename terminals_type>
  void terminals(terminals_type& terms) const {
    if (children.empty())
      terms.push_back(specialize().category());
    else 
      for (typename ptrs_type::const_iterator it = children.begin(); it != children.end(); ++it)
	(*it)->terminals(terms);
  } // xtree_type::terminals()

  std::ostream& write_tree(std::ostream& os) const {
      if (children.empty()) 
          return specialize().write_label(os);
      else {
          os << '(';
          specialize().write_label(os);
          for (typename ptrs_type::const_iterator it = children.begin();
                  it != children.end(); ++it) 
              (*it)->specialize().write_tree(os << ' ');
          return os << ')';
      }
  }  // xtree_type::write_tree()

  std::ostream& write_label(std::ostream& os) const {
    return os << cat;
  }  // xtree_type::write_label()

  bool operator== (const general_type& t) const {
    return cat == t.cat && specialize().equal_children(t.specialize());
  }  // xtree_type::operator==

  bool equal_children(const general_type& t) const {
    typename ptrs_type::const_iterator it0 = children.begin();
    typename ptrs_type::const_iterator it1 = t.children.begin();
    for ( ; it0 != children.end(); ++it0, ++it1) {
        if (it1 == t.children.end())
            return false;
        if (!(**it0 == **it1))
            return false;
    }
    return it1 == t.children.end();
  }  // xtree_type::equal_children()

};

template<typename special_type>
std::ostream& operator<< (std::ostream& os, const xtree_type<special_type>& xt) {
  return xt.specialize().write_tree(os);
}

template<typename special_type>
std::ostream& operator<< (std::ostream& os, const xtree_type<special_type>* xtp) {
  return xtp->specialize().write_tree(os);
}

//! cattree_type{} labels only contains a category label
//
struct cattree_type : public xtree_type<cattree_type> {
  cattree_type(symbol cat=symbol()) : xtree_type<cattree_type>(cat) { }
};  // cattree_type{}

//! catcounttree_type{} labels contain a category label and an integer count
//
class catcounttree_type : public xtree_type<catcounttree_type> {
 
public:

  typedef unsigned int count_type;

  count_type count;

  /* added by Jackie */
  typedef std::map<unsigned int, std::vector<symbol> > U_Ss;

  catcounttree_type(symbol cat=symbol(), count_type count=0) 
    : xtree_type<catcounttree_type>(cat), count(count) { }

  //! terminal_ptrs() returns the terminal yield of the tree
  //
  template <typename terminals_type>
  void terminal_ptrs(std::vector<terminals_type*>& leaf_nodes) {
    if (children.empty())
        leaf_nodes.push_back(this);
    else 
        for (typename ptrs_type::iterator it = children.begin(); it != children.end(); ++it)
            (*it) -> terminal_ptrs(leaf_nodes);
  } // xtree_type::terminal_ptrs()

  bool operator== (const catcounttree_type& t) const {
    return cat == t.cat && count == t.count && equal_children(t);
  } 

  typedef std::set<ptr_type> sptr_type;

  //! selective_delete() deletes all nodes from the top of the tree that have a count
  //! of zero
  //
  void selective_delete() {
    sptr_type todelete;
    selective_delete_helper(todelete);
    assert((count != 0) == todelete.empty());
    for (sptr_type::iterator it = todelete.begin(); it != todelete.end(); ++it)
      delete *it;
  }  // catcounttree_type::selective_delete()

private:

  void selective_delete_helper(sptr_type& todelete) {
    if (count == 0) {
        todelete.insert(this);
        for (ptrs_type::iterator it = children.begin(); it != children.end(); ++it) {
            (*it)->selective_delete_helper(todelete);
        }
    }
  }  // catcounttree_type::selective_delete_helper()
  
public:

  //! swap() swaps the contents of two catcounttrees
  //
  void swap(catcounttree_type& t) {
    std::swap(count, t.count);
    this->generalize().swap(t.generalize());
  }

  void write_label(std::vector<std::vector<symbol> >& terminals) { 
      (compact_trees || count == 0)  ? (std::cerr << cat) : (std::cerr << cat << '#' << count); 
      if (is_plu_top()) {
          Ss plu_btms = terminals[0]; 
          Ss::const_iterator s_iter = plu_btms.begin();
          for (; s_iter != plu_btms.end(); ++s_iter) {
              std::cerr << " [" << (*s_iter) << "] ";
          }
          terminals.erase(terminals.begin());
      }
  }

  std::ostream& write_label(std::ostream& os, \
      std::vector<std::vector<symbol> >& terminals) { 
      (compact_trees || count == 0)  ? (os << cat) : (os << cat << '#' << count); 
      if (is_plu_top()) {
          Ss plu_btms = terminals[0];
          Ss::iterator s_iter = plu_btms.begin();
          for (; s_iter != plu_btms.end(); ++s_iter) {
              os << " [" << (*s_iter) << "] ";
          }
          terminals.erase(terminals.begin());
      }
      return os; 
  }

  std::ostream& write_label(std::ostream& os) const { 
      return (compact_trees || count == 0)  ? (os << cat) : (os << cat << '#' << count); }

  std::ostream& write_tree(std::ostream& os, \
      std::vector<std::vector<symbol> >& terminals) { 
      if (children.empty()) 
          return specialize().write_label(os, terminals);
      else {
          os << '(';
          specialize().write_label(os);
          for (typename ptrs_type::const_iterator it = children.begin();
                  it != children.end(); ++it) {
              os << ' ';
              (*it)->specialize().write_tree(os, terminals);
          }
          return os << ')';
      }
  }  // xtree_type::write_tree()

  void write_tree(std::vector<std::vector<symbol> >& terminals) { 
      if (children.empty()) 
          return specialize().write_label(terminals);
      else {
          std::cerr << '(';
          specialize().write_label(terminals);
          for (typename ptrs_type::const_iterator it = children.begin();
                  it != children.end(); ++it) {
              std::cerr << ' ';
              (*it)->specialize().write_tree(terminals);
          }
          std::cerr << ')';
      }
  }  // xtree_type::write_tree()
  

  std::ostream& write_tree(std::ostream& os) const {
    if (children.empty()) 
      return write_label(os);
    else if (compact_trees && count == 0) {
      for (ptrs_type::const_iterator it = children.begin();
	   it != children.end(); ++it) {
	if (it != children.begin())
	  os << ' ';
	(*it)->specialize().write_tree(os);
      }
      return os;
    }
    else {
      os << '(';
      specialize().write_label(os);
      for (ptrs_type::const_iterator it = children.begin();
	   it != children.end(); ++it) 
	(*it)->specialize().write_tree(os << ' ');
      return os << ')';
    }
  }  // xtree_type::write_tree()

};  //catcounttree_type{}

#endif // XTREE_H
