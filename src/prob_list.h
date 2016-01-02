#ifndef PROB_LIST_H
#define PROB_LIST_H

#include <vector>
#include <map>
#include "toolkit.h"

#define MIN_PROB_VALUE -70000000

using namespace std;

template <class T>
class ProbList {
 public:
  ProbList() {
      _need_to_update = true;
  }
  ProbList(int size) {
      _need_to_update = true;
      _index.resize(size);
      _probs.resize(size);
  }
  void resize(int size) {
      _index.resize(size);
      _probs.resize(size);
      _need_to_update = true;
  }
  void reserve(int size) {
      _index.reserve(size);
      _probs.reserve(size);
      _need_to_update = true;
  }
  ~ProbList() {};
  void push_back(float prob, T index) {
      _probs.push_back(prob);
      _index.push_back(index);
      _need_to_update = true;
  }
  void clear() {
      _probs.clear();
      _index.clear();
      _need_to_update = true;
  }
  void assign(vector<float> probs) {
      _probs = probs;
      _need_to_update = true;
  }
  float value() {
      if (_probs.size() == 0) {
        return MIN_PROB_VALUE;
      }
      if (_need_to_update) {
        _value = ToolKit::SumLogs(_probs);
        _need_to_update = false;
      }
      return _value;
  }
  vector<float> probs() {
      return _probs;
  }
  vector<T> index() {
      return _index;
  }
  void add_item(float prob, T index, int id) {
    _probs[id] = prob;
    _index[id] = index;
    _need_to_update = true;
  }
  void assign(vector<float> probs, vector<T> index) {
      _probs = probs;
      _index = index;
      _need_to_update = true;
  }
  T index(int i) {
     return _index[i];
  }
 private:
  vector<T> _index;
  vector<float> _probs;
  float _value;
  bool _need_to_update;
};

#endif
