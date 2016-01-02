#ifndef TOOLKIT_H
#define TOOLKIT_H

#include <vector>

using namespace std;

namespace ToolKit {
  float FindLogMax(vector<float>&);
  void MaxRemovedLogDist(vector<float>&);
  int NormalizeDist(vector<float>&);
  float SumLogs(vector<float>&);
  float SumLogs(float, float);
};

#endif
