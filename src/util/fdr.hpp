#ifndef __FDR_H_
#define __FDR_H_


#include <vector>
#include <utility> // std::pair
#include <algorithm>


void calc_fdr(const std::vector <double> & p_values, std::vector <double> & fdrs, int test_num) {
  // 生成一个带index的p_values, index用于还原fdr的位置
  int len = p_values.size();
  std::vector <std::pair <double, int> > p_values_index (len);
  for (int i = 0; i < len; ++i) {
    p_values_index[i].first = p_values[i];
    p_values_index[i].second = i;
  }
  // 从大到小排序
  std::sort(p_values_index.begin(), p_values_index.end(),
            [&](const std::pair <double, int> & a, const std::pair <double, int> & b) {
                return (a.first > b.first);
            });
  // 计算fdr
  double prev_fdr = 1.0;
  double fdr = 0.0;
  for (int i = 0; i < len; ++i) {
    fdr = p_values_index[i].first * test_num / (len - i);
    if (fdr < prev_fdr) {
      prev_fdr = fdr;
    }
    fdrs[p_values_index[i].second] = prev_fdr;
  }
}


#endif