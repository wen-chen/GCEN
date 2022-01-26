#ifndef __NORM_H_
#define __NORM_H_

#include <algorithm>
#include <cmath>
#include <limits>
#include "base.hpp"

double calc_weighted_factor(const std::vector<double>& obs,
                            const std::vector<double>& ref,
                            double logratio_trim = 0.3, double sum_trim = 0.05,
                            double cutoff = -1e+10) {
  double obs_sum = 0.0;
  for (auto& n : obs) obs_sum += n;

  double ref_sum = 0.0;
  for (auto& n : ref) ref_sum += n;

  std::vector<double> log_r;
  std::vector<double> abs_e;
  std::vector<double> v;
  double n = 0.0;
  for (unsigned int i = 0; i < obs.size(); ++i) {
    if ((obs[i] - 0.0 > 1.0e-30) && (ref[i] - 0.0 > 1.0e-30)) {
      if ((log2(obs[i] / obs_sum) + log2(ref[i] / ref_sum)) / 2 > cutoff) {
        log_r.push_back(log2((obs[i] / obs_sum) / (ref[i] / ref_sum)));
        abs_e.push_back((log2(obs[i] / obs_sum) + log2(ref[i] / ref_sum)) / 2);
        v.push_back((obs_sum - obs[i]) / obs_sum / obs[i] +
                    (ref_sum - ref[i]) / ref_sum / ref[i]);
        n += 1.0;
      }
    }
  }

  double low_l = floor(n * logratio_trim) + 1.0;
  double high_l = n + 1.0 - low_l;
  double low_s = floor(n * sum_trim) + 1.0;
  double high_s = n + 1.0 - low_s;

  std::vector<double> log_r_rank = get_rank(log_r);
  std::vector<double> abs_e_rank = get_rank(abs_e);
  double log_r_sum = 0.0;
  double v_sum = 0.0;
  for (unsigned int i = 0; i < log_r.size(); ++i) {
    if ((log_r_rank[i] >= low_l) && (log_r_rank[i] <= high_l) &&
        (abs_e_rank[i] >= low_s) && (abs_e_rank[i] <= high_s)) {
      log_r_sum += log_r[i] / v[i];
      v_sum += 1 / v[i];
    }
  }

  double f = pow(2.0, log_r_sum / v_sum);
  return f;
}

void norm_tmm(std::vector<std::vector<double>>& matrix) {
  std::vector<std::vector<double>> matrix_t = matrix_transpose(matrix);

  // select reference column
  std::vector<double> quantile_vec;
  for (unsigned i = 0; i < matrix_t.size(); ++i) {
    quantile_vec.push_back(upper_quartile(matrix_t[i]));
  }
  double quantile_mean = mean(quantile_vec);
  int ref_col = 0;
  double diff = std::numeric_limits<double>::max();
  for (unsigned i = 0; i < quantile_vec.size(); ++i) {
    if (fabs(quantile_vec[i] - quantile_mean) < diff) {
      ref_col = i;
      diff = fabs(quantile_vec[i] - quantile_mean);
    }
  }

  // std::cout << ref_col << std::endl;

  // calc norm factors
  std::vector<double> factor_vec;
  for (unsigned i = 0; i < matrix_t.size(); ++i) {
    factor_vec.push_back(calc_weighted_factor(matrix_t[i], matrix_t[ref_col]));
  }

  std::vector<double> log_factor_vec;
  for (unsigned i = 0; i < matrix_t.size(); ++i) {
    log_factor_vec.push_back(log(factor_vec[i]));
  }

  double _d = exp(mean(log_factor_vec));
  for (unsigned i = 0; i < matrix_t.size(); ++i) {
    factor_vec[i] = factor_vec[i] / _d;
    // std::cout << factor_vec[i] << std::endl;
  }

  // norm
  for (unsigned i = 0; i < matrix.size(); ++i) {
    for (unsigned j = 0; j < matrix[i].size(); ++j) {
      matrix[i][j] = matrix[i][j] / factor_vec[j];
    }
  }
}

void norm_upqt(std::vector<std::vector<double>>& matrix) {
  std::vector<std::vector<double>> matrix_t = matrix_transpose(matrix);

  std::vector<double> upqt_vec;
  for (unsigned int i = 0; i < matrix_t.size(); ++i) {
    std::vector<double> non_zeros;
    for (unsigned int j = 0; j < matrix_t[i].size(); ++j) {
      if (matrix_t[i][j] > 0) {
        non_zeros.push_back(matrix_t[i][j]);
      }
    }
    upqt_vec.push_back(upper_quartile(non_zeros));
  }
  double upqt_mean = mean(upqt_vec);

  // calc norm factors
  std::vector<double> factor_vec;
  for (unsigned int i = 0; i < upqt_vec.size(); ++i) {
    factor_vec.push_back(upqt_vec[i] / upqt_mean);
  }

  // norm
  for (unsigned i = 0; i < matrix.size(); ++i) {
    for (unsigned j = 0; j < matrix[i].size(); ++j) {
      matrix[i][j] = matrix[i][j] / factor_vec[j];
    }
  }
}

void norm_median(std::vector<std::vector<double>>& matrix) {
  std::vector<std::vector<double>> matrix_t = matrix_transpose(matrix);

  std::vector<double> median_vec;
  for (unsigned int i = 0; i < matrix_t.size(); ++i) {
    std::vector<double> non_zeros;
    for (unsigned int j = 0; j < matrix_t[i].size(); ++j) {
      if (matrix_t[i][j] > 0) {
        non_zeros.push_back(matrix_t[i][j]);
      }
    }
    median_vec.push_back(median(non_zeros));
  }
  double median_mean = mean(median_vec);

  // calc norm factors
  std::vector<double> factor_vec;
  for (unsigned int i = 0; i < median_vec.size(); ++i) {
    factor_vec.push_back(median_vec[i] / median_mean);
  }

  // norm
  for (unsigned i = 0; i < matrix.size(); ++i) {
    for (unsigned j = 0; j < matrix[i].size(); ++j) {
      matrix[i][j] = matrix[i][j] / factor_vec[j];
    }
  }
}

void norm_deseq(std::vector<std::vector<double>>& matrix) {
  //对每一行的数据求几何平均数，如果该行的几何平均数不为零（即该行所有值都不为零），则该行每个数据除以几何平均数，得到一个比例数的矩阵。
  //比例数矩阵和原数据矩阵，列数一样。
  std::vector<std::vector<double>> matrix_ratio;
  for (unsigned int i = 0; i < matrix.size(); ++i) {
    double geo_mean = geometric_mean(matrix[i]);
    if (geo_mean > 0.0) {
      std::vector<double> item;
      for (unsigned int j = 0; j < matrix[i].size(); ++j) {
        item.push_back((matrix[i][j] / geo_mean));
      }
      matrix_ratio.push_back(item);
    }
  }

  // 对比例数矩阵的每一列求中位数。
  std::vector<std::vector<double>> matrix_ratio_t =
      matrix_transpose(matrix_ratio);
  std::vector<double> median_vec;
  for (unsigned int i = 0; i < matrix_ratio_t.size(); ++i) {
    median_vec.push_back(median(matrix_ratio_t[i]));
  }
  double median_mean = mean(median_vec);

  // 原数据矩阵的每一个值除以该列的(几何平均数的中位数)
  std::vector<double> factor_vec;
  for (unsigned int i = 0; i < median_vec.size(); ++i) {
    factor_vec.push_back(median_vec[i] / median_mean);
  }

  // norm
  for (unsigned i = 0; i < matrix.size(); ++i) {
    for (unsigned j = 0; j < matrix[i].size(); ++j) {
      matrix[i][j] = matrix[i][j] / factor_vec[j];
    }
  }
}

void norm_hkg(std::vector<std::vector<double>>& matrix,
              const std::vector<bool>& hkg_vec) {
  // calc median of housekeeping gene expression
  std::vector<double> hkg_val_median_vec;
  unsigned int col_n = matrix[0].size();
  for (unsigned int j = 0; j < col_n; ++j) {
    std::vector<double> hkg_val_vec;
    for (unsigned int i = 0; i < hkg_vec.size(); ++i) {
      if (hkg_vec[i]) {
        double hkg_val = matrix[i][j];
        hkg_val_vec.push_back(hkg_val);
      }
    }
    double hkg_val_median = median(hkg_val_vec);
    hkg_val_median_vec.push_back(hkg_val_median);
  }

  // calc norm factors
  double mean_of_hkg_val_median = mean(hkg_val_median_vec);

  std::vector<double> factor_vec;
  for (unsigned int i = 0; i < hkg_val_median_vec.size(); ++i) {
    double factor = mean_of_hkg_val_median / hkg_val_median_vec[i];
    factor_vec.push_back(factor);
  }

  // norm
  for (unsigned i = 0; i < matrix.size(); ++i) {
    for (unsigned j = 0; j < matrix[i].size(); ++j) {
      matrix[i][j] = matrix[i][j] * factor_vec[j];
    }
  }
}

#endif
