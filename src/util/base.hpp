#ifndef __BASE_H_
#define __BASE_H_

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>  // needed to use the exit() function
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>  // exception handling
#include <vector>
#include "strim.hpp"

std::string version =
    "GCEN 0.6.1 by Wen Chen (chenwen@biochen.org, "
    "https://www.biochen.org/gcen)\n";

void display_version() { std::cout << version; }

double string_to_log(const std::string &a_str, size_t *idx = 0) {
  double a_double = std::stod(a_str, idx);
  double b_double = std::log(a_double + 1.0);
  return b_double;
}

double string_to_log2(const std::string &a_str, size_t *idx = 0) {
  double a_double = std::stod(a_str, idx);
  double b_double = std::log2(a_double + 1.0);
  return b_double;
}

double string_to_log10(const std::string &a_str, size_t *idx = 0) {
  double a_double = std::stod(a_str, idx);
  double b_double = std::log10(a_double + 1.0);
  return b_double;
}

void load(std::string &in_file_name, std::vector<std::string> &annotation_vec,
          std::vector<std::string> &GeneNameVector,
          std::vector<std::vector<double>> &GeneDataFrame, bool if_log = false,
          bool if_log2 = false, bool if_log10 = false) {
  double (*str_to_double)(const std::string &, size_t *) = std::stod;
  if (if_log) {
    str_to_double = string_to_log;
  }
  if (if_log2) {
    str_to_double = string_to_log2;
  }
  if (if_log10) {
    str_to_double = string_to_log10;
  }

  std::ifstream in_file(in_file_name, std::ios::in);
  if (!in_file.good()) {
    std::cerr << "Error while opening " << in_file_name << ".\n";
    exit(-1);
  }

  std::string line;
  while (getline(in_file, line)) {
    strim(line);
    if (line[0] == '#') {
      annotation_vec.push_back(line);
    }
    std::vector<double> line_list;
    std::stringstream slineString;
    slineString << line;
    std::string singleString;
    try {
      getline(slineString, singleString, '\t');
      std::string gene_name = singleString;

      while (getline(slineString, singleString, '\t')) {
        double item = str_to_double(singleString, 0);
        line_list.push_back(item);
      }

      GeneDataFrame.push_back(line_list);
      GeneNameVector.push_back(gene_name);
      slineString.clear();
    } catch (std::invalid_argument &) {
      slineString.clear();
    }
  }

  in_file.close();
}

std::vector<std::vector<double>> matrix_transpose(
    const std::vector<std::vector<double>> &matrix) {
  int matrix_row_num = matrix.size();
  int matrix_col_num = matrix[0].size();
  std::vector<std::vector<double>> matrix_t(
      matrix_col_num, std::vector<double>(matrix_row_num));
  for (int i = 0; i < matrix_row_num; ++i) {
    for (int j = 0; j < matrix_col_num; ++j) {
      matrix_t[j][i] = matrix[i][j];
    }
  }
  return matrix_t;
}

double upper_quartile(std::vector<double> vec) {
  std::sort(vec.begin(), vec.end());
  double upqt;
  int len = vec.size();
  if (len % 2) {
    int half_len = len / 2 + 1;
    if (half_len % 2) {
      upqt = vec[(half_len + len) / 2];
    } else {
      upqt = (vec[(half_len + len) / 2] + vec[(half_len + len) / 2 - 1]) / 2;
    }
  } else {
    int half_len = len / 2;
    if (half_len % 2) {
      upqt = vec[(half_len + len) / 2];
    } else {
      upqt = (vec[(half_len + len) / 2] + vec[(half_len + len) / 2 - 1]) / 2;
    }
  }
  return upqt;
}

double median(std::vector<double> vec) {
  std::sort(vec.begin(), vec.end());
  double m;
  int len = vec.size();
  if (len % 2) {
    m = vec[len / 2];
  } else {
    m = (vec[len / 2] + vec[len / 2 - 1]) / 2;
  }
  return m;
}

double geometric_mean(const std::vector<double> &a_vector) {
  double geometric = 1.0;
  int len = a_vector.size();
  for (int i = 0; i < len; ++i) {
    geometric = geometric * a_vector[i];
  }
  if (geometric) {
    geometric = std::pow(geometric, 1.0 / len);
  } else {
    geometric = 0.0;
  }
  return geometric;
}

double mean(std::vector<double> &a_vector) {
  double sum = 0.0;
  for (auto item : a_vector) {
    sum = sum + item;
  }
  return sum / a_vector.size();
}

double deviation(std::vector<double> &a_vector) {
  double mean = 0.0;
  for (auto item : a_vector) {
    mean = mean + item;
  }
  mean = mean / a_vector.size();

  double variance = 0.0;
  for (auto item : a_vector) {
    variance = variance + (item - mean) * (item - mean);
  }
  variance = variance / a_vector.size();

  return std::sqrt(variance);
}

void split_string(const std::string &s, std::vector<std::string> &v,
                  const std::string &c) {
  std::string::size_type last_pos = s.find_first_not_of(c, 0);
  std::string::size_type pos = s.find_first_of(c, last_pos);
  while (std::string::npos != pos || std::string::npos != last_pos) {
    v.push_back(s.substr(last_pos, pos - last_pos));
    last_pos = s.find_first_not_of(c, pos);
    pos = s.find_first_of(c, last_pos);
  }
}

void join_vector(const std::vector<std::string> &v, const std::string c,
                 std::string &s) {
  s.clear();
  for (unsigned int i = 0; i < v.size(); ++i) {
    s = s + v[i];
    if (i != v.size() - 1) {
      s = s + c;
    }
  }
}

std::string double_to_string(double d) {
  char tmp[20];
  std::sprintf(tmp, "%e", d);
  std::string s = tmp;
  return s;
}

std::vector<double> get_rank(std::vector<double> &data) {
  int len = data.size();
  std::vector<int> index(len, 0);
  std::vector<double> rank(len, 0);
  for (int i = 0; i < len; ++i) {
    index[i] = i;
  }
  std::sort(index.begin(), index.end(),
            [&](const int &a, const int &b) { return (data[a] > data[b]); });
  int sumranks = 0;
  int dupcount = 0;
  for (int i = 0; i < len; ++i) {
    sumranks = sumranks + i;
    dupcount = dupcount + 1;
    if ((i == len - 1) || (data[index[i]] != data[index[i + 1]])) {
      double averank = double(sumranks) / double(dupcount) + 1;
      for (int j = i - dupcount + 1; j < i + 1; ++j) {
        rank[index[j]] = averank;
      }
      sumranks = 0;
      dupcount = 0;
    }
  }
  return rank;
}

#endif
