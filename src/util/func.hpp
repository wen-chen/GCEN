#ifndef __FUNC_H_
#define __FUNC_H_


#include <stdexcept> // exception handling
#include <cstdlib> // needed to use the exit() function
#include <iostream>
#include <fstream> 
#include <sstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <cstdio>
#include "strim.hpp"


void display_version() {
  std::cout << "GCEN 0.5.0 by Wen Chen (chenwen@biochen.com, httpw://www.biochen.com/gcen)\n";
}


double string_to_log2(const std::string & a_str, size_t * idx = 0) {
  double a_double = std::stod(a_str, idx);
  double b_double = std::log2(a_double + 1.0);
  return b_double;
}


double string_to_log10(const std::string & a_str, size_t * idx = 0) {
  double a_double = std::stod(a_str, idx);
  double b_double = std::log10(a_double + 1.0);
  return b_double;
}


void load(std::string & InFileName, std::vector <std::string> & GeneNameVector,
    std::vector <std::vector <double> > & GeneDataFrame, bool if_log2 = false, bool if_log10 = false) {
  double (* str_to_double) (const std::string &, size_t *) = std::stod;
  if (if_log2) {
    str_to_double = string_to_log2;
  }
  if (if_log10) {
    str_to_double = string_to_log10;
  }
  std::ifstream inFile(InFileName, std::ios::in);
  if (!inFile.good()) {
    std::cerr << "Error while opening " << InFileName << ".\n";
    exit(-1);
  }

  std::string lineString;
  while (getline(inFile, lineString)) {
    strim(lineString);
    std::vector <double> line_list;
    std::stringstream slineString;
    slineString << lineString;
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
    } catch (std::invalid_argument) {
      slineString.clear();
    }
  }
  
  inFile.close();
}


void matrix_transpose(const std::vector <std::vector <double> > & Matrix1,
    std::vector <std::vector <double> > & Matrix2) {
  int Matrix1_RowNum = Matrix1.size();
  int Matrix1_ColNum = Matrix1[0].size();
  int Matrix2_RowNum = Matrix2.size();
  int Matrix2_ColNum = Matrix2[0].size();
  if (Matrix1_RowNum == Matrix2_ColNum && Matrix1_ColNum == Matrix2_RowNum) {
    for (int i = 0; i < Matrix1_RowNum; ++i) {
      for (int j = 0; j < Matrix1_ColNum; ++j) {
        Matrix2[j][i] = Matrix1[i][j];
      }
    }
  } else {
    std::cerr << "Error in function matrix_transpose.\n";
    exit(-1);
  }
}


double sorted_vector_median(const std::vector <double> & sorted_vector) {
  double median;
  int len = sorted_vector.size();
  if (len % 2) {
    median = sorted_vector[len / 2];
  } else {
    median = (sorted_vector[len / 2] + sorted_vector[len / 2 - 1]) / 2;
  }
  return median;
}


double sorted_vector_upperquartile(const std::vector <double> & sorted_vector) {
  double upperquartile;
  int len = sorted_vector.size();
  if (len % 2) {
    int half_len = len / 2 + 1;
    if (half_len % 2) {
      upperquartile = sorted_vector[(half_len + len) / 2];
    } else {
      upperquartile = (sorted_vector[(half_len + len) / 2] + sorted_vector[(half_len + len) / 2 - 1]) / 2;
    }
  } else {
    int half_len = len / 2;
    if (half_len % 2) {
      upperquartile = sorted_vector[(half_len + len) / 2];
    } else {
      upperquartile = (sorted_vector[(half_len + len) / 2] + sorted_vector[(half_len + len) / 2 - 1]) / 2;
    }
  }
  return upperquartile;
}


double geometric_mean(const std::vector <double> & a_vector) {
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


double mean(std::vector <double> & a_vector) {
  double sum = 0.0;
  for (auto item : a_vector) {
    sum = sum + item;
  }
  return sum / a_vector.size();
}


double deviation(std::vector <double> & a_vector) {
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


void split_string(const std::string & s, std::vector <std::string> & v, const std::string & c) {
  std::string::size_type last_pos = s.find_first_not_of(c, 0);
  std::string::size_type pos = s.find_first_of(c, last_pos);
  while (std::string::npos != pos || std::string::npos != last_pos) {
    v.push_back(s.substr(last_pos, pos - last_pos));
    last_pos = s.find_first_not_of(c, pos);
    pos = s.find_first_of(c, last_pos);        
  }
}


void join_vector(const std::vector <std::string>& v, const std::string c, std::string& s) {
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


std::vector <double> GetRanks(std::vector <double> & data) {
  int len = data.size();
  std::vector <int> index(len, 0);
  std::vector <double> rank(len, 0);
  for (int i = 0; i < len; ++i) {
    index[i] = i;
  }
  std::sort(index.begin(), index.end(),
            [&](const int & a, const int & b) {
                return (data[a] > data[b]);
            }
  );
  int sumranks = 0;
  int dupcount = 0;
  for (int i = 0; i < len; ++i) {
    sumranks = sumranks + i;
    dupcount = dupcount + 1;
    if ((i == len -1) || (data[index[i]] != data[index[i + 1]]) ) {
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