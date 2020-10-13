#ifndef __RWR_UTIL_H_
#define __RWR_UTIL_H_


#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <cstdlib> // needed to use the exit() function
#include "../Eigen3.3.7/Dense"
#include "func.hpp"
#include "strim.hpp"

using Dynamic_Matrix = Eigen::Matrix <double, Eigen::Dynamic, Eigen::Dynamic>;
using Dynamic_Vector = Eigen::Matrix <double, Eigen::Dynamic, 1>;


class RWR_result {
 public:
  std::string gene;
  double prob;

  bool operator < (const RWR_result & another) {
    if (prob > another.prob) {
      return true;
    }
    return false;
  }
};


int rwr(Dynamic_Matrix & W, Dynamic_Vector & P0, double gamma, std::vector <std::string> & gene_vec, Dynamic_Vector & PT) {
  // network normalization
  for (int i = 0; i < W.cols(); ++i) {
    double col_sum = W.col(i).sum();
    if (col_sum > 0) {
      W.col(i) = W.col(i) / col_sum;
    } else {
      std::string gene = gene_vec[i];
      std::cerr << "Error: " << gene << " is not linked to any other gene.\n";
      exit(-1);
    }
  }

  // Random Walk with Restart
  PT = P0;
  int step = 0;
  double delta = 1.0;
  while (delta > 1e-10) {
    Dynamic_Vector PT1 = (1 - gamma) * W * PT + gamma * P0;
    delta = (PT1 - PT).array().abs().sum();    
    PT = PT1;
    step = step + 1;
  }

  return step;
}

void load_network(std::string & network_file_name, 
    std::vector <std::string> & gene_vec, 
    std::unordered_map <std::string, std::unordered_set <std::string>> & network) {
  // open network file
  std::ifstream network_file(network_file_name, std::ios::in);
  if (!network_file.good()) {
    std::cerr << "Error while opening " << network_file_name << ".\n";
  }

  // read network file
  std::string lineString;
  std::unordered_set <std::string> gene_set;

  while (getline(network_file, lineString)) {
    strim(lineString);
    std::vector <std::string> item;
    split_string(lineString, item, "\t");
    std::string geneA = item[0];
    std::string geneB = item[1];

    gene_set.insert(geneA);
    gene_set.insert(geneB);
    
    if (network.find(geneA) != network.end()) {
      network[geneA].insert(geneB);
    } else {
      network[geneA] = std::unordered_set < std::string > {};
      network[geneA].insert(geneB);
    }

    if (network.find(geneB) != network.end()) {
      network[geneB].insert(geneA);
    } else {
      network[geneB] = std::unordered_set <std::string> {};
      network[geneB].insert(geneA);
    }    
  }

  gene_vec.assign(gene_set.begin(), gene_set.end());
}

void network_2_matrix(std::vector <std::string> & gene_vec,
    std::unordered_map <std::string, std::unordered_set <std::string>> & network, 
    Dynamic_Matrix & matrix) {
  for (int i = 0; i < gene_vec.size(); ++i) {
    matrix(i, i) = 0.0;
    for (int j = i + 1; j < gene_vec.size(); ++j) {
      std::string geneA = gene_vec[i];
      std::string geneB = gene_vec[j];
      if (network[geneA].find(geneB) != network[geneA].end()) {
        matrix(i, j) = 1.0;
        matrix(j, i) = 1.0;
      } else {
        matrix(i, j) = 0.0;
        matrix(j, i) = 0.0;
      }
    }
  }
}

void load_gene(std::string & gene_file_name, 
    std::vector <std::string> & gene_vec,
    Dynamic_Vector & dynamic_vec) {
  std::ifstream gene_file(gene_file_name, std::ios::in);
  if (!gene_file.good()) {
    std::cerr << "Error while opening " << gene_file_name << ".\n";
  }

  std::string gene;
  std::unordered_set <std::string> gene_set;
  while (getline(gene_file, gene)) {
    strim(gene);
    gene_set.insert(gene);
  }

  for (int i = 0; i < gene_vec.size(); ++i) {
    gene = gene_vec[i];
    if (gene_set.find(gene) != gene_set.end()) {
      dynamic_vec[i] = 1.0;
    } else {
      dynamic_vec[i] = 0.0;
    }
  }
}


#endif
