#ifndef __RWR_H_
#define __RWR_H_

#include <cstdlib>  // needed to use the exit() function
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "../third_party/Eigen3.3.7/Dense"
#include "../third_party/robin_hood.h"
#include "base.hpp"
#include "strim.hpp"

using Dynamic_Matrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
using Dynamic_Vector = Eigen::Matrix<double, Eigen::Dynamic, 1>;

class RWR_result {
 public:
  std::string gene;
  double prob;
};

bool operator<(const RWR_result &a, const RWR_result &b) {
  if (a.prob > b.prob) {
    return true;
  }
  return false;
}

int rwr(Dynamic_Matrix &W, Dynamic_Vector &P0, double gamma,
        Dynamic_Vector &PT) {
  // network normalization
  for (int i = 0; i < W.cols(); ++i) {
    double col_sum = W.col(i).cwiseAbs().sum();
    if (col_sum > 0) {
      W.col(i) = W.col(i) / col_sum;
    }
  }

  // Random Walk with Restart
  PT = P0;
  int step = 0;
  double delta = 1.0;
  while ((delta > 1e-10) && (step < 1000000)) {
    Dynamic_Vector PT1 = (1 - gamma) * W * PT + gamma * P0;
    delta = (PT1 - PT).array().abs().sum();
    PT = PT1;
    step = step + 1;
  }

  return step;
}

void load_network(
    std::string &network_file_name, bool if_directed_network,
    bool if_weighted_network, std::vector<std::string> &gene_vec,
    robin_hood::unordered_map<
        std::string, robin_hood::unordered_map<std::string, double>> &network) {
  // open network file
  std::ifstream network_file(network_file_name, std::ios::in);
  if (!network_file.good()) {
    std::cerr << "Error while opening " << network_file_name << ".\n";
  }

  // read network file
  std::string line;
  robin_hood::unordered_set<std::string> gene_set;

  while (getline(network_file, line)) {
    strim(line);
    if (line[0] == '#') {
      continue;
    }
    std::vector<std::string> str_vec;
    split_string(line, str_vec, "\t");
    std::string gene_a = str_vec[0];
    std::string gene_b = str_vec[1];

    double weight = 1.0;
    if ((if_weighted_network) && (str_vec.size() > 2)) {
      try {
        weight = std::stod(str_vec[2]);
      } catch (std::invalid_argument &) {
        std::cerr << "Error: the edge weights of network must be numerical.\n";
        exit(-1);
      }
    }

    gene_set.insert(gene_a);
    gene_set.insert(gene_b);

    if (network.find(gene_a) != network.end()) {
      network[gene_a][gene_b] = weight;
    } else {
      network[gene_a] = robin_hood::unordered_map<std::string, double>{};
      network[gene_a][gene_b] = weight;
    }

    if (!if_directed_network) {
      if (network.find(gene_b) != network.end()) {
        network[gene_b][gene_a] = weight;
      } else {
        network[gene_b] = robin_hood::unordered_map<std::string, double>{};
        network[gene_b][gene_a] = weight;
      }
    }
  }

  gene_vec.assign(gene_set.begin(), gene_set.end());
}

void network_2_matrix(
    std::vector<std::string> &gene_vec,
    robin_hood::unordered_map<
        std::string, robin_hood::unordered_map<std::string, double>> &network,
    Dynamic_Matrix &matrix) {
  for (unsigned int i = 0; i < gene_vec.size(); ++i) {
    matrix(i, i) = 0.0;
    for (unsigned int j = i + 1; j < gene_vec.size(); ++j) {
      std::string gene_a = gene_vec[i];
      std::string gene_b = gene_vec[j];

      if ((network.find(gene_a) != network.end()) &&
          (network[gene_a].find(gene_b) != network[gene_a].end())) {
        matrix(i, j) = network[gene_a][gene_b];
      } else {
        matrix(i, j) = 0.0;
      }

      if ((network.find(gene_b) != network.end()) &&
          (network[gene_b].find(gene_a) != network[gene_b].end())) {
        matrix(j, i) = network[gene_b][gene_a];
      } else {
        matrix(j, i) = 0.0;
      }
    }
  }
}

void load_gene(std::string &gene_file_name, std::vector<std::string> &gene_vec,
               bool if_weighted_gene, Dynamic_Vector &dynamic_vec) {
  std::ifstream gene_file(gene_file_name, std::ios::in);
  if (!gene_file.good()) {
    std::cerr << "Error while opening " << gene_file_name << ".\n";
  }

  std::string line;
  robin_hood::unordered_map<std::string, double> gene_map;
  while (getline(gene_file, line)) {
    strim(line);
    if (line[0] == '#') {
      continue;
    }
    std::vector<std::string> str_vec;
    split_string(line, str_vec, "\t");
    std::string gene = str_vec[0];

    double weight = 1.0;
    if ((if_weighted_gene) && (str_vec.size() > 1)) {
      try {
        weight = std::stod(str_vec[1]);
      } catch (std::invalid_argument &) {
        std::cerr << "Error: the weight of seed genes must be numerical.\n";
        exit(-1);
      }
    }

    gene_map[gene] = weight;
  }

  for (unsigned int i = 0; i < gene_vec.size(); ++i) {
    std::string gene = gene_vec[i];
    if (gene_map.find(gene) != gene_map.end()) {
      dynamic_vec[i] = gene_map[gene];
    } else {
      dynamic_vec[i] = 0.0;
    }
  }
}

#endif
