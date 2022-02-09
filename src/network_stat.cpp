#include <getopt.h>
#include <fstream>                   // std::ifstream
#include <iostream>                  // std::cout, std::cerr
#include <string>                    // std::string
#include <vector>                    // std::vector
#include "third_party/robin_hood.h"  // robin_hood::unordered_map, robin_hood::unordered_set
#include "util/base.hpp"  // version, display_version(), strim(), split_string()

void network_stat_help() {
  std::cout << version;
  std::cout << "network_stat usage:\n";
  std::cout << "  network_stat -i input_file\n";
  std::cout << "options:\n";
  std::cout << "  -i --input <input file>\n";
  std::cout << "  -v --version display GCEN version\n";
  std::cout << "  -h --help print help information\n";
  std::cout << "example:\n";
  std::cout << "  network_stat -i ../sample_data/gene_co_expr.network\n";
}

void dfs(const std::string& node,
         robin_hood::unordered_map<
             std::string, robin_hood::unordered_set<std::string>>& network,
         robin_hood::unordered_map<std::string, bool>& node_state,
         robin_hood::unordered_set<std::string>& component) {
  if (node_state[node]) {
    component.insert(network[node].begin(), network[node].end());
    node_state[node] = false;
    for (auto& neighbor : network[node]) {
      dfs(neighbor, network, node_state, component);
    }
  }
}

int main(int argc, char* argv[]) {
  if (argc < 2) {
    network_stat_help();
    return 0;
  }

  // get option
  std::string in_file_name;

  const char* const short_opts = "hvi:";
  const struct option long_opts[] = {{"help", 0, NULL, 'h'},
                                     {"version", 0, NULL, 'v'},
                                     {"input", 1, NULL, 'i'},
                                     {NULL, 0, NULL, 0}};
  int opt = getopt_long(argc, argv, short_opts, long_opts, NULL);
  while (opt != -1) {
    switch (opt) {
      case 'i':
        in_file_name = optarg;
        break;
      case 'h':
        network_stat_help();
        return 0;
      case 'v':
        display_version();
        return 0;
      case '?':
        network_stat_help();
        return 0;
      case -1:
        break;
      default:
        return -1;
    }
    opt = getopt_long(argc, argv, short_opts, long_opts, NULL);
  }

  // check options
  if (in_file_name.empty()) {
    std::cerr << "Error: -i/--input is required but not specified!\n";
    exit(-1);
  }

  // stat
  robin_hood::unordered_map<std::string, robin_hood::unordered_set<std::string>>
      network;
  robin_hood::unordered_map<std::string, bool> node_state;
  size_t edge_num = 0;

  std::ifstream in_file(in_file_name, std::ios::in);
  if (!in_file.good()) {
    std::cerr << "Error while opening " << in_file_name << ".\n";
  }
  std::string line;
  while (getline(in_file, line)) {
    strim(line);
    if (line[0] == '#') {
      continue;
    }
    std::vector<std::string> str_vec;
    split_string(line, str_vec, "\t");
    std::string node_a = str_vec[0];
    std::string node_b = str_vec[1];

    node_state[node_a] = true;
    node_state[node_b] = true;

    if (network.find(node_a) == network.end()) {
      network[node_a] = robin_hood::unordered_set<std::string>{};
    }
    if (network[node_a].find(node_b) == network[node_a].end()) {
      network[node_a].insert(node_b);
      edge_num += 1;
    }

    if (network.find(node_b) == network.end()) {
      network[node_b] = robin_hood::unordered_set<std::string>{};
    }
    if (network[node_b].find(node_a) == network[node_b].end()) {
      network[node_b].insert(node_a);
    }
  }

  in_file.close();

  size_t node_num = node_state.size();
  if (node_num < 3) {
    std::cerr << "Error: the network is too small!\n";
    exit(-1);
  }

  size_t component_num = 0;
  size_t max_component = 0;
  size_t degree_sum = 0;

  for (auto& node_item : node_state) {
    std::string node = node_item.first;
    bool state = node_item.second;

    if (state) {
      component_num += 1;
      robin_hood::unordered_set<std::string> component;
      dfs(node, network, node_state, component);
      size_t component_size = component.size();
      if (component_size > max_component) {
        max_component = component_size;
      }
    }

    int degree = network[node].size();
    degree_sum += degree;
  }

  double density = edge_num * 2.0 / node_num / (node_num - 1);
  double degree_mean = ((double)degree_sum) / node_num;

  double clustering_coefficient_sum = 0.0;
  size_t triangle_sum = 0;
  size_t triangle_ctrl = 0;
  for (auto& node_item : network) {
    size_t triangle = 0;
    robin_hood::unordered_set<std::string> neighbour_set = node_item.second;
    size_t neighbour_num = neighbour_set.size();
    if (neighbour_num < 2) {
      continue;
    }
    std::vector<std::string> neighbour_vec(neighbour_set.begin(),
                                           neighbour_set.end());
    for (size_t i = 0; i < neighbour_num; ++i) {
      for (size_t j = i + 1; j < neighbour_num; ++j) {
        std::string neighbour_a = neighbour_vec[i];
        std::string neighbour_b = neighbour_vec[j];
        if (network[neighbour_a].find(neighbour_b) !=
            network[neighbour_a].end()) {
          triangle += 1;
        }
      }
    }

    triangle_sum += triangle;
    triangle_ctrl += neighbour_num * (neighbour_num - 1);

    double clustering_coefficient =
        triangle * 2.0 / ((neighbour_num - 1) * neighbour_num);
    clustering_coefficient_sum += clustering_coefficient;
  }

  double average_clustering_coefficient = clustering_coefficient_sum / node_num;
  double transitivity = triangle_sum * 2.0 / triangle_ctrl;

  // output
  std::cout << "Number of nodes: " << node_num
            << "; Number of edges: " << edge_num
            << ";\nConnected components: " << component_num
            << "; The giant component: " << max_component
            << ";\nNetwork density: " << density
            << "; Average node degree: " << degree_mean
            << ";\nAverage clustering coefficient: "
            << average_clustering_coefficient
            << "; Transitivity: " << transitivity << ".\n";

  return 0;
}
