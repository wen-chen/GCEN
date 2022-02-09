#include <getopt.h>
#include <fstream>   // std::ifstream, std::ofstream
#include <iterator>  // std::distance, std::advance
#include <random>  // std::random_device, std::mt19937_64, std::uniform_int_distribution
#include <string>                    // std::string
#include <vector>                    // std::vector
#include "third_party/robin_hood.h"  // robin_hood::unordered_map, robin_hood::unordered_set
#include "util/base.hpp"  // version, display_version(), strim(), split_string()

void network_shuffle_help() {
  std::cout << version;
  std::cout << "network_shuffle usage:\n";
  std::cout << "  network_shuffle -i input.network -o output.network\n";
  std::cout << "options:\n";
  std::cout << "  -i --input <input network file>\n";
  std::cout << "  -o --output <output network file>\n";
  std::cout << "  -s --swap <number> multiples of edges number for double-edge "
               "swaps to perform (default 100)\n";
  std::cout << "  -v --version display GCEN version\n";
  std::cout << "  -h --help print help information\n";
  std::cout << "example:\n";
  std::cout << "  network_shuffle -i ../sample_data/test_1.network -o "
               "../sample_data/random_shuffled.network\n";
}

template <typename Iter>
Iter select_randomly(Iter start, Iter end) {
  static std::random_device rd;
  static std::mt19937_64 gen(rd());
  std::uniform_int_distribution<> dis(0, std::distance(start, end) - 1);
  std::advance(start, dis(gen));
  return start;
}

int main(int argc, char* argv[]) {
  if (argc < 2) {
    network_shuffle_help();
    return 0;
  }

  // get option
  std::string in_file_name;
  std::string out_file_name;
  size_t swap_num = 100;

  const char* const short_opts = "hvi:o:s:";
  const struct option long_opts[] = {
      {"help", 0, NULL, 'h'},  {"version", 0, NULL, 'v'},
      {"input", 1, NULL, 'i'}, {"output", 1, NULL, 'o'},
      {"swap", 1, NULL, 's'},  {NULL, 0, NULL, 0}};
  int opt = getopt_long(argc, argv, short_opts, long_opts, NULL);
  while (opt != -1) {
    switch (opt) {
      case 'i':
        in_file_name = optarg;
        break;
      case 'o':
        out_file_name = optarg;
        break;
      case 's':
        swap_num = std::stoi(optarg);
        break;
      case 'h':
        network_shuffle_help();
        return 0;
      case 'v':
        display_version();
        return 0;
      case '?':
        network_shuffle_help();
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

  if (out_file_name.empty()) {
    std::cerr << "Error: -o/--output is required but not specified!\n";
    exit(-1);
  }

  // read network file
  robin_hood::unordered_map<std::string, robin_hood::unordered_set<std::string>>
      network;
  robin_hood::unordered_set<std::string> node_set;
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

    node_set.insert(node_a);
    node_set.insert(node_b);

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

  std::vector<std::string> node_vec(node_set.begin(), node_set.end());

  // swap
  swap_num *= edge_num;
  size_t swap_count = 0;
  while (swap_count < swap_num) {
    std::string node_a = *select_randomly(node_vec.begin(), node_vec.end());
    std::string node_b = *select_randomly(node_vec.begin(), node_vec.end());

    if (node_a == node_b) {  // same source, skip
      continue;
    }

    std::string node_c =
        *select_randomly(network[node_a].begin(), network[node_a].end());
    std::string node_d =
        *select_randomly(network[node_b].begin(), network[node_b].end());

    if (node_c == node_d) {  // same target, skip
      continue;
    }

    if ((network[node_a].find(node_b) == network[node_a].end()) &&
        (network[node_c].find(node_d) ==
         network[node_c].end())) {  // don't create parallel edges
      network[node_a].insert(node_b);
      network[node_b].insert(node_a);
      network[node_c].insert(node_d);
      network[node_d].insert(node_c);

      network[node_a].erase(node_c);
      network[node_c].erase(node_a);
      network[node_b].erase(node_d);
      network[node_d].erase(node_b);
      swap_count += 1;
    }
  }

  // output
  std::ofstream out_file(out_file_name, std::ios::out);
  if (!out_file.good()) {
    std::cerr << "Error while opening " << out_file_name << ".\n";
    exit(-1);
  }

  for (auto& item : network) {
    std::string node_a = item.first;
    for (auto& node_b : item.second) {
      if (node_a.compare(node_b) < 0) {
        out_file << node_a << '\t' << node_b << '\n';
      }
    }
  }

  out_file.close();

  return 0;
}