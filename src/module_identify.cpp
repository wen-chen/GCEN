#include <getopt.h>
#include <iostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <fstream>
#include <sstream>
#include <chrono>
#include <thread>
#include <mutex>
#include "util/func.hpp"
#include "util/strim.hpp"


static std::mutex mutex_lock;


void module_identify_help() {
  std::cout << "GCEN 0.5.0 by Wen Chen (chenwen@biochen.com, https://www.biochen.com/gcen)\n";
  std::cout << "module_identify usage:\n";
  std::cout << "  module_identify -i input_file -o output_file\n";
  std::cout << "options:\n";
  std::cout << "  -i --input <input file>\n";
  std::cout << "  -o --output <output file>\n";
  std::cout << "  -s --similarity <number> similarity cutoff (default: 0.5)\n";
  std::cout << "  -t --thread <number> cpu cores (default: 2)\n";
  std::cout << "  -v --version display GCEN version\n";
  std::cout << "  -h --help print help information\n";
  std::cout << "example:\n";
  std::cout << "  module_identify -i ../sample_data/gene_co_expr.network -o ../sample_data/module.txt\n";
}


int set_intersect(const std::unordered_set <std::string> & set_a, const std::unordered_set <std::string> & set_b) {
  int n = 0;

  if (set_a.size() < set_b.size()) {
    for (auto & item : set_a) {
      if (set_b.find(item) != set_b.end()) {
        n += 1;
      }
    }
  } else {
    for (auto & item : set_b) {
      if (set_a.find(item) != set_a.end()) {
        n += 1;
      }
    }
  }

  return n;
}


void dfs(const std::string & node, std::unordered_map <std::string, std::unordered_set <std::string>> & similarity_network,
         std::unordered_map <std::string, int> & node_state, std::unordered_set <std::string> & module) {
  if (node_state[node] == 0) {
    module.insert(similarity_network[node].begin(), similarity_network[node].end());
    node_state[node] = 1;
    for (auto & neighbor : similarity_network[node]) {
      dfs(neighbor, similarity_network, node_state, module);
    }
  }
}


struct edeg_similarity {
  std::string node_a;
  std::string node_b;
  double similarity;
};


void ThreadFunc(int n, int thread_num, std::vector <std::pair <std::string, std::string>> & edge_vec, 
  std::unordered_map <std::string, std::unordered_set <std::string>> & network, 
  std::vector <edeg_similarity> & similarity_vec) {

  thread_local std::vector <edeg_similarity> local_similarity_vec;
  int edge_num = edge_vec.size();
  int start = n * edge_num / thread_num;
	int stop = (n + 1) * edge_num / thread_num;
  for (int i = start; i < stop; ++i) {
    std::string node_a = edge_vec[i].first;
    std::string node_b = edge_vec[i].second;
    int node_a_degree = network[node_a].size();
    int node_b_degree = network[node_b].size();
    int shared = set_intersect(network[node_a], network[node_b]);
    double similarity = (shared + 1.0) / (node_a_degree + node_b_degree - shared - 1.0);
    local_similarity_vec.push_back(edeg_similarity {node_a, node_b, similarity});
  }

  mutex_lock.lock();
  similarity_vec.insert(similarity_vec.end(), local_similarity_vec.begin(), local_similarity_vec.end());
  mutex_lock.unlock();
}


int main(int argc, char* argv[]) {
  if (argc < 2) {
    module_identify_help();
    return -1;
  }

  // get option
  std::string in_file_name = "";
  std::string out_file_name = "";
  double similarity_cutoff = 0.5;
  int thread_num = 2;

  const char * const short_opts = "hvi:o:s:p:t:";
  const struct option long_opts[] =  {
    { "help", 0, NULL, 'h' },
    { "version", 0, NULL, 'v' },
    { "input", 1, NULL, 'i' },
    { "output", 1, NULL, 'o' },
    { "similarity", 1, NULL, 's' },
    { "pval", 1, NULL, 'p' },
    { "thread", 1, NULL, 't' },
    { NULL, 0, NULL, 0 }
  };
  int opt = getopt_long(argc, argv, short_opts, long_opts, NULL);
  while (opt != -1) {
    switch (opt) {
      case 'h':
        module_identify_help();
        return 0;
      case 'v':
        display_version();
        return 0;
      case 'o':
        out_file_name = optarg;
        break;
      case 'i':
        in_file_name = optarg;
        break;
      case 's':
        similarity_cutoff = std::stod(optarg);
        break;
      case 't':
        thread_num = std::stoi(optarg);
        break;
      case '?':
        module_identify_help();
        return 0;
      case -1:
        break;
      default:
        return -1;
    }
    opt = getopt_long( argc, argv, short_opts, long_opts, NULL );
  }

  // check options
  if (in_file_name.empty()) {
    std::cerr << "Error: -i/--input is required but not specified!\n";
    return -1;
  }

  if (out_file_name.empty()) {
    std::cerr << "Error: -o/--output is required but not specified!\n";
    return -1;
  }

  // open network file
  std::ifstream network_file(in_file_name, std::ios::in);
  if (!network_file.good()) {
    std::cerr << "Error while opening " << in_file_name << ".\n";
  }

  // read network file
  std::unordered_map <std::string, std::unordered_set <std::string>> network;
  std::vector <std::pair <std::string, std::string>> edge_vec;

  std::string line;
  while (getline(network_file, line)) {
    strim(line);
    std::stringstream line_stream;
    line_stream << line;
    std::string node_a, node_b;
    getline(line_stream, node_a, '\t');
    getline(line_stream, node_b, '\t');

    edge_vec.push_back(std::pair <std::string, std::string> {node_a, node_b});

    if (network.find(node_a) != network.end()) {
      network[node_a].insert(node_b);
    } else {      
      network[node_a] = std::unordered_set <std::string> {};
      network[node_a].insert(node_b);
    }

    if (network.find(node_b) != network.end()) {
      network[node_b].insert(node_a);
    } else {      
      network[node_b] = std::unordered_set <std::string> {};
      network[node_b].insert(node_a);
    }
  }

  // calc similarity
  std::vector <edeg_similarity> similarity_vec;
  std::vector <std::thread> threads;
  for (int i = 0; i < thread_num; ++i) {
    threads.push_back(std::thread{ ThreadFunc, i, thread_num, std::ref(edge_vec), std::ref(network), std::ref(similarity_vec)});
  }

  for (auto & t : threads) {
    t.join();
  }

  // initialize dfs
  std::unordered_map <std::string, std::unordered_set <std::string>> similarity_network;
  std::unordered_map <std::string, int> node_state;

  for (auto & item : similarity_vec) {
    if (item.similarity > similarity_cutoff) {
      std::string node_a = item.node_a;
      std::string node_b = item.node_b;
      if (similarity_network.find(node_a) != similarity_network.end()) {
        similarity_network[node_a].insert(node_b);
      } else {
        node_state[node_a] = 0;
        similarity_network[node_a] = std::unordered_set <std::string> {};
        similarity_network[node_a].insert(node_b);
      }
      if (similarity_network.find(node_b) != similarity_network.end()) {
        similarity_network[node_b].insert(node_a);
      } else {
        node_state[node_b] = 0;
        similarity_network[node_b] = std::unordered_set <std::string> {};
        similarity_network[node_b].insert(node_a);
      }
    }
  }

  // run dfs
  std::vector <std::unordered_set <std::string>> module_vec;
  for (auto & node_item : node_state) {
    std::string node = node_item.first;
    int state = node_item.second;
    if (state == 0) {
      std::unordered_set <std::string> module;
      dfs(node, similarity_network, node_state, module);
      module_vec.push_back(module);
    }
  }

  // output
  std::ofstream out_file(out_file_name, std::ios::out);
  if (!out_file.good()) {
    std::cerr << "Error while opening " << out_file_name << ".\n";
    return -1;
  }

  for (auto & module : module_vec) {
    bool flag = false;
    for (auto & node : module) {
      if (flag) {
        out_file << '\t';
      }
      out_file << node;
      flag = true;
    }    
    out_file << '\n';
  }

  return 0;
}
