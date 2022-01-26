#include <getopt.h>
#include <fstream>                   // std::ifstream, std::ofstream
#include <string>                    // std::string
#include "third_party/robin_hood.h"  // robin_hood::unordered_map, robin_hood::unordered_set
#include "util/base.hpp"  // version, strim(), split_string()

void network_extract_help() {
  std::cout << version;
  std::cout << "network_extract usage:\n";
  std::cout << "  network_extract -i input.network -g gene_list.txt -o "
               "output.network\n";
  std::cout << "options:\n";
  std::cout << "  -i --input <input network file>\n";
  std::cout << "  -o --output <output network file>\n";
  std::cout << "  -g --gene <gene list file>\n";
  std::cout << "  -v --version display GCEN version\n";
  std::cout << "  -h --help print help information\n";
  std::cout << "example:\n";
  std::cout << "  network_extract -i ../sample_data/gene_co_expr.network -g "
               "../sample_data/gene_list.txt -o "
               "../sample_data/sub.network\n";
}

int main(int argc, char* argv[]) {
  if (argc < 2) {
    network_extract_help();
    return 0;
  }

  // get option
  std::string in_file_name;
  std::string out_file_name;
  std::string gene_file_name;

  const char* const short_opts = "hi:o:g:";
  const struct option long_opts[] = {{"help", 0, NULL, 'h'},
                                     {"input", 1, NULL, 'i'},
                                     {"output", 1, NULL, 'o'},
                                     {"gene", 1, NULL, 'g'},
                                     {NULL, 0, NULL, 0}};
  int opt = getopt_long(argc, argv, short_opts, long_opts, NULL);
  while (opt != -1) {
    switch (opt) {
      case 'i':
        in_file_name = optarg;
        break;
      case 'o':
        out_file_name = optarg;
        break;
      case 'g':
        gene_file_name = optarg;
        break;
      case 'h':
        network_extract_help();
        return 0;
      case '?':
        network_extract_help();
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

  if (gene_file_name.empty()) {
    std::cerr << "Error: -g/--gene is required but not specified!\n";
    exit(-1);
  }

  // read gene list
  robin_hood::unordered_set<std::string> gene_set;
  std::ifstream gene_file(gene_file_name, std::ios::in);
  if (!gene_file.good()) {
    std::cerr << "Error while opening " << gene_file_name << ".\n";
  }
  std::string line;
  while (getline(gene_file, line)) {
    strim(line);
    if (line[0] == '#') {
      continue;
    }
    std::vector<std::string> str_vec;
    split_string(line, str_vec, "\t");
    std::string gene = str_vec[0];
    gene_set.insert(gene);
  }
  gene_file.close();

  // read input network and output
  std::ifstream in_file(in_file_name, std::ios::in);
  if (!in_file.good()) {
    std::cerr << "Error while opening " << in_file_name << ".\n";
  }
  std::ofstream out_file(out_file_name, std::ios::out);
  if (!out_file.good()) {
    std::cerr << "Error while opening " << out_file_name << ".\n";
    exit(-1);
  }

  while (getline(in_file, line)) {
    strim(line);
    if (line[0] == '#') {
      out_file << line << '\n';
    }
    std::vector<std::string> str_vec;
    split_string(line, str_vec, "\t");
    std::string node_a = str_vec[0];
    std::string node_b = str_vec[1];

    if ((gene_set.find(node_a) != gene_set.end()) ||
        (gene_set.find(node_b) != gene_set.end())) {
      out_file << line << '\n';
    }
  }

  in_file.close();
  out_file.close();

  return 0;
}
