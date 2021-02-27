#include <getopt.h>
#include <string>
#include <unordered_map>
#include <fstream>
#include "util/func.hpp"


void network_merge_help() {
  std::cout << "GCEN 0.5.1 by Wen Chen (chenwen@biochen.com, https://www.biochen.com/gcen)\n";
  std::cout << "network_merge usage:\n";
  std::cout << "  network_merge -i input_files -o output_file\n";
  std::cout << "options:\n";
  std::cout << "  -i --input <input files> multiple files are separated by commas\n";
  std::cout << "  -o --output <output file>\n";
  std::cout << "  -c --cor <number> correlation coefficient cutoff (default: 0.5)\n";
  std::cout << "  -h --help print help information\n";
  std::cout << "example:\n";
  std::cout << "  network_merge -i ../sample_data/test_1.network,../sample_data/test_2.network -o ../sample_data/test_merge.network\n";
}


int main(int argc, char* argv[]) {
  if (argc < 2) {
    network_merge_help();
    return -1;
  }

  // get option
  std::string in_file_names = "";
  std::string out_file_name = "";
  double cor_cutoff = 0.5;

  const char * const short_opts = "hi:o:c:";
  const struct option long_opts[] = {
    { "help", 0, NULL, 'h' },
    { "input", 1, NULL, 'i' },
    { "output", 1, NULL, 'o' },
    { "cor", 1, NULL, 'c' },
    { NULL, 0, NULL, 0 }
  };
  int opt = getopt_long(argc, argv, short_opts, long_opts, NULL);
  while (opt != -1) {
    switch (opt) {
      case 'h':
        network_merge_help();
        return 0;
      case 'o':
        out_file_name = optarg;
        break;
      case 'i':
        in_file_names = optarg;
        break;
      case 'c':
        cor_cutoff = std::stod(optarg);
        break;
      case '?':
        network_merge_help();
        return 0;
      case -1:
        break;
      default:
        return -1;
    }
    opt = getopt_long( argc, argv, short_opts, long_opts, NULL );		
  }

  if (in_file_names.empty()) {
    std::cerr << "Error: -i/--input is required but not specified!\n";
    return -1;
  }

  if (out_file_name.empty()) {
    std::cerr << "Error: -o/--output is required but not specified!\n";
    return -1;
  }

  // split multiple input file
  std::vector <std::string> in_file_name_vec;
  split_string(in_file_names, in_file_name_vec, ",");

  // read network
  std::unordered_map <std::string, std::unordered_map <std::string, std::vector<double>>> network;
  for (std::string& in_file_name : in_file_name_vec) {
    std::ifstream network_file(in_file_name, std::ios::in);
    if (!network_file.good()) {
      std::cerr << "Error while opening " << in_file_name << ".\n";
      return -1;
    }

    std::string line;
    while (getline(network_file, line)) {
      strim(line);
      if (line[0] == '#') {
        continue;
      }
      std::vector <std::string> str_vec;
      split_string(line, str_vec, "\t");
      std::string gene_a = str_vec[0];
      std::string gene_b = str_vec[1];
      double cor = std::stod(str_vec[2]);

      if (gene_a.compare(gene_b) > 0) {
        std::string tmp = gene_a;
        gene_a = gene_b;
        gene_b = tmp;
      }

      if (network.find(gene_a) != network.end() and network[gene_a].find(gene_b) != network[gene_a].end()) {
        network[gene_a][gene_b].push_back(cor);
      } else {
        network[gene_a][gene_b] = std::vector <double> {cor};
      }
    }
  }

  // write network
  std::ofstream out_file(out_file_name, std::ios::out);
  if (!out_file.good()) {
    std::cerr << "Error while opening " << out_file_name << ".\n";
    return -1;
  }

  for (auto const & x : network) {
    for (auto const & y : x.second) {
      double positive_cor = 0;
      double negative_cor = 0;
      for (auto const & z : y.second) {
        if (z < 0.0) {
          negative_cor = negative_cor + z;
        } else {
          positive_cor = positive_cor + z;
        }
      }

      if (positive_cor >= std::abs(negative_cor) and positive_cor >= cor_cutoff) {
        out_file << x.first << "\t"  << y.first << "\t" << positive_cor << "1" << "\n";
      } else if (std::abs(negative_cor) > positive_cor and std::abs(negative_cor) >= cor_cutoff) {
        out_file << x.first << "\t"  << y.first << "\t" << negative_cor << "-1" << "\n";
      }
    }
  }

  return 0;
}
