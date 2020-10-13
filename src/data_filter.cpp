#include <getopt.h>
#include <iostream>
#include "util/func.hpp"
#include "util/strim.hpp"


void data_filter_help() {
  std::cout << "GCEN 0.5.0 by Wen Chen (chenwen@biochen.com, https://www.biochen.com/gcen)\n";
  std::cout << "data_filter usage:\n";
  std::cout << "  data_filter -i input_file -o output_file\n";
  std::cout << "options:\n";
  std::cout << "  -i --input <input file>\n";
  std::cout << "  -o --output <output file>\n";
  std::cout << "  -m --mean <number> mean cutoff (default: 0.0)\n";
  std::cout << "  -s --std <number> standard deviation cutoff (default: 0.0)\n";
  std::cout << "  -v --version display GCEN version\n";
  std::cout << "  -h --help print help information\n";
  std::cout << "example:\n";
  std::cout << "  data_filter -i ../sample_data/gene_expr.tsv -o ../sample_data/gene_expr_filter.tsv\n";
}


int main(int argc, char* argv[]) {
  if (argc < 2) {
    data_filter_help();
    return -1;
  }

  // get option
  std::string in_file_name = "";
  std::string out_file_name = "";
  double mean_cutoff = 0.0;
  double std_cutoff = 0.0;

  const char * const short_opts = "hvi:o:m:s:";
  const struct option long_opts[] = {
    { "help", 0, NULL, 'h' },
    { "version", 0, NULL, 'v' },
    { "input", 1, NULL, 'i' },
    { "output", 1, NULL, 'o' },
    { "mean", 1, NULL, 'm' },
    { "std", 1, NULL, 's' },
    { NULL, 0, NULL, 0 }
  };
  int opt = getopt_long(argc, argv, short_opts, long_opts, NULL);
  while (opt != -1) {
    switch (opt) {
      case 'h':
        data_filter_help();
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
      case 'm':
        mean_cutoff = std::stod(optarg);
        break;
      case 's':
        std_cutoff = std::stod(optarg);
        break;
      case '?':
        data_filter_help();
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

  // open file
  std::ifstream in_file(in_file_name, std::ios::in);
  if (!in_file.good()) {
    std::cerr << "Error while opening " << in_file_name << ".\n";
    return -1;
  }

  std::ofstream out_file(out_file_name, std::ios::out);
  if (!out_file.good()) {
    std::cerr << "Error while opening " << out_file_name << ".\n";
    return -1;
  }

  // filter
  std::string line;
  while (getline(in_file, line)) {
    strim(line);
    std::vector <std::string> str_vec;
    split_string(line, str_vec, "\t");
    try {
      std::vector <double> value_vector;
      for (unsigned int i = 1; i < str_vec.size(); ++i) {
        value_vector.push_back(std::stod(str_vec[i]));
      }
      if (mean(value_vector) >= mean_cutoff && deviation(value_vector) >= std_cutoff) {
        out_file << line << '\n';
      }
    } catch (std::invalid_argument) {
      continue;
    }
  }

  return 0;
}
