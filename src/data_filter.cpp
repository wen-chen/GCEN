#include <getopt.h>
#include <iostream>
#include <algorithm>
#include "util/func.hpp"
#include "util/strim.hpp"


void data_filter_help() {
  std::cout << "GCEN 0.5.1 by Wen Chen (chenwen@biochen.com, https://www.biochen.com/gcen)\n";
  std::cout << "data_filter usage:\n";
  std::cout << "  data_filter -i input_file -o output_file\n";
  std::cout << "options:\n";
  std::cout << "  -i --input <input file>\n";
  std::cout << "  -o --output <output file>\n";
  std::cout << "  -c --cutoff_mean <number> mean cutoff of gene expression (default: 0.0)\n";
  std::cout << "  -C --cutoff_sd <number> standard deviation cutoff of gene expression (default: 0.0)\n";
  std::cout << "  -p --percent_mean <number> keep a proportion of total genes based mean of gene expression (default: 1.0)\n";
  std::cout << "  -P --percent_sd <number> keep a proportion of total genes based standard deviation of gene expression (default: 1.0)\n";
  std::cout << "  -v --version display GCEN version\n";
  std::cout << "  -h --help print help information\n";
  std::cout << "example:\n";
  std::cout << "  data_filter -i ../sample_data/gene_expr.tsv -o ../sample_data/gene_expr_filter.tsv -p 0.75\n";
}


int main(int argc, char* argv[]) {
  if (argc < 2) {
    data_filter_help();
    return -1;
  }

  // get option
  std::string in_file_name = "";
  std::string out_file_name = "";
  double cutoff_mean = 0.0;
  double cutoff_sd = 0.0;
  double percent_mean = 1.0;
  double percent_sd = 1.0;

  const char * const short_opts = "hvi:o:c:C:p:P";
  const struct option long_opts[] = {
    { "help", 0, NULL, 'h' },
    { "version", 0, NULL, 'v' },
    { "input", 1, NULL, 'i' },
    { "output", 1, NULL, 'o' },
    { "cutoff_mean", 1, NULL, 'c' },
    { "cutoff_sd", 1, NULL, 'C' },
    { "percent_mean", 1, NULL, 'p' },
    { "percent_sd", 1, NULL, 'P' },
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
      case 'c':
        cutoff_mean = std::stod(optarg);
        break;
      case 'C':
        cutoff_sd = std::stod(optarg);
        break;
      case 'p':
        percent_mean = std::stod(optarg);
        break;
      case 'P':
        percent_sd = std::stod(optarg);
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

  // stat
  std::vector <double> mean_vec;
  std::vector <double> sd_vec;

  std::ifstream in_file(in_file_name, std::ios::in);
  if (!in_file.good()) {
    std::cerr << "Error while opening " << in_file_name << ".\n";
    return -1;
  }

  std::string line;
  while (getline(in_file, line)) {
    strim(line);
    if (line[0] == '#') {
      continue;
    }
    std::vector <std::string> str_vec;
    split_string(line, str_vec, "\t");
    try {
      std::vector <double> val_vec;
      for (unsigned int i = 1; i < str_vec.size(); ++i) {
        val_vec.push_back(std::stod(str_vec[i]));
      }
      mean_vec.push_back(mean(val_vec));
      sd_vec.push_back(deviation(val_vec));
    } catch (std::invalid_argument) {
      continue;
    }
  }

  in_file.close();

  // determine the cutoff
  std::sort(mean_vec.begin(), mean_vec.end());
  if (cutoff_mean < mean_vec[int (mean_vec.size() * (1.0 - percent_mean))]) {
    cutoff_mean = mean_vec[int (mean_vec.size() * (1.0 - percent_mean))];
  }

  std::sort(sd_vec.begin(), sd_vec.end());
  if (cutoff_sd < sd_vec[int (sd_vec.size() * (1.0 - percent_sd))]) {
    cutoff_sd = sd_vec[int (sd_vec.size() * (1.0 - percent_sd))];
  }

  // filter
  std::ifstream in_file_2(in_file_name, std::ios::in);
  if (!in_file_2.good()) {
    std::cerr << "Error while opening " << in_file_name << ".\n";
    return -1;
  }

  std::ofstream out_file(out_file_name, std::ios::out);
  if (!out_file.good()) {
    std::cerr << "Error while opening " << out_file_name << ".\n";
    return -1;
  }

  while (getline(in_file_2, line)) {
    strim(line);
    if (line[0] == '#') {
      continue;
    }
    std::vector <std::string> str_vec;
    split_string(line, str_vec, "\t");
    try {
      std::vector <double> value_vector;
      for (unsigned int i = 1; i < str_vec.size(); ++i) {
        value_vector.push_back(std::stod(str_vec[i]));
      }
      if (mean(value_vector) >= cutoff_mean && deviation(value_vector) >= cutoff_sd) {
        out_file << line << '\n';
      }
    } catch (std::invalid_argument) {
      continue;
    }
  }

  in_file_2.close();
  out_file.close();

  return 0;
}
