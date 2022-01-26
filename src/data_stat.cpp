#include <getopt.h>
#include <algorithm>
#include <iostream>
#include "util/base.hpp"
#include "util/strim.hpp"

void data_stat_help() {
  std::cout << version;
  std::cout << "data_stat usage:\n";
  std::cout << "  data_stat -i input_file\n";
  std::cout << "options:\n";
  std::cout << "  -i --input <input file>\n";
  std::cout << "  -v --version display GCEN version\n";
  std::cout << "  -h --help print help information\n";
  std::cout << "example:\n";
  std::cout << "  data_stat -i ../sample_data/gene_expr.tsv\n";
}

int main(int argc, char* argv[]) {
  if (argc < 2) {
    data_stat_help();
    return 0;
  }

  // get option
  std::string in_file_name = "";

  const char* const short_opts = "hvi:o:m:s:";
  const struct option long_opts[] = {{"help", 0, NULL, 'h'},
                                     {"version", 0, NULL, 'v'},
                                     {"input", 1, NULL, 'i'},
                                     {NULL, 0, NULL, 0}};
  int opt = getopt_long(argc, argv, short_opts, long_opts, NULL);
  while (opt != -1) {
    switch (opt) {
      case 'h':
        data_stat_help();
        return 0;
      case 'v':
        display_version();
        return 0;
      case 'i':
        in_file_name = optarg;
        break;
      case '?':
        data_stat_help();
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

  // open file
  std::ifstream in_file(in_file_name, std::ios::in);
  if (!in_file.good()) {
    std::cerr << "Error while opening " << in_file_name << ".\n";
    exit(-1);
  }

  // stat
  std::string line;
  std::vector<double> mean_vec;
  std::vector<double> sd_vec;
  while (getline(in_file, line)) {
    strim(line);
    if (line[0] == '#') {
      continue;
    }
    std::vector<std::string> str_vec;
    split_string(line, str_vec, "\t");
    try {
      std::vector<double> val_vec;
      for (unsigned int i = 1; i < str_vec.size(); ++i) {
        val_vec.push_back(std::stod(str_vec[i]));
      }
      mean_vec.push_back(mean(val_vec));
      sd_vec.push_back(deviation(val_vec));
    } catch (std::invalid_argument&) {
      continue;
    }
  }

  in_file.close();

  double mean_of_mean = mean(mean_vec);
  double sd_of_mean = deviation(mean_vec);
  std::cout << "The average gene expression mean is " << mean_of_mean << " +- "
            << sd_of_mean << ".\n";

  std::sort(mean_vec.begin(), mean_vec.end());
  double per5_mean = mean_vec[int(mean_vec.size() * 0.05)];
  double per10_mean = mean_vec[int(mean_vec.size() * 0.10)];
  double per20_mean = mean_vec[int(mean_vec.size() * 0.20)];
  double per25_mean = mean_vec[int(mean_vec.size() * 0.25)];
  double per33_mean = mean_vec[int(mean_vec.size() * 0.33)];
  double per50_mean = mean_vec[int(mean_vec.size() * 0.50)];
  std::cout << "(5%: " << per5_mean << "; ";
  std::cout << "10%: " << per10_mean << "; ";
  std::cout << "20%: " << per20_mean << "; ";
  std::cout << "25%: " << per25_mean << "; ";
  std::cout << "33%: " << per33_mean << "; ";
  std::cout << "50%: " << per50_mean << ").\n";

  double mean_of_sd = mean(sd_vec);
  double sd_of_sd = deviation(sd_vec);
  std::cout << "The average gene expression sd is " << mean_of_sd << " +- "
            << sd_of_sd << ".\n";

  std::sort(sd_vec.begin(), sd_vec.end());
  double per5_sd = sd_vec[int(sd_vec.size() * 0.05)];
  double per10_sd = sd_vec[int(sd_vec.size() * 0.10)];
  double per20_sd = sd_vec[int(sd_vec.size() * 0.20)];
  double per25_sd = sd_vec[int(sd_vec.size() * 0.25)];
  double per33_sd = sd_vec[int(sd_vec.size() * 0.33)];
  double per50_sd = sd_vec[int(sd_vec.size() * 0.50)];
  std::cout << "(5%: " << per5_sd << "; ";
  std::cout << "10%: " << per10_sd << "; ";
  std::cout << "20%: " << per20_sd << "; ";
  std::cout << "25%: " << per25_sd << "; ";
  std::cout << "33%: " << per33_sd << "; ";
  std::cout << "50%: " << per50_sd << ").\n";

  return 0;
}
