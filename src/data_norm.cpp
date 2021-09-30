#include <getopt.h>
#include <vector>
#include "./util/func.hpp"
#include "./util/norm_util.hpp"


void data_norm_help() {
  std::cout << version;
  std::cout << "data_norm usage:\n";
  std::cout << "  data_norm -i input_file -o output_file -m normalization_method\n";
  std::cout << "options:\n";
  std::cout << "  -i --input <input file>\n";
  std::cout << "  -o --output <output file>\n";
  std::cout << "  -m --method <upqt/median/deseq/tmm> normalization method (default: upqt)\n";
  std::cout << "  -v --version display GCEN version\n";
  std::cout << "  -h --help print help information\n";
  std::cout << "example:\n";
  std::cout << "  data_norm -i ../sample_data/gene_expr.tsv -o ../sample_data/gene_expr_norm.tsv -m upqt\n";
}


int main(int argc, char* argv[]) {
  if (argc < 2) {
    data_norm_help();
    return 0;
  }

  // get option
  std::string in_file_name = "";
  std::string out_file_name = "";
  std::string method = "upqt";

  const char * const short_opts = "hvi:o:m:";
  const struct option long_opts[] = {
    { "help", 0, NULL, 'h' },
    { "version", 0, NULL, 'v' },
    { "input", 1, NULL, 'i' },
    { "output", 1, NULL, 'o' },
    { "method", 1, NULL, 'm' },
    { NULL, 0, NULL, 0 }
  };
  int opt = getopt_long(argc, argv, short_opts, long_opts, NULL);
  while (opt != -1) {
    switch (opt) {
      case 'h':
        data_norm_help();
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
        method = optarg;
        break;
      case '?':
        data_norm_help();
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

  if (!((method == "upqt") || (method == "median") || (method == "deseq") || (method == "tmm"))) {
    std::cerr << "Error: unrecognized -m/--method parameter!\n";
    return -1;
  }

  // load file
  std::vector <std::string> gene_name_vec;
  std::vector <std::vector <double>> expr_matrix;
  load(in_file_name, gene_name_vec, expr_matrix);

  // norm
  if (method == "upqt") {
    norm_upqt(expr_matrix);
  } else if (method == "median") {
    norm_median(expr_matrix);
  } else if (method == "deseq") {
    norm_deseq(expr_matrix);
  } else if (method == "deseq") {
    norm_tmm(expr_matrix);
  }

  // output
  std::ofstream out_file(out_file_name, std::ios::out);
  if (!out_file.good()) {
    std::cerr << "Error while opening " << out_file_name << ".\n";
    return -1;
  }
  
  std::vector <double> mean_vec;
  std::vector <double> sd_vec;

  for (unsigned int i = 0; i < expr_matrix.size(); ++i) {
    out_file << gene_name_vec[i];
    mean_vec.push_back(mean(expr_matrix[i]));
    sd_vec.push_back(deviation(expr_matrix[i]));
    for (unsigned int j = 0; j < expr_matrix[i].size(); ++j) {
      out_file << '\t' << expr_matrix[i][j];
    }
    out_file << '\n';
  }

  std::sort(mean_vec.begin(), mean_vec.end());
  double per5_mean = mean_vec[int (mean_vec.size() * 0.05)];
  double per10_mean = mean_vec[int (mean_vec.size() * 0.10)];
  double per20_mean = mean_vec[int (mean_vec.size() * 0.20)];
  double per25_mean = mean_vec[int (mean_vec.size() * 0.25)];
  double per33_mean = mean_vec[int (mean_vec.size() * 0.33)];
  double per50_mean = mean_vec[int (mean_vec.size() * 0.50)];
  std::cout << "(5%: " << per5_mean << "; ";
  std::cout << "10%: " << per10_mean << "; ";
  std::cout << "20%: " << per20_mean << "; ";
  std::cout << "25%: " << per25_mean << "; ";
  std::cout << "33%: " << per33_mean << "; ";
  std::cout << "50%: " << per50_mean << ").\n";

  double mean_of_sd = mean(sd_vec);
  double sd_of_sd = deviation(sd_vec);
  std::cout << "The average gene expression sd is " << mean_of_sd << " +- " << sd_of_sd << ".\n";

  std::sort(sd_vec.begin(), sd_vec.end());
  double per5_sd = sd_vec[int (sd_vec.size() * 0.05)];
  double per10_sd = sd_vec[int (sd_vec.size() * 0.10)];
  double per20_sd = sd_vec[int (sd_vec.size() * 0.20)];
  double per25_sd = sd_vec[int (sd_vec.size() * 0.25)];
  double per33_sd = sd_vec[int (sd_vec.size() * 0.33)];
  double per50_sd = sd_vec[int (sd_vec.size() * 0.50)];
  std::cout << "(5%: " << per5_sd << "; ";
  std::cout << "10%: " << per10_sd << "; ";
  std::cout << "20%: " << per20_sd << "; ";
  std::cout << "25%: " << per25_sd << "; ";
  std::cout << "33%: " << per33_sd << "; ";
  std::cout << "50%: " << per50_sd << ").\n";

  return 0;
}
