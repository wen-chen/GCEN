#include <getopt.h>
#include <vector>
#include <unordered_set>
#include "./util/base.hpp"
#include "./util/norm.hpp"


void data_norm_help() {
  std::cout << version;
  std::cout << "data_norm usage:\n";
  std::cout << "  data_norm -i input_file -o output_file -m normalization_method\n";
  std::cout << "options:\n";
  std::cout << "  -i --input <input file>\n";
  std::cout << "  -o --output <output file>\n";
  std::cout << "  -m --method <upqt | median | deseq | tmm | hkg> normalization method (default: upqt)\n";
  std::cout << "  -g --gene <housekeeping gene list>  only for '--method hkg'\n";
  std::cout << "  -v --version display GCEN version\n";
  std::cout << "  -h --help print help information\n";
  std::cout << "example:\n";
  std::cout << "  data_norm -i ../sample_data/gene_expr.tsv -o ../sample_data/gene_expr_norm.tsv -m tmm\n";
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
  std::string housekeeping_gene_file_name = "";

  const char * const short_opts = "hvi:o:m:g:";
  const struct option long_opts[] = {
    { "help", 0, NULL, 'h' },
    { "version", 0, NULL, 'v' },
    { "input", 1, NULL, 'i' },
    { "output", 1, NULL, 'o' },
    { "method", 1, NULL, 'm' },
    { "gene", 1, NULL, 'g' },
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
      case 'g':
        housekeeping_gene_file_name = optarg;
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

  if (!((method == "upqt") || (method == "median") || (method == "deseq") || (method == "tmm") || (method == "hkg"))) {
    std::cerr << "Error: unrecognized -m/--method parameter!\n";
    return -1;
  }

  if ((method == "hkg") && (housekeeping_gene_file_name.empty())) {
    std::cerr << "Error: -g/--gene is required but not specified!\n";
    return -1;
  }

  // load file
  std::vector <std::string> annotation_vec;
  std::vector <std::string> gene_name_vec;
  std::vector <std::vector <double>> expr_matrix;
  load(in_file_name, annotation_vec, gene_name_vec, expr_matrix);

  // norm
  if (method == "upqt") {
    norm_upqt(expr_matrix);
  } else if (method == "median") {
    norm_median(expr_matrix);
  } else if (method == "deseq") {
    norm_deseq(expr_matrix);
  } else if (method == "tmm") {
    norm_tmm(expr_matrix);
  } else if (method == "hkg") {
    std::ifstream housekeeping_gene_file(housekeeping_gene_file_name, std::ios::in);
    if (!housekeeping_gene_file.good()) {
      std::cerr << "Error while opening " << housekeeping_gene_file_name << ".\n";
      exit(-1);
    }
    
    std::unordered_set <std::string> hkg_set;
    std::string line;
    while (getline(housekeeping_gene_file, line)) {
      strim(line);
      hkg_set.insert(line);
    }

    housekeeping_gene_file.close();

    std::vector <bool> hkg_vec(gene_name_vec.size(), false);
    int hkg_num = 0;
    for (unsigned int i = 0; i < hkg_vec.size(); ++i) {
      if (hkg_set.find(gene_name_vec[i]) != hkg_set.end()) {
        hkg_vec[i] = true;
        hkg_num += 1;
      }
    }

    if (hkg_num ==0) {
      std::cerr << "Error: None of the housekeeping genes were identified!\n";
      exit(-1);
    }

    norm_hkg(expr_matrix, hkg_vec);
  }

  // output
  std::ofstream out_file(out_file_name, std::ios::out);
  if (!out_file.good()) {
    std::cerr << "Error while opening " << out_file_name << ".\n";
    return -1;
  }

  for (auto & line : annotation_vec) {
    out_file << line << '\n';
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

  out_file.close();

  double mean_of_mean = mean(mean_vec);
  double sd_of_mean = deviation(mean_vec);
  std::cout << "The average gene expression mean is " << mean_of_mean << " +- " << sd_of_mean << ".\n";

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
