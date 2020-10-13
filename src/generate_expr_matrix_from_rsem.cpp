#include <getopt.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include "util/strim.hpp"
#include "util/func.hpp"


void rsem_help() {
  std::cout << "GCEN 0.5.0 by Wen Chen (chenwen@biochen.com, https://www.biochen.com/gcen)\n";
  std::cout << "generate_expr_matrix_from_rsem usage:\n";
  std::cout << "  generate_expr_matrix_from_rsem -i input_file -o output_file\n";
  std::cout << "options:\n";
  std::cout << "  -i --input <input file> a text file with sample ID and path to its RSEM result file on each line\n";
  std::cout << "  -o --output <output file>\n";
  std::cout << "  -t --tpm output TPM value instead of FPKM vaule\n";
  std::cout << "  -v --version display GCEN version\n";
  std::cout << "  -h --help print help information\n";
  std::cout << "example:\n";
  std::cout << "  generate_expr_matrix_from_rsem -i ../sample_data/rsem/rsem_sample.txt -o ../sample_data/rsem/rsem_gene_expr.tsv\n";
}


void read_sample_file(std::string sample_file_name, int i, std::string flag,
    std::unordered_map<std::string, std::vector <std::string>> &gene_map);


int main(int argc, char* argv[]) {
  if (argc < 2) {
    rsem_help();
    return -1;
  }

  // get option
  std::string in_file_name = "";
  std::string out_file_name = "";
  std::string flag = "FPKM";

  const char * const short_opts = "hvti:o:";
  const struct option long_opts[] = {
    { "help", 0, NULL, 'h' },
    { "version", 0, NULL, 'v' },
    { "input", 1, NULL, 'i' },
    { "output", 1, NULL, 'o' },
    { "tpm", 0, NULL, 't' },
    { NULL, 0, NULL, 0 }
  };
  int opt = getopt_long(argc, argv, short_opts, long_opts, NULL);
  while (opt != -1) {
    switch (opt) {
      case 'h':
        rsem_help();
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
      case 't':
        flag = "TPM";
        break;
      case '?':
        rsem_help();
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

  // read input file
  std::vector <std::string> sample_id_vec;
  std::vector <std::string> sample_file_vec;

  std::ifstream in_file(in_file_name, std::ios::in);
  if (!in_file.good()) {
    std::cerr << "Error while opening " << in_file_name << ".\n";
    exit(-1);
  }

  std::string line;
  while(getline(in_file, line)) {
    strim(line);
    std::vector <std::string> str_vec;
    split_string(line, str_vec, "\t");
    std::string sample_id = str_vec[0];
    std::string sample_file_name = str_vec[1];
    sample_id_vec.push_back(sample_id);
    sample_file_vec.push_back(sample_file_name);
  }

  in_file.close();

  // read sample files
  std::unordered_map<std::string, std::vector <std::string>> gene_map;
  for (unsigned int i = 0; i < sample_file_vec.size(); ++i) {
    std::string sample_file_name = sample_file_vec[i];
    read_sample_file(sample_file_name, i, flag, gene_map);
  }

  // output
  std::ofstream out_file(out_file_name, std::ios::out);
  if (!out_file.good()) {
    std::cerr << "Error while opening " << out_file_name << ".\n";
    exit(-1);
  }

  join_vector(sample_id_vec, "\t", line);
  out_file << "#gene\t" << line << '\n';

  for (auto item : gene_map) {
    std::string gene = item.first;
    std::vector <std::string> expr = item.second;
    join_vector(expr, "\t", line);
    out_file << gene << '\t' << line << '\n';
  }

  out_file.close();

  return 0;
}


void read_sample_file(std::string sample_file_name, int i, std::string flag,
    std::unordered_map<std::string, std::vector <std::string>> &gene_map) {
  std::ifstream sample_file(sample_file_name, std::ios::in);
  if (!sample_file.good()) {
    std::cerr << "Error while opening " << sample_file_name << ".\n";
    exit(-1);
  }

  std::string line;
  getline(sample_file, line); // skip the first line
  while(getline(sample_file, line)) {
    strim(line);
    std::vector <std::string> str_vec;
    split_string(line, str_vec, "\t");
    std::string gene = str_vec[0];
    std::string TPM = str_vec[5];
    std::string FPKM = str_vec[6];

    if (i > 0) { // check file consistency
      if (gene_map.find(gene) == gene_map.end()) {
        std::cerr << "Error: different gene size detected!";
        exit(-1);
      } else {
        if (gene_map[gene].size() != i) {
          std::cerr << "Error: different gene size detected!";
          exit(-1);
        }
      }
    } else {
      gene_map[gene] = std::vector <std::string> {};
    }

    if (flag == "FPKM") {
      gene_map[gene].push_back(FPKM);
    } else if (flag == "TPM") {
      gene_map[gene].push_back(TPM);
    }
  }

  sample_file.close();
}
