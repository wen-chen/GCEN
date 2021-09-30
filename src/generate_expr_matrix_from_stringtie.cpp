#include <getopt.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include "util/strim.hpp"
#include "util/func.hpp"


void stringtie_help() {
  std::cout << version;
  std::cout << "generate_expr_matrix_from_stringtie usage:\n";
  std::cout << "  generate_expr_matrix_from_stringtie -i input_file -o output_file\n";
  std::cout << "options:\n";
  std::cout << "  -i --input <input file> a text file with sample ID and path to its GTF file on each line\n";
  std::cout << "  -o --output <output file>\n";
  std::cout << "  -t --tpm output TMP value instead of FPKM vaule\n";
  std::cout << "  -v --version display GCEN version\n";
  std::cout << "  -h --help print help information\n";
  std::cout << "example:\n";
  std::cout << "  generate_expr_matrix_from_stringtie -i ../sample_data/stringtie/stringtie_sample.txt -o ../sample_data/stringtie/stringtie_gene_expr.tsv\n";
}


void read_sample_file(std::string sample_file_name, int i, std::string flag,
    std::unordered_map<std::string, std::unordered_set <std::string>> &gene_map,
    std::unordered_map<std::string, std::vector <std::string>> &transcript_map);


int main(int argc, char* argv[]) {
  if (argc < 2) {
    stringtie_help();
    return 0;
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
        stringtie_help();
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
        stringtie_help();
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
  std::unordered_map<std::string, std::unordered_set <std::string>> gene_map;
  std::unordered_map<std::string, std::vector <std::string>> transcript_map;
  
  for (unsigned int i = 0; i < sample_file_vec.size(); ++i) {
    std::string sample_file_name = sample_file_vec[i];
    read_sample_file(sample_file_name, i, flag, gene_map, transcript_map);
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
    std::string gene_id = item.first;
    std::unordered_set <std::string> transcript_set = item.second;
    std::vector <double> gene_expr_vec(sample_file_vec.size(), 0.0);
    for (auto transcript_id : transcript_set) {
      std::vector <std::string> transcript_expr_vec = transcript_map[transcript_id];
      for (unsigned int i = 0; i < transcript_expr_vec.size(); ++i) {
        double expr = std::stod(transcript_expr_vec[i]);
        gene_expr_vec[i] = gene_expr_vec[i] + expr;
      }
    }
    
    out_file << gene_id;
    for (auto expr : gene_expr_vec) {
      out_file << '\t' << expr;
    }
    out_file << '\n';
  }

  out_file.close();

  return 0;
}


std::string pretend_regex(const std::string &str, const std::string &pattern) {
  std::vector <std::string> str_vec_1;
  split_string(str, str_vec_1, ";");
  for (auto & str_1 : str_vec_1) {
    if (str_1.find(pattern) != str_1.npos) {
      std::vector <std::string> str_vec_2;
      split_string(str_1, str_vec_2, "\"");
      return str_vec_2[1];
    }
  }
  return "";
}


void push_record(int i, std::string &gene_id, std::string &expr_val, 
    std::unordered_map<std::string, std::vector <std::string>> &gene_map) {
  if (i > 0) { // check file consistency
    if (gene_map.find(gene_id) == gene_map.end()) {
      std::cerr << "Error: different gene size detected.";
      exit(-1);
    } else {
      if (gene_map[gene_id].size() != i) {
        std::cerr << "Error: different gene size detected.";
        exit(-1);
      }
    }
  } else {
    gene_map[gene_id] = std::vector <std::string> {};
  }
  gene_map[gene_id].push_back(expr_val);
}


void read_sample_file(std::string sample_file_name, int i, std::string flag,
    std::unordered_map<std::string, std::unordered_set <std::string>> &gene_map,
    std::unordered_map<std::string, std::vector <std::string>> &transcript_map) {
  std::ifstream sample_file(sample_file_name, std::ios::in);
  if (!sample_file.good()) {
    std::cerr << "Error while opening " << sample_file_name << ".\n";
    exit(-1);
  }

  std::string line;
  while(getline(sample_file, line)) {
    if (line.at(0) == '#') {
      continue;
    }

    strim(line);
    std::vector <std::string> str_vec;
    split_string(line, str_vec, "\t");
    std::string record_type = str_vec[2];
    std::string attribute = str_vec[8];
    
    if (record_type == "transcript") {
      std::string gene_id = pretend_regex(attribute, "gene_id");
      std::string transcript_id = pretend_regex(attribute, "transcript_id");

      std::string expr_val;
      if (flag == "FPKM") {
        expr_val = pretend_regex(attribute, "FPKM");
      } else if (flag == "TPM") {
        expr_val = pretend_regex(attribute, "TPM");
      }

      if (gene_id.empty()) {
        std::cerr << "Error while read " << sample_file_name << ".\n";
        exit(-1);
      }
      if (transcript_id.empty()) {
        std::cerr << "Error while read " << sample_file_name << ".\n";
        exit(-1);
      }
      if (expr_val.empty()) {
        std::cerr << "Error while read " << sample_file_name << ".\n";
        exit(-1);
      }

      push_record(i, transcript_id, expr_val, transcript_map);
      
      if (i == 0) {
        if (gene_map.find(gene_id) == gene_map.end()) {
          gene_map[gene_id] = std::unordered_set <std::string> {};
        }
        gene_map[gene_id].insert(transcript_id);
      }
    }
  }

  sample_file.close();
}
