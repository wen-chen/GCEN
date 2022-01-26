#include <getopt.h>

#include <algorithm>  // std::replace

#include "util/base.hpp"

void csv_to_tsv_help() {
  std::cout << version;
  std::cout << "csv_to_tsv usage:\n";
  std::cout << "  csv_to_tsv -i input.csv -o output.tsv\n";
  std::cout << "options:\n";
  std::cout << "  -i --input <input csv file>\n";
  std::cout << "  -o --output <output tsv file>\n";
  std::cout << "  -v --version display GCEN version\n";
  std::cout << "  -h --help print help information\n";
  std::cout << "example:\n";
  std::cout << "  csv_to_tsv -i ../sample_data/gene_expr.csv -o "
               "../sample_data/gene_expr.tsv\n";
}

int main(int argc, char* argv[]) {
  if (argc < 2) {
    csv_to_tsv_help();
    return 0;
  }

  // get option
  char input_separator = ',';
  char output_separator = '\t';

  std::string in_file_name = "";
  std::string out_file_name = "";
  const char* const short_opts = "hvi:o:";
  const struct option long_opts[] = {{"help", 0, NULL, 'h'},
                                     {"version", 0, NULL, 'v'},
                                     {"input", 1, NULL, 'i'},
                                     {"output", 1, NULL, 'o'},
                                     {NULL, 0, NULL, 0}};
  int opt = getopt_long(argc, argv, short_opts, long_opts, NULL);
  while (opt != -1) {
    switch (opt) {
      case 'h':
        csv_to_tsv_help();
        return 0;
      case 'v':
        display_version();
        return 0;
      case 'i':
        in_file_name = optarg;
        break;
      case 'o':
        out_file_name = optarg;
        break;
      case '?':
        csv_to_tsv_help();
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

  // tsv to csv
  std::ifstream in_file(in_file_name, std::ios::in);
  if (!in_file.good()) {
    std::cerr << "Error while opening " << in_file_name << ".\n";
    exit(-1);
  }

  std::ofstream out_file(out_file_name, std::ios::out);
  if (!out_file.good()) {
    std::cerr << "Error while opening " << out_file_name << ".\n";
    exit(-1);
  }

  std::string line;
  while (getline(in_file, line)) {
    strim(line);
    if (line[0] == '#') {
      out_file << line << '\n';
      continue;
    }
    std::replace(line.begin(), line.end(), input_separator, output_separator);
    out_file << line << '\n';
  }

  in_file.close();
  out_file.close();

  return 0;
}
