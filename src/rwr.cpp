#include <getopt.h>
#include <algorithm> // sort rwr result
#include "util/func.hpp"
#include "util/rwr_util.hpp"


void rwr_help() {
  std::cout << "GCEN 0.5.1 by Wen Chen (chenwen@biochen.com, https://www.biochen.com/gcen)\n";
  std::cout << "rwr usage:\n";
  std::cout << "  rwr -n input_network -g gene_list -o output_result\n";
  std::cout << "options:\n";
  std::cout << "  -n --network <network file>\n";
  std::cout << "  -g --gene <seed genes list file>\n";
  std::cout << "  -r --gamma <number> restart probability (default: 0.5)\n";
  std::cout << "  -p --prob <number> probability cutoff (defalut: 0.01)\n";
  std::cout << "  -o --output <output file>\n";
  std::cout << "  -d --directed_network the input network is directed (defalut: the input network is undirected)\n";
  std::cout << "  -w --weighted_network the edge weights of network will be considered (defalut: all edge weights of network are set to 1.0)\n";
  std::cout << "  -W --weighted_gene the weights of seed genes will be considered (defalut: all weights of seed genes are set to 1.0)\n";
  std::cout << "  -v --version display GCEN version\n";
  std::cout << "  -h --help print help information\n";
  std::cout << "example:\n";
  std::cout << "  rwr -n ../sample_data/gene_co_expr.network -g ../sample_data/rwr_seed_genes.list -o ../sample_data/rwr_ranked_gene.tsv\n";
}


int main(int argc, char* argv[]) {
  if (argc < 2) {
    rwr_help();
    return -1;
  }

  // get option
  std::string network_file_name = "";
  std::string gene_file_name = "";
  std::string out_file_name = "";
  double gamma = 0.5;
  double prob_cutoff = 0.01;
  bool if_directed_network = false;
  bool if_weighted_network = false;
  bool if_weighted_gene = false;

  const char * const short_opts = "hvdWwn:g:r:p:o:";
  const struct option long_opts[] =  {
    { "help", 0, NULL, 'h' },
    { "version", 0, NULL, 'v' },
    { "directed_network", 0, NULL, 'd' },
    { "weighted_network", 0, NULL, 'w' },
    { "weighted_gene", 0, NULL, 'W' },
    { "network", 1, NULL, 'n' },
    { "gene", 1, NULL, 'g' },
    { "gamma", 1, NULL, 'r' },
    { "prob", 1, NULL, 'p' },
    { "output", 1, NULL, 'o'},
    { NULL, 0, NULL, 0 }
  };

  int opt = getopt_long(argc, argv, short_opts, long_opts, NULL);
  while (opt != -1) {
    switch (opt) {
      case 'h':
        rwr_help();
        return 0;
      case 'v':
        display_version();
        return 0;
      case 'd':
        if_directed_network = true;
        break;
      case 'w':
        if_weighted_network = true;
        break;
      case 'W':
        if_weighted_gene = true;
        break;
      case 'n':
        network_file_name = optarg;
        break;
      case 'g':
        gene_file_name = optarg;
        break;      
      case 'r':
        gamma = std::stod(optarg);
        break;
      case 'p':
        prob_cutoff = std::stod(optarg);
        break;
      case 'o':
        out_file_name = optarg;
        break;
      case '?':
        rwr_help();
        return 0;
      case -1:
        break;
      default:
        return -1;
    }
    opt = getopt_long( argc, argv, short_opts, long_opts, NULL );
  }

  // check options
  if (network_file_name.empty()) {
    std::cerr << "Error: -n/--network is required but not specified!\n";
    return -1;
  }

  if (gene_file_name.empty()) {
    std::cerr << "Error: -g/--gene is required but not specified!\n";
    return -1;
  }

  if (out_file_name.empty()) {
    std::cerr << "Error: -o/--output is required but not specified!\n";
    return -1;
  }

  // load network file
  std::vector <std::string> gene_vec;
  robin_hood::unordered_map <std::string, robin_hood::unordered_map <std::string, double>> network;
  load_network(network_file_name, if_directed_network, if_weighted_network, gene_vec, network);

  // network to matrix
  Dynamic_Matrix matrix(gene_vec.size(), gene_vec.size());
  network_2_matrix(gene_vec, network, matrix);

  // load gene list
  Dynamic_Vector P0(gene_vec.size());
  load_gene(gene_file_name, gene_vec, if_weighted_gene, P0);

  // Random Walk with Restart
  Dynamic_Vector PT(gene_vec.size());
  int step = rwr(matrix, P0, gamma, gene_vec, PT);

  // output
  std::cout << "Random walk " << step << " steps." << std::endl;

  std::vector <RWR_result> rwr_result_vec;
  for (int i = 0; i < gene_vec.size(); ++i) {
    if (PT[i] > prob_cutoff) {
      RWR_result rwr_result;
      rwr_result.gene = gene_vec[i];
      rwr_result.prob = PT[i];
      rwr_result_vec.push_back(rwr_result);
    }
  }

  std::sort(rwr_result_vec.begin(), rwr_result_vec.end());

  std::ofstream out_file(out_file_name, std::ios::out);
  if (!out_file.good()) {
    std::cerr << "Error while opening " << out_file_name << ".\n";
    return -1;
  }

  for (auto & rwr_result : rwr_result_vec) {
    out_file <<  rwr_result.gene << '\t' << rwr_result.prob  << '\n';
  }

  return 0;
}
