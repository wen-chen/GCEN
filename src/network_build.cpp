#include <getopt.h>
#include <iostream>
#include <fstream> 
#include <sstream>
#include <thread>
#include <mutex>
#include <utility>
#include <ctime>
#include <chrono>
#include "util/Pearson.hpp"
#include "util/func.hpp"


void network_build_help() {
  std::cout << "GCEN 0.5.0 by Wen Chen (chenwen@biochen.com, https://www.biochen.com/gcen)\n";
  std::cout << "network_build usage:\n";
  std::cout << "  network_build -i gene_expression_file -o co_expression_network_file\n";
  std::cout << "options:\n";
  std::cout << "  -i --input <input file>\n";
  std::cout << "  -o --output <output file>\n";
  std::cout << "  -m --method <pearson or spearman> correlation coefficient method (default: spearman)\n";
  std::cout << "  -l --log <log2 or log10> make a log transformation (default: not transform)\n";
  std::cout << "  -t --thread <number> cpu cores (default: 2)\n";
  std::cout << "  -p --pval <number> p value cutoff (default: 0.001)\n";
  std::cout << "  -c --cor <number> correlation coefficient cutoff (default: 0.1)\n";
  std::cout << "  -s --signed <y or n> singed network (default: n)\n";
  std::cout << "  -f --fdr <y or n> calculate FDR (default: n)\n";
  std::cout << "  -v --version display GCEN version\n";
  std::cout << "  -h --help print help information\n";
  std::cout << "example:\n";
  std::cout << "  network_build -i ../sample_data/gene_expr.tsv -o ../sample_data/gene_co_expr.network -m spearman -p 0.001 -f y\n";
}


void ThreadFunc(int n, int thread_num, int GeneNum, int SampleNum, std::vector <std::vector <double>> & GeneDataFrame,
    std::vector <std::string> & GeneNameVector, double PvalCutoff, double CorCutoff, bool signed_network,
    std::ofstream & OutFile);


void ThreadFunc_FDR(int n, int thread_num, int GeneNum, int SampleNum, std::vector <std::vector <double>> & GeneDataFrame,
    std::vector <std::string> & GeneNameVector, double PvalCutoff, double CorCutoff, bool signed_network,
    std::vector <std::pair <int, int>> & node_id_pairs, std::vector <double> & corrs, std::vector <double> & p_values);


int main(int argc, char* argv[]) {
  if (argc < 2) {
    network_build_help();
    return -1;
  }

  // get option
  std::string in_file_name = "";
  std::string out_file_name = "";
  std::string method = "spearman";
  std::string log_transform = "";
  int thread_num = 2;
  double CorCutoff = 0.1;
  double PvalCutoff = 0.001;
  std::string fdr_opt = "";
  std::string signed_opt = "";
    
  const char * const short_opts = "hvi:o:m:l:t:p:c:f:s:";
  const struct option long_opts[] =  {
    { "help", 0, NULL, 'h' },
    { "version", 0, NULL, 'v' },
    { "input", 1, NULL, 'i' },
    { "output", 1, NULL, 'o' },
    { "method", 1, NULL, 'm' },
    { "log", 1, NULL, 'l' },
    { "thread", 1, NULL, 't' },
    { "pval", 1, NULL, 'p' },
    { "cor", 1, NULL, 'c' },
    { "fdr", 1, NULL, 'f' },
    { "signed", 1, NULL, 's' },
    { NULL, 0, NULL, 0 }
  };
  int opt = getopt_long(argc, argv, short_opts, long_opts, NULL);
  while (opt != -1) {
    switch (opt) {
      case 'h':
        network_build_help();
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
      case 'l':
        log_transform = optarg;
        break;
      case 't':
        thread_num = std::stoi(optarg);
        break;
      case 'p':
        PvalCutoff = std::stod(optarg);
        break;
      case 'c':
        CorCutoff = std::stod(optarg);
        break;
      case 'f':
        fdr_opt = optarg;
        break;
      case 's':
        signed_opt = optarg;
        break;
      case '?':
        network_build_help();
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

  bool fdr = false;
  if (fdr_opt == "y" or fdr_opt == "Y") {
    fdr = true;
  }

  bool if_log2 = false;
  bool if_log10 = false;
  if (log_transform == "log2") {
    if_log2 = true;
  } else if (log_transform == "log10") {
    if_log10 = true;
  }

  bool signed_network = false;
  if (signed_opt == "y" or signed_opt == "Y") {
    signed_network = true;
  }

  // read file	
  std::vector <std::string> GeneNameVector;
  std::vector <std::vector <double> > GeneDataFrame;
  load(in_file_name, GeneNameVector, GeneDataFrame, if_log2, if_log10);

  // get ranks
  int GeneNum = GeneNameVector.size();
  int SampleNum = GeneDataFrame[0].size();
  if (method == "spearman") {
    for (int i = 0; i < GeneNum; ++i) {
      GeneDataFrame[i] = GetRanks(GeneDataFrame[i]);
    }
  }

  // open output file
  std::ofstream OutFile(out_file_name, std::ios::out);
  if (!OutFile.good()) {
    std::cerr << "Error while opening " << out_file_name << ".\n";
    return -2;
  }

  // calc
  if (fdr) { // fdr
    // declaration of results variables
    std::vector <std::pair <int, int>> node_id_pairs;
    std::vector <double> corrs;
    std::vector <double> p_values;

    // multi-thread

    std::vector <std::thread> threads;
    for (int i = 0; i < thread_num; ++i) {
      threads.push_back(std::thread{ ThreadFunc_FDR, i, thread_num, GeneNum, SampleNum, std::ref(GeneDataFrame),
                                     std::ref(GeneNameVector), PvalCutoff, CorCutoff, signed_network, std::ref(node_id_pairs),
                                     std::ref(corrs), std::ref(p_values)});
    }

    for (auto & t : threads) {
      t.join();
    }

    // calc fdr
    std::vector <double> fdrs(p_values.size(), 0.0);
    double max_edges = GeneNum * (GeneNum - 1) * 0.5;
    calc_fdr(p_values, fdrs, max_edges);

    // output
    for (unsigned int i = 0; i <  p_values.size(); ++i) {
      std::string line = GeneNameVector[node_id_pairs[i].first] + '\t' + GeneNameVector[node_id_pairs[i].second]
                         + '\t' + std::to_string(corrs[i]) + '\t' + double_to_string(p_values[i]) + '\t'
                         + double_to_string(fdrs[i]) + '\n';
      OutFile << line;
    }
  } else { // non-fdr
    //multi-thread
    std::vector <std::thread> threads;
    for (int i = 0; i < thread_num; ++i) {
      threads.push_back(std::thread{ ThreadFunc, i, thread_num, GeneNum, SampleNum, std::ref(GeneDataFrame),
                                     std::ref(GeneNameVector), PvalCutoff, CorCutoff, signed_network, std::ref(OutFile)});
    }

    for (auto & t : threads) {
      t.join();
    }
  }

  return 0;
}


static std::mutex mutex_lock;


void ThreadFunc(int n, int thread_num, int GeneNum, int SampleNum, std::vector <std::vector <double>> & GeneDataFrame,
    std::vector <std::string> & GeneNameVector, double PvalCutoff, double CorCutoff, bool signed_network,
    std::ofstream & OutFile) {
  std::vector <std::string> tmp;

  // section I
  int start = n * GeneNum / (2 * thread_num);
  int stop = (n + 1) * GeneNum / (2 * thread_num);
  for (int i = start; i < stop; ++i) {
    for (int j = i + 1; j < GeneNum; ++j) {
      double corr = pearson(GeneDataFrame[i], GeneDataFrame[j]);
      if (signed_network and corr < 0) {
        continue;
      }
      double prob = get_p_value(corr, SampleNum);
      if (prob < PvalCutoff and std::abs(corr) > CorCutoff) {
        std::string item = GeneNameVector[i] + '\t' + GeneNameVector[j] + '\t' + std::to_string(corr) + '\t'
                           + double_to_string(prob) + '\n';
        tmp.push_back(item);
      }
      if (tmp.size() > 9999) {
        mutex_lock.lock();
        for (auto line : tmp) {
          OutFile << line;
        }
        mutex_lock.unlock();
        tmp.clear();
      }
    }
  }

  // section II
  start = (2 * thread_num - n - 1) * GeneNum / (2 * thread_num);
  stop = (2 * thread_num - n) * GeneNum / (2 * thread_num);
  for (int i = start; i < stop; ++i) {
    for (int j = i + 1; j < GeneNum; ++j) {
      double corr = pearson(GeneDataFrame[i], GeneDataFrame[j]);
      if (signed_network and corr < 0) {
        continue;
      }
      double prob = get_p_value(corr, SampleNum);
      if (prob < PvalCutoff and std::abs(corr) > CorCutoff) {
        std::string item = GeneNameVector[i] + '\t' + GeneNameVector[j] + '\t' + std::to_string(corr) + '\t'
                           + double_to_string(prob) + '\n';
        tmp.push_back(item);
      }
      if (tmp.size() > 9999) {
        mutex_lock.lock();
        for (auto line : tmp) {
          OutFile << line;
        }
        mutex_lock.unlock();
        tmp.clear();
      }
    }
  }

  if (tmp.size() > 0) {
    mutex_lock.lock();
    for (auto line : tmp) {
      OutFile << line;
    }
    mutex_lock.unlock();
    tmp.clear();
  }
}


void ThreadFunc_FDR(int n, int thread_num, int GeneNum, int SampleNum, std::vector <std::vector <double>> & GeneDataFrame,
    std::vector <std::string> & GeneNameVector, double PvalCutoff, double CorCutoff, bool signed_network,
    std::vector <std::pair <int, int>> & node_id_pairs, std::vector <double> & corrs, std::vector <double> & p_values) {
  std::vector <std::pair <int, int>> local_node_id_pairs;
  std::vector <double> local_corrs;
  std::vector <double> local_p_values;

  // section I
  int start = n * GeneNum / (2 * thread_num);
  int stop = (n + 1) * GeneNum / (2 * thread_num);
  for (int i = start; i < stop; ++i) {
    for (int j = i + 1; j < GeneNum; ++j) {
      double corr = pearson(GeneDataFrame[i], GeneDataFrame[j]);
      if (signed_network and corr < 0) {
        continue;
      }
      double prob = get_p_value(corr, SampleNum);
      if (prob < PvalCutoff and std::abs(corr) > CorCutoff) {
        local_node_id_pairs.push_back(std::pair <int, int> {i, j});
        local_corrs.push_back(corr);
        local_p_values.push_back(prob);
      }
    }
  }

  //section II
  start = (2 * thread_num - n - 1) * GeneNum / (2 * thread_num);
  stop = (2 * thread_num - n) * GeneNum / (2 * thread_num);
  for (int i = start; i < stop; ++i) {
    for (int j = i + 1; j < GeneNum; ++j) {
      double corr = pearson(GeneDataFrame[i], GeneDataFrame[j]);
      if (signed_network and corr < 0) {
        continue;
      }
      double prob = get_p_value(corr, SampleNum);
      if (prob < PvalCutoff and std::abs(corr) > CorCutoff) {
        local_node_id_pairs.push_back(std::pair <int, int> {i, j});
        local_corrs.push_back(corr);
        local_p_values.push_back(prob);
      }
    }
  }

  // local to global
  mutex_lock.lock();
  node_id_pairs.insert(node_id_pairs.end(), local_node_id_pairs.begin(), local_node_id_pairs.end());
  corrs.insert(corrs.end(), local_corrs.begin(), local_corrs.end());
  p_values.insert(p_values.end(), local_p_values.begin(), local_p_values.end());
  mutex_lock.unlock();
}
