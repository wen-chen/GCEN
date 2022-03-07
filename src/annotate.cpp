#include <getopt.h>
#include <chrono>
#include <fstream>
#include <iostream>
#include <string>
#include <thread>  // Multithreading
#include <unordered_map>
#include <vector>
#include "third_party/ghc/filesystem.hpp"
#include "util/base.hpp"
#include "util/enrich.hpp"
#include "util/hypergeometric_p_value.hpp"
#include "util/thread_pool.hpp"

void annotate_help() {
  std::cout << version;
  std::cout << "annotate usage:\n";
  std::cout << "  annotate -g go-basic.obo -a gene_go_association_file -n "
               "input_network -o out_dir\n";
  std::cout << "options:\n";
  std::cout << "  -g --go <go-basic.obo file>\n";
  std::cout << "  -k --kegg <kegg information> (if the -g/--go is specified, "
               "the -k/--kegg are ignored)\n";
  std::cout << "  -a --assoc <gene-GO/KEGG association file>\n";
  std::cout << "  -n --network <network file>\n";
  std::cout << "  -m --module <module file> (if -n is specified, the -m is "
               "ignored)\n";
  std::cout << "  -p --pval <number> p value cutoff (default: 0.05)\n";
  std::cout << "  -o --output <output directory>\n";
  std::cout << "  -t --thread <number> cpu cores (default: 2)\n";
  std::cout << "  -v --version display GCEN version\n";
  std::cout << "  -h --help print help information\n";
  std::cout << "examples:\n";
  std::cout << "  ./annotate -g ../sample_data/go-basic.obo -a "
               "../sample_data/gene_go.assoc "
               "-n ../sample_data/gene_co_expr.network -o "
               "../sample_data/network_go_annotation\n";
  std::cout << "  ./annotate -g ../sample_data/go-basic.obo -a "
               "../sample_data/gene_go.assoc "
               "-m ../sample_data/module.txt -o "
               "../sample_data/module_go_annotation\n";
  std::cout << "  ./annotate -k ../sample_data/K2ko.tsv -a "
               "../sample_data/gene_kegg.assoc "
               "-n ../sample_data/gene_co_expr.network -o "
               "../sample_data/network_kegg_annotation\n";
  std::cout << "  ./annotate -k ../sample_data/K2ko.tsv -a "
               "../sample_data/gene_kegg.assoc "
               "-m ./sample_data/module.txt -o "
               "../sample_data/module_kegg_annotation\n";
}

static std::atomic<unsigned int> completed_num;

void network_go_annotate(
    std::string gene,
    std::unordered_map<std::string, std::unordered_set<std::string>>& network,
    std::unordered_map<std::string, std::unordered_set<std::string>>& assoc_map,
    std::unordered_map<std::string, GO_term>& go_term_map,
    std::unordered_map<std::string, int>& background_count_map, int pop_n,
    double pval_cutoff, std::string out_dir);

void module_go_annotate(
    int n, int thread_num,
    std::vector<std::unordered_set<std::string>>& module_vec,
    std::unordered_map<std::string, std::unordered_set<std::string>>& assoc_map,
    std::unordered_map<std::string, GO_term>& go_term_map,
    std::unordered_map<std::string, int>& background_count_map, int pop_n,
    double pval_cutoff, std::string& out_dir);

void network_kegg_annotate(
    int n, int thread_num, std::vector<std::string>& gene_vec,
    std::unordered_map<std::string, std::unordered_set<std::string>>& network,
    std::unordered_map<std::string, std::unordered_set<std::string>>& assoc_map,
    std::unordered_map<std::string, std::string> ko_map,
    std::unordered_map<std::string, int>& background_count_map, int pop_n,
    double pval_cutoff, std::string& out_dir);

void module_kegg_annotate(
    int n, int thread_num,
    std::vector<std::unordered_set<std::string>>& module_vec,
    std::unordered_map<std::string, std::unordered_set<std::string>>& assoc_map,
    std::unordered_map<std::string, std::string> ko_map,
    std::unordered_map<std::string, int>& background_count_map, int pop_n,
    double pval_cutoff, std::string& out_dir);

int main(int argc, char* argv[]) {
  if (argc < 2) {
    annotate_help();
    return 0;
  }

  // get options
  std::string obo_file_name = "";
  std::string K2ko_file_name = "";
  std::string assoc_file_name = "";
  std::string network_file_name = "";
  std::string module_file_name = "";
  std::string out_dir = "";
  double pval_cutoff = 0.05;
  int thread_num = 2;

  const char* const short_opts = "hvg:k:a:n:m:o:t:";
  const struct option long_opts[] = {
      {"help", 0, NULL, 'h'},   {"version", 0, NULL, 'v'},
      {"go", 1, NULL, 'g'},     {"kegg", 1, NULL, 'k'},
      {"assoc", 1, NULL, 'a'},  {"network", 1, NULL, 'n'},
      {"module", 1, NULL, 'm'}, {"output", 1, NULL, 'o'},
      {"thread", 1, NULL, 't'}, {"pval", 1, NULL, 'p'},
      {NULL, 0, NULL, 0}};
  int opt = getopt_long(argc, argv, short_opts, long_opts, NULL);
  while (opt != -1) {
    switch (opt) {
      case 'h':
        annotate_help();
        return 0;
      case 'v':
        display_version();
        return 0;
      case 'g':
        obo_file_name = optarg;
        break;
      case 'k':
        K2ko_file_name = optarg;
        break;
      case 'a':
        assoc_file_name = optarg;
        break;
      case 'n':
        network_file_name = optarg;
        break;
      case 'm':
        module_file_name = optarg;
        break;
      case 'o':
        out_dir = optarg;
        break;
      case 'p':
        pval_cutoff = std::stod(optarg);
        break;
      case 't':
        thread_num = std::stoi(optarg);
        break;
      case '?':
        annotate_help();
        return 0;
      case -1:
        break;
      default:
        return -1;
    }
    opt = getopt_long(argc, argv, short_opts, long_opts, NULL);
  }

  // check options
  if (obo_file_name.empty() and K2ko_file_name.empty()) {
    std::cerr << "Error: -g/--go or -k/--kegg is required but not specified!\n";
    exit(-1);
  }

  if (network_file_name.empty() and module_file_name.empty()) {
    std::cerr << "Error: -n/--network or -m/--module is required but not "
                 "specified!\n";
    exit(-1);
  }

  if (assoc_file_name.empty()) {
    std::cerr << "Error: -a/--assoc is required but not specified!\n";
    exit(-1);
  }

  if (out_dir.empty()) {
    std::cerr << "Error: -o --output is required but not specified!\n";
    exit(-1);
  }

  if (!obo_file_name.empty()) {  // go annotate
    // load obo
    std::unordered_map<std::string, GO_term> go_term_map;
    obo_parser(obo_file_name, go_term_map);

    // load assoc
    std::unordered_map<std::string, std::unordered_set<std::string>> assoc_map;
    assoc_parser(assoc_file_name, assoc_map);

    // propagate
    propagate(assoc_map, go_term_map);

    // mkdir
    ghc::filesystem::create_directories(out_dir);

    // annotation
    std::unordered_set<std::string> background_gene_set;

    if (!network_file_name.empty()) {
      // load network
      std::unordered_map<std::string, std::unordered_set<std::string>> network;
      load_network(network_file_name, assoc_map, network, background_gene_set);

      std::unordered_map<std::string, int> background_count_map;
      count(background_gene_set, assoc_map, background_count_map);
      int pop_n = background_gene_set.size();

      ThreadPool thead_pool(thread_num);

      unsigned int task_num = network.size();
      for (auto& network_item : network) {
        std::string gene = network_item.first;
        thead_pool.submit(new Task(network_go_annotate, gene, std::ref(network),
                                   std::ref(assoc_map), std::ref(go_term_map),
                                   std::ref(background_count_map), pop_n,
                                   pval_cutoff, out_dir));
      }

      double progress = 0.0;
      int barWidth = 70;
      while (true) {
        progress = (double)completed_num / task_num;
        std::cout << "[";
        int pos = barWidth * progress;
        for (int i = 0; i < barWidth; ++i) {
          if (i < pos)
            std::cout << "=";
          else if (i == pos)
            std::cout << ">";
          else
            std::cout << " ";
        }
        std::cout << "] " << int(progress * 100.0) << " %\r";
        std::cout.flush();
        if (!(progress < 1.0)) {
          break;
        }
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
      }
      std::cout << std::endl;

      thead_pool.wait();

    } else {
      // load module
      std::vector<std::unordered_set<std::string>> module_vec;
      load_module(module_file_name, module_vec, background_gene_set);

      std::unordered_map<std::string, int> background_count_map;
      count(background_gene_set, assoc_map, background_count_map);
      int pop_n = background_gene_set.size();

      std::vector<std::thread> threads;
      for (int i = 0; i < thread_num; ++i) {
        threads.push_back(std::thread{module_go_annotate, i, thread_num,
                                      std::ref(module_vec), std::ref(assoc_map),
                                      std::ref(go_term_map),
                                      std::ref(background_count_map), pop_n,
                                      pval_cutoff, std::ref(out_dir)});
      }

      for (auto& t : threads) {
        t.join();
      }
    }
  } else {  // kegg annotate
    // load K2ko
    std::unordered_map<std::string, std::unordered_set<std::string>> K_map;
    std::unordered_map<std::string, std::string> ko_map;
    K2ko_parser(K2ko_file_name, K_map, ko_map);

    // load assoc
    std::unordered_map<std::string, std::unordered_set<std::string>> assoc_map;
    assoc_parser(assoc_file_name, K_map, assoc_map);

    // mkdir
    ghc::filesystem::create_directories(out_dir);

    // annotate
    std::unordered_set<std::string> background_gene_set;

    if (!network_file_name.empty()) {
      // load network
      std::unordered_map<std::string, std::unordered_set<std::string>> network;

      load_network(network_file_name, assoc_map, network, background_gene_set);

      std::unordered_map<std::string, int> background_count_map;
      count(background_gene_set, assoc_map, background_count_map);
      int pop_n = background_gene_set.size();

      std::vector<std::string> gene_vec;
      for (auto network_item : network) {
        std::string gene = network_item.first;
        gene_vec.push_back(gene);
      }

      std::vector<std::thread> threads;
      for (int i = 0; i < thread_num; ++i) {
        threads.push_back(std::thread{network_kegg_annotate, i, thread_num,
                                      std::ref(gene_vec), std::ref(network),
                                      std::ref(assoc_map), std::ref(ko_map),
                                      std::ref(background_count_map), pop_n,
                                      pval_cutoff, std::ref(out_dir)});
      }

      for (auto& t : threads) {
        t.join();
      }

    } else {
      // load module
      std::vector<std::unordered_set<std::string>> module_vec;
      load_module(module_file_name, module_vec, background_gene_set);

      std::unordered_map<std::string, int> background_count_map;
      count(background_gene_set, assoc_map, background_count_map);
      int pop_n = background_gene_set.size();

      std::vector<std::thread> threads;
      for (int i = 0; i < thread_num; ++i) {
        threads.push_back(std::thread{module_kegg_annotate, i, thread_num,
                                      std::ref(module_vec), std::ref(assoc_map),
                                      std::ref(ko_map),
                                      std::ref(background_count_map), pop_n,
                                      pval_cutoff, std::ref(out_dir)});
      }

      for (auto& t : threads) {
        t.join();
      }
    }
  }

  return 0;
}

void network_go_annotate(
    std::string gene,
    std::unordered_map<std::string, std::unordered_set<std::string>>& network,
    std::unordered_map<std::string, std::unordered_set<std::string>>& assoc_map,
    std::unordered_map<std::string, GO_term>& go_term_map,
    std::unordered_map<std::string, int>& background_count_map, int pop_n,
    double pval_cutoff, std::string out_dir) {
  std::unordered_set<std::string> neighbor_gene_set = network[gene];

  std::unordered_map<std::string, int> enrichment_count_map;
  count(neighbor_gene_set, assoc_map, enrichment_count_map);
  int study_n = neighbor_gene_set.size();

  std::vector<GO_result> GO_result_vec;
  for (auto count_item : enrichment_count_map) {
    std::string go = count_item.first;
    int study_count = enrichment_count_map[go];
    int pop_count = background_count_map[go];
    double p_val =
        calc_p_hypergeometric(pop_n, pop_count, study_n, study_count);
    if (p_val < pval_cutoff) {
      GO_result result;
      result.id = go;
      result.name = go_term_map[go].name;
      result.name_space = go_term_map[go].name_space;
      if ((study_count / (double)study_n) > (pop_count / (double)pop_n)) {
        result.enrichment = 'e';
      } else {
        result.enrichment = 'p';
      }
      result.study_count = study_count;
      result.study_n = study_n;
      result.pop_count = pop_count;
      result.pop_n = pop_n;
      result.p_val = p_val;
      GO_result_vec.push_back(result);
    }
  }

  if (!GO_result_vec.empty()) {
    std::sort(GO_result_vec.begin(), GO_result_vec.end());
    std::ofstream result_file(out_dir + "/" + gene + ".go", std::ios::out);

    result_file << "GO\tname\tname_space\tenrichment\tstudy_count\t"
                   "study_n\tpop_count\tpop_n\tp_val\n";

    for (auto& result : GO_result_vec) {
      result_file << result.id << '\t' << result.name << '\t'
                  << result.name_space << '\t' << result.enrichment << '\t'
                  << result.study_count << '\t' << result.study_n << '\t'
                  << result.pop_count << '\t' << result.pop_n << '\t'
                  << result.p_val << '\n';
    }
    result_file.close();
  }
  completed_num += 1;
}

void module_go_annotate(
    int n, int thread_num,
    std::vector<std::unordered_set<std::string>>& module_vec,
    std::unordered_map<std::string, std::unordered_set<std::string>>& assoc_map,
    std::unordered_map<std::string, GO_term>& go_term_map,
    std::unordered_map<std::string, int>& background_count_map, int pop_n,
    double pval_cutoff, std::string& out_dir) {
  int module_num = module_vec.size();
  int start = n * module_num / thread_num;
  int stop = (n + 1) * module_num / thread_num;
  for (int i = start; i < stop; ++i) {
    std::unordered_set<std::string> module = module_vec[i];

    std::vector<GO_result> GO_result_vec;

    std::unordered_map<std::string, int> enrichment_count_map;
    int study_n = count(module, assoc_map, enrichment_count_map);

    for (auto count_item : enrichment_count_map) {
      std::string go = count_item.first;
      int study_count = enrichment_count_map[go];
      int pop_count = background_count_map[go];
      double p_val =
          calc_p_hypergeometric(pop_n, pop_count, study_n, study_count);
      if (p_val < pval_cutoff) {
        GO_result result;
        result.id = go;
        result.name = go_term_map[go].name;
        result.name_space = go_term_map[go].name_space;
        result.study_count = study_count;
        result.study_n = study_n;
        result.pop_count = pop_count;
        result.pop_n = pop_n;
        result.p_val = p_val;
        GO_result_vec.push_back(result);
      }
    }

    if (!GO_result_vec.empty()) {
      std::sort(GO_result_vec.begin(), GO_result_vec.end());

      std::ofstream module_file(
          out_dir + "/module_" + std::to_string(i) + ".txt", std::ios::out);
      for (auto& gene : module) {
        module_file << gene << '\n';
      }
      module_file.close();

      std::ofstream result_file(
          out_dir + "/module_" + std::to_string(i) + ".go", std::ios::out);

      result_file << "id\tname\tname_space\tstudy_count\tstudy_n\tpop_"
                     "count\tpop_n\tp_val\n";

      for (auto& result : GO_result_vec) {
        result_file << result.id << '\t' << result.name << '\t'
                    << result.name_space << '\t' << result.study_count << '\t'
                    << result.study_n << '\t' << result.pop_count << '\t'
                    << result.pop_n << '\t' << result.p_val << '\n';
      }

      result_file.close();
    }
  }
}

void network_kegg_annotate(
    int n, int thread_num, std::vector<std::string>& gene_vec,
    std::unordered_map<std::string, std::unordered_set<std::string>>& network,
    std::unordered_map<std::string, std::unordered_set<std::string>>& assoc_map,
    std::unordered_map<std::string, std::string> ko_map,
    std::unordered_map<std::string, int>& background_count_map, int pop_n,
    double pval_cutoff, std::string& out_dir) {
  int gene_num = gene_vec.size();
  int start = n * gene_num / thread_num;
  int stop = (n + 1) * gene_num / thread_num;

  for (int i = start; i < stop; ++i) {
    std::string gene = gene_vec[i];
    std::unordered_set<std::string> neighbor_gene_set = network[gene];

    std::unordered_map<std::string, int> enrichment_count_map;
    count(neighbor_gene_set, assoc_map, enrichment_count_map);
    int study_n = neighbor_gene_set.size();

    std::vector<KO_result> ko_result_vec;

    for (auto count_item : enrichment_count_map) {
      std::string ko = count_item.first;
      int study_count = enrichment_count_map[ko];
      int pop_count = background_count_map[ko];
      double p_val =
          calc_p_hypergeometric(pop_n, pop_count, study_n, study_count);
      if (p_val < pval_cutoff) {
        KO_result result;
        result.id = ko;
        result.name = ko_map[ko];
        if ((study_count / (double)study_n) > (pop_count / (double)pop_n)) {
          result.enrichment = 'e';
        } else {
          result.enrichment = 'p';
        }
        result.study_count = study_count;
        result.study_n = study_n;
        result.pop_count = pop_count;
        result.pop_n = pop_n;
        result.p_val = p_val;
        ko_result_vec.push_back(result);
      }
    }

    if (!ko_result_vec.empty()) {
      std::sort(ko_result_vec.begin(), ko_result_vec.end());

      std::ofstream result_file(out_dir + "/" + gene + ".kegg", std::ios::out);
      result_file << "ko\tname\tenrichment\tstudy_count\tstudy_n\tpop_count\t"
                     "pop_n\tp_val\n";
      for (auto& result : ko_result_vec) {
        result_file << result.id << '\t' << result.name << '\t'
                    << result.enrichment << '\t' << result.study_count << '\t'
                    << result.study_n << '\t' << result.pop_count << '\t'
                    << result.pop_n << '\t' << result.p_val << '\n';
      }

      result_file.close();
    }
  }
}

void module_kegg_annotate(
    int n, int thread_num,
    std::vector<std::unordered_set<std::string>>& module_vec,
    std::unordered_map<std::string, std::unordered_set<std::string>>& assoc_map,
    std::unordered_map<std::string, std::string> ko_map,
    std::unordered_map<std::string, int>& background_count_map, int pop_n,
    double pval_cutoff, std::string& out_dir) {
  int module_num = module_vec.size();
  int start = n * module_num / thread_num;
  int stop = (n + 1) * module_num / thread_num;
  for (int i = start; i < stop; ++i) {
    std::unordered_set<std::string> module = module_vec[i];

    std::vector<KO_result> ko_result_vec;

    std::unordered_map<std::string, int> enrichment_count_map;
    size_t study_n = count(module, assoc_map, enrichment_count_map);

    for (auto count_item : enrichment_count_map) {
      std::string ko = count_item.first;
      int study_count = enrichment_count_map[ko];
      int pop_count = background_count_map[ko];
      double p_val =
          calc_p_hypergeometric(pop_n, pop_count, study_n, study_count);
      if (p_val < pval_cutoff) {
        KO_result result;
        result.id = ko;
        result.name = ko_map[ko];
        if ((study_count / (double)study_n) > (pop_count / (double)pop_n)) {
          result.enrichment = 'e';
        } else {
          result.enrichment = 'p';
        }
        result.study_count = study_count;
        result.study_n = study_n;
        result.pop_count = pop_count;
        result.pop_n = pop_n;
        result.p_val = p_val;
        ko_result_vec.push_back(result);
      }
    }

    if (!ko_result_vec.empty()) {
      std::sort(ko_result_vec.begin(), ko_result_vec.end());

      std::ofstream module_file(
          out_dir + "/module_" + std::to_string(i) + ".txt", std::ios::out);
      for (auto& gene : module) {
        module_file << gene << '\n';
      }
      module_file.close();

      std::ofstream result_file(
          out_dir + "/module_" + std::to_string(i) + ".kegg", std::ios::out);

      result_file << "ko\tname\tenrichment\tstudy_count\tstudy_n\tpop_count\t"
                     "pop_n\tp_val\n";
      for (auto& result : ko_result_vec) {
        result_file << result.id << '\t' << result.name << '\t'
                    << result.enrichment << '\t' << result.study_count << '\t'
                    << result.study_n << '\t' << result.pop_count << '\t'
                    << result.pop_n << '\t' << result.p_val << '\n';
      }

      result_file.close();
    }
  }
}
