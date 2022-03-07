#include <getopt.h>
#include <chrono>
#include <fstream>
#include <iostream>
#include <mutex>
#include <string>
#include <thread>  // Multithreading
#include <vector>
#include "util/base.hpp"
#include "util/enrich.hpp"
#include "util/hypergeometric_p_value.hpp"
#include "util/thread_pool.hpp"

void calculate_accuracy_help() {
  std::cout << version;
  std::cout << "calculate_accuracy usage:\n";
  std::cout << "  calculate_accuracy -g go-basic.obo -a "
               "gene_go_association_file -n input_network\n";
  std::cout << "options:\n";
  std::cout << "  -g --go <go-basic.obo file>\n";
  std::cout << "  -a --assoc <gene-GO/KEGG association file>\n";
  std::cout << "  -n --network <network file>\n";
  std::cout << "  -p --pval <number> p value cutoff (default: 0.05)\n";
  std::cout << "  -t --thread <number> cpu cores (default: 2)\n";
  std::cout << "  -v --version display GCEN version\n";
  std::cout << "  -h --help print help information\n";
  std::cout << "examples:\n";
  std::cout << "  ./calculate_accuracy -g ../sample_data/go-basic.obo -a "
               "../sample_data/gene_go.assoc -n "
               "../sample_data/gene_co_expr.network\n";
}

void obo_parser_2(std::string& obo_file_name,
                  std::unordered_map<std::string, GO_term>&
                      go_term_map) {  // for calculate accuracy
  // open obo file
  std::ifstream obo_file(obo_file_name, std::ios::in);
  if (!obo_file.good()) {
    std::cerr << "Error while opening " << obo_file_name << ".\n";
    exit(-1);
  }

  // read obo file
  std::vector<GO_term> go_term_vector;
  bool flag = false;
  int i = -1;
  std::string line;

  while (getline(obo_file, line)) {
    strim(line);
    if (line[0] == '#') {
      continue;
    }
    if ((!flag) && (line.substr(0, 6) == "[Term]")) {
      go_term_vector.push_back(GO_term());
      flag = true;
      i = i + 1;
    } else if (flag && (line.substr(0, 4) == "id: ")) {
      go_term_vector[i].id = line.substr(4);
    } else if (flag && (line.substr(0, 8) == "alt_id: ")) {
      go_term_vector[i].alt_ids.insert(line.substr(8));
    } else if (flag && (line.substr(0, 6) == "name: ")) {
      go_term_vector[i].name = line.substr(6);
    } else if (flag && (line.substr(0, 11) == "namespace: ")) {
      go_term_vector[i].name_space = line.substr(11);
    } else if (flag && (line.substr(0, 6) == "is_a: ")) {
      go_term_vector[i].parents.insert(line.substr(6, 10));
    } else if (flag && (line.substr(0, 13) == "is_obsolete: ") &&
               (line.substr(13) == "true")) {
      go_term_vector[i].is_obsolete = true;
    } else if (line.empty()) {
      flag = false;
    }
  }

  // vector to map
  for (GO_term go_term : go_term_vector) {
    if (go_term.name_space == "biological_process") {
      go_term_map[go_term.id] = go_term;
    }
  }
}

void assoc_parser_2(
    std::string& assoc_file_name,
    std::unordered_map<std::string, GO_term>& go_term_map,
    std::unordered_map<std::string, std::unordered_set<std::string>>&
        assoc_map) {  // for calculate accuracy
  // open assoc file
  std::ifstream assoc_file(assoc_file_name, std::ios::in);
  if (!assoc_file.good()) {
    std::cerr << "Error while opening " << assoc_file_name << ".\n";
    exit(-1);
  }

  // read assoc file
  std::string lineString;
  while (getline(assoc_file, lineString)) {
    strim(lineString);
    if (lineString[0] == '#') {
      continue;
    }
    std::stringstream slineString;
    slineString << lineString;
    std::string gene_name;
    getline(slineString, gene_name, '\t');
    std::unordered_set<std::string> go_id_set;
    std::string go_id;
    while (getline(slineString, go_id, '\t')) {
      if (go_term_map.find(go_id) != go_term_map.end()) {
        go_id_set.insert(go_id);
      }
    }
    if (assoc_map.find(gene_name) != assoc_map.end()) {
      assoc_map[gene_name].insert(go_id_set.begin(), go_id_set.end());
    } else {
      assoc_map[gene_name] = go_id_set;
    }
  }
}

void load_network_2(
    std::string& network_file_name,
    std::unordered_map<std::string, std::unordered_set<std::string>>& assoc_map,
    std::unordered_map<std::string, std::unordered_set<std::string>>& network,
    std::unordered_set<std::string>&
        background_gene_set) {  // for calculate accuracy
  // open network file
  std::ifstream network_file(network_file_name, std::ios::in);
  if (!network_file.good()) {
    std::cerr << "Error while opening " << network_file_name << ".\n";
    exit(-1);
  }

  // read file
  std::string line;
  while (getline(network_file, line)) {
    strim(line);
    if (line[0] == '#') {
      continue;
    }
    std::vector<std::string> str_vec;
    split_string(line, str_vec, "\t");
    std::string gene_a = str_vec[0];
    std::string gene_b = str_vec[1];
    if ((assoc_map.find(gene_a) != assoc_map.end()) &&
        (assoc_map.find(gene_b) != assoc_map.end())) {
      if (network.find(gene_a) != network.end()) {
        network[gene_a].insert(gene_b);
      } else {
        network[gene_a] = std::unordered_set<std::string>{};
        network[gene_a].insert(gene_b);
      }
      if (network.find(gene_b) != network.end()) {
        network[gene_b].insert(gene_a);
      } else {
        network[gene_b] = std::unordered_set<std::string>{};
        network[gene_b].insert(gene_a);
      }
      background_gene_set.insert(gene_a);
      background_gene_set.insert(gene_b);
    }
  }
}

static std::mutex mutex_lock;
static std::atomic<unsigned int> completed_num;

void network_go_annotate_2(  // for calculate accuracy
    std::string gene,
    std::unordered_map<std::string, std::unordered_set<std::string>>& network,
    std::unordered_map<std::string, std::unordered_set<std::string>>& assoc_map,
    std::unordered_map<std::string, GO_term>& go_term_map,
    std::unordered_map<std::string, int>& background_count_map, int pop_n,
    double pval_cutoff,
    std::unordered_map<std::string, std::unordered_set<std::string>>&
        results_map) {
  std::unordered_set<std::string> neighbor_gene_set = network[gene];
  int study_n = neighbor_gene_set.size();

  std::unordered_map<std::string, int> enrichment_count_map;
  count(neighbor_gene_set, assoc_map, enrichment_count_map);

  std::unordered_set<std::string> go_set;
  for (auto count_item : enrichment_count_map) {
    std::string go = count_item.first;
    int study_count = enrichment_count_map[go];
    int pop_count = background_count_map[go];
    double p_val =
        calc_p_hypergeometric(pop_n, pop_count, study_n, study_count);
    if ((p_val < pval_cutoff) &&
        ((study_count / (double)study_n) > (pop_count / (double)pop_n))) {
      go_set.insert(go);
    }
  }

  if (go_set.size() > 0) {
    mutex_lock.lock();
    results_map[gene] = go_set;
    mutex_lock.unlock();
  }
  completed_num += 1;
}

int main(int argc, char* argv[]) {
  if (argc < 2) {
    calculate_accuracy_help();
    return 0;
  }

  // get options
  std::string obo_file_name = "";
  std::string assoc_file_name = "";
  std::string network_file_name = "";
  double pval_cutoff = 0.05;
  int thread_num = 2;

  const char* const short_opts = "hvg:a:n:p:t:";
  const struct option long_opts[] = {
      {"help", 0, NULL, 'h'},   {"version", 0, NULL, 'v'},
      {"go", 1, NULL, 'g'},     {"assoc", 1, NULL, 'a'},
      {"pval", 1, NULL, 'p'},   {"network", 1, NULL, 'n'},

      {"thread", 1, NULL, 't'}, {NULL, 0, NULL, 0}};
  int opt = getopt_long(argc, argv, short_opts, long_opts, NULL);
  while (opt != -1) {
    switch (opt) {
      case 'h':
        calculate_accuracy_help();
        return 0;
      case 'v':
        display_version();
        return 0;
      case 'g':
        obo_file_name = optarg;
        break;
      case 'a':
        assoc_file_name = optarg;
        break;
      case 'n':
        network_file_name = optarg;
        break;
      case 'p':
        pval_cutoff = std::stod(optarg);
        break;
      case 't':
        thread_num = std::stoi(optarg);
        break;
      case '?':
        calculate_accuracy_help();
        return 0;
      case -1:
        break;
      default:
        return -1;
    }
    opt = getopt_long(argc, argv, short_opts, long_opts, NULL);
  }

  // check options
  if (obo_file_name.empty()) {
    std::cerr << "Error: -g/--go is required but not specified!\n";
    exit(-1);
  }

  if (network_file_name.empty()) {
    std::cerr << "Error: -n/--network is required but not "
                 "specified!\n";
    exit(-1);
  }

  if (assoc_file_name.empty()) {
    std::cerr << "Error: -a/--assoc is required but not specified!\n";
    exit(-1);
  }

  // load obo
  std::unordered_map<std::string, GO_term> go_term_map;
  obo_parser(obo_file_name, go_term_map);

  // load assoc
  std::unordered_map<std::string, std::unordered_set<std::string>> assoc_map;
  assoc_parser_2(assoc_file_name, go_term_map, assoc_map);

  // propagate
  std::unordered_map<std::string, std::unordered_set<std::string>>
      assoc_map_raw = assoc_map;
  propagate(assoc_map, go_term_map);

  // load network
  std::unordered_set<std::string> background_gene_set;
  std::unordered_map<std::string, std::unordered_set<std::string>> network;
  load_network_2(network_file_name, assoc_map, network, background_gene_set);

  // count background
  std::unordered_map<std::string, int> background_count_map;
  count(background_gene_set, assoc_map, background_count_map);
  int pop_n = background_gene_set.size();

  // go enrichment
  std::unordered_map<std::string, std::unordered_set<std::string>> results_map;

  ThreadPool thead_pool(thread_num);

  unsigned int task_num = network.size();
  for (auto& network_item : network) {
    std::string gene = network_item.first;
    thead_pool.submit(new Task(network_go_annotate_2, gene, std::ref(network),
                               std::ref(assoc_map), std::ref(go_term_map),
                               std::ref(background_count_map), pop_n,
                               pval_cutoff, std::ref(results_map)));
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

  // calculate accuracy
  int n = 0;
  for (auto& results : results_map) {
    std::string gene = results.first;
    std::unordered_set<std::string> go_set = results.second;
    for (auto& go : assoc_map_raw[gene]) {
      if (go_set.find(go) != go_set.end()) {
        n += 1;
        break;
      }
    }
  }
  std::cout << "Prediction accuracy: " << ((double)n) / results_map.size()
            << "\n";

  return 0;
}
