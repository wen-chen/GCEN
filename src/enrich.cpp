#include <getopt.h>
#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <algorithm> // sort result
#include "util/func.hpp"
#include "util/enrich_util.hpp"
#include "util/hypergeometric_p_value.hpp"


void enrich_help() {
  std::cout << "GCEN 0.5.0 by Wen Chen (chenwen@biochen.com, https://www.biochen.com/gcen)\n";
  std::cout << "enrich usage:\n";
  std::cout << "  enrich -e enrichment_gene_list_file -b background_gene_list_file -g go-basic.obo "
               "-a gene_go_association_file -p p_value_cutoff -o out_put_file\n";
  std::cout << "options:\n";
  std::cout << "  -e --enrich <enrichment gene list file>\n";
  std::cout << "  -b --background <background gene list file>\n";
  std::cout << "  -g --go <go-basic.obo file>\n";
  std::cout << "  -k --kegg <kegg information> (if the -g/--go is specified, the -k/--kegg are ignored)\n";
  std::cout << "  -a --assoc <gene/go association file>\n";
  std::cout << "  -p --pval <number> p value cutoff (default: 0.05)\n";
  std::cout << "  -o --output <output file>\n";
  std::cout << "  -v --version display GCEN version\n";
  std::cout << "  -h --help print help information\n";
  std::cout << "examples:\n";
  std::cout << "  enrich -e ../sample_data/enrichment_gene.list -b ../sample_data/background_gene.list "
               "-g ../sample_data/go-basic.obo -a ../sample_data/gene_go.assoc -p 0.05 -o ../sample_data/enrichment.go\n";
  std::cout << "  enrich -e ../sample_data/enrichment_gene.list -b ../sample_data/background_gene.list "
               "-k ../sample_data/K2ko.tsv -a ../sample_data/gene_kegg.assoc -p 0.05 -o ../sample_data/enrichment.kegg\n";
}


int main(int argc, char* argv[]) {
  if (argc < 2) {
    enrich_help();
    return -1;
  }

  // get options
  std::string enrichment_gene_file_name = "";
  std::string background_gene_file_name = "";
  std::string obo_file_name = "";
  std::string K2ko_file_name = "";
  std::string assoc_file_name = "";
  std::string out_file_name = "";
  double pval_cutoff = 0.5;

  const char * const short_opts = "hve:b:g:k:a:p:o:";
  const struct option long_opts[] =  {
    { "help", 0, NULL, 'h' },
    { "version", 0, NULL, 'v' },
    { "enrich", 1, NULL, 'e' },
    { "background", 1, NULL, 'b' },
    { "go", 1, NULL, 'g' },
    { "kegg", 1, NULL, 'k' },
    { "assoc", 1, NULL, 'a' },
    { "pval", 1, NULL, 'p' },
    { "output", 1, NULL, 'o'},
    { NULL, 0, NULL, 0 }
  };
  int opt = getopt_long(argc, argv, short_opts, long_opts, NULL);
  while (opt != -1) {
    switch (opt) {
      case 'h':
        enrich_help();
        return 0;
      case 'v':
        display_version();
        return 0;
      case 'e':
        enrichment_gene_file_name = optarg;
        break;
      case 'b':
        background_gene_file_name = optarg;
        break;
      case 'g':
        obo_file_name = optarg;
        break;
      case 'k':
        K2ko_file_name = optarg;
        break;
      case 'a':
        assoc_file_name = optarg;
        break;
      case 'p':
        pval_cutoff = std::stod(optarg);
        break;
      case 'o':
        out_file_name = optarg;
        break;
      case '?':
        enrich_help();
        return 0;
      case -1:
        break;
      default:
        return -1;
    }
    opt = getopt_long( argc, argv, short_opts, long_opts, NULL );
  }
  
  // check options
  if (obo_file_name.empty() and K2ko_file_name.empty()) {
    std::cerr << "Error: -g/--go is required but not specified!\n";
    return -1;
  }

  if (enrichment_gene_file_name.empty()) {
    std::cerr << "Error: -e/--enrich is required but not specified!\n";
    return -1;
  }

  if (background_gene_file_name.empty()) {
    std::cerr << "Error: -b/--background is required but not specified!\n";
    return -1;
  }

  if (assoc_file_name.empty()) {
    std::cerr << "Error: -a/--assoc is required but not specified!\n";
    return -1;
  }

  if (out_file_name.empty()) {
    std::cerr << "Error: -o --output is required but not specified!\n";
    return -1;
  }

  // enrich
  if (!obo_file_name.empty()) { // go enrich
    // load obo
    std::unordered_map <std::string, GO_term> go_term_map;
    obo_parser(obo_file_name, go_term_map);

    // load assoc
    std::unordered_map <std::string, std::unordered_set <std::string>> assoc_map;
    assoc_parser(assoc_file_name, assoc_map);

    // propagate
    propagate(assoc_map, go_term_map);

    // load gene list
    std::unordered_map <std::string, std::unordered_set <std::string>> enrichment_gene_go_map;
    load_gene_list(enrichment_gene_file_name, assoc_map, enrichment_gene_go_map);

    std::unordered_map <std::string, std::unordered_set <std::string>> background_gene_go_map;
    load_gene_list(background_gene_file_name, assoc_map, background_gene_go_map);

    //count go number
    std::unordered_map <std::string, int> enrichment_count_map;
    count(enrichment_gene_go_map, enrichment_count_map);

    std::unordered_map <std::string, int> background_count_map;
    count(background_gene_go_map, background_count_map);
    
    // enrich
    std::ofstream result_file(out_file_name, std::ios::out);

    int study_n = enrichment_gene_go_map.size();
    int pop_n = background_gene_go_map.size();
    std::vector <GO_result> GO_result_vec;

    for (auto count_item : enrichment_count_map) {
      std::string go = count_item.first;
      int study_count = enrichment_count_map[go];
      int pop_count = background_count_map[go];
      double p_val = calc_p_hypergeometric(pop_n, pop_count, study_n, study_count);
      if (p_val < pval_cutoff) {
        GO_result result;
        result.id = go;
        result.name = go_term_map[go].name;
        result.name_space = go_term_map[go].name_space;
        if ((study_count / (double) study_n) > (pop_count / (double) pop_n)) {
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

    std::sort(GO_result_vec.begin(), GO_result_vec.end());
  
    result_file << "GO\tname\tname_space\tenrichment\tstudy_count/study_n\tpop_count/pop_n\tp_val\n";
    for (auto & result : GO_result_vec) {
      result_file << result.id << '\t' << result.name << '\t' << result.name_space << '\t' <<
      result.enrichment << '\t' << result.study_count << '/' << result.study_n << '\t' << 
      result.pop_count << '/' << result.pop_n << '\t' << result.p_val << '\n';
    }

    result_file.close();
  } else { // kegg enrich
    // load K2ko
    std::unordered_map <std::string, std::unordered_set<std::string>> K_map;
    std::unordered_map <std::string, std::string> ko_map;
    K2ko_parser(K2ko_file_name, K_map, ko_map);

    // load assoc
    std::unordered_map <std::string, std::unordered_set <std::string>> assoc_map;
    assoc_parser(assoc_file_name, K_map, assoc_map);

    // load gene list
    std::unordered_map <std::string, std::unordered_set <std::string>> enrichment_gene_map;
    load_gene_list(enrichment_gene_file_name, assoc_map, enrichment_gene_map);

    std::unordered_map <std::string, std::unordered_set <std::string>> background_gene_map;
    load_gene_list(background_gene_file_name, assoc_map, background_gene_map);

    // count ko
    std::unordered_map <std::string, int> enrichment_count_map;
    count(enrichment_gene_map, enrichment_count_map);

    std::unordered_map <std::string, int> background_count_map;
    count(background_gene_map, background_count_map);
    
    // enrichment
    std::ofstream result_file(out_file_name, std::ios::out);

    int study_n = enrichment_gene_map.size();
    int pop_n = background_gene_map.size();
    std::vector <KO_result> ko_result_vec;

    for (auto count_item : enrichment_count_map) {
      std::string ko = count_item.first;
      int study_count = enrichment_count_map[ko];
      int pop_count = background_count_map[ko];
      double p_val = calc_p_hypergeometric(pop_n, pop_count, study_n, study_count);
      if (p_val < pval_cutoff) {
        KO_result result;
        result.id = ko;
        result.name = ko_map[ko];
        if ((study_count / (double) study_n) > (pop_count / (double) pop_n)) {
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

    std::sort(ko_result_vec.begin(), ko_result_vec.end());

    result_file << "ko\tname\tenrichment\tstudy_count/study_n\tpop_count/pop_n\tp_val\n";
    for (auto & result : ko_result_vec) {
      result_file << result.id << '\t' << result.name << '\t' <<
      result.enrichment << '\t' << result.study_count << '/' << result.study_n << '\t' << 
      result.pop_count << '/' << result.pop_n << '\t' << result.p_val << '\n';
    }

    result_file.close();
  }

  return 0;
}
