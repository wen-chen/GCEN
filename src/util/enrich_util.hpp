#ifndef __GO_UTIL_H_
#define __GO_UTIL_H_


#include <cstdlib> // Needed to use the exit function
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include "func.hpp"
#include "strim.hpp"


class GO_term {
 public:
  std::string id;
  std::string name;
  std::string name_space;
  std::unordered_set <std::string> parents;
  std::unordered_set <std::string> children;
  int level; // shortest distance from root node
  int depth; // longest distance from root node
  bool is_obsolete = false;
  std::unordered_set <std::string> alt_ids;

  std::unordered_set <std::string> get_all_parents(std::unordered_map <std::string, GO_term> & go_term_map);
};


std::unordered_set <std::string> GO_term::get_all_parents(std::unordered_map <std::string, GO_term> & go_term_map) {
  std::unordered_set <std::string> all_parents;
  for (std::string parent: this -> parents) {
    all_parents.insert(parent);
    std::unordered_set <std::string> all_parents_of_parent = go_term_map[parent].get_all_parents(go_term_map);
    all_parents.insert(all_parents_of_parent.begin(), all_parents_of_parent.end());
  }
  return all_parents;
}


class GO_result {
 public:
  std::string id;
  std::string name;
  std::string name_space;
  char enrichment;
  int study_count;
  int study_n;
  int pop_count;
  int pop_n;
  double p_val;

  bool operator < (const GO_result & another) {
    int rank1, rank2;
    if (name_space == "biological_process") {
      rank1 = 1;
    } else if (name_space == "molecular_function") {
      rank1 = 2;
    } else if (name_space == "cellular_component") {
      rank1 = 3;
    } else {
      rank1 = 4;
    }
    if (another.name_space == "biological_process") {
      rank2 = 1;
    } else if (another.name_space == "molecular_function") {
      rank2 = 2;
    } else if (another.name_space == "cellular_component") {
      rank2 = 3;
    } else {
      rank2 = 4;
    }

    if (rank1 < rank2) {
      return true;
    } else if (rank1 == rank2) {
      if (p_val < another.p_val) {
        return true;
      }
    } 

    return false;
  }
  
};

class KO_result {
 public:
  std::string id;
  std::string name;
  char enrichment;
  int study_count;
  int study_n;
  int pop_count;
  int pop_n;
  double p_val;

  bool operator < (const KO_result & another) {
    if (p_val < another.p_val) {
      return true;
    }
    return false;
  }
};


void obo_parser(std::string & obo_file_name, std::unordered_map <std::string, GO_term> & go_term_map) {
  // open obo file
  std::ifstream obo_file(obo_file_name, std::ios::in);
  if (!obo_file.good()) {
    std::cerr << "Error while opening " << obo_file_name << ".\n";
    exit(-1);
  }

  // read obo file
  std::vector <GO_term> go_term_vector;
  bool flag = false;
  int i = -1;
  std::string lineString;

  while (getline(obo_file, lineString)) {
    strim(lineString);
    if (lineString[0] == '#') {
        continue;
      }
    if ((!flag) && (lineString.substr(0,6) == "[Term]")) {
      go_term_vector.push_back(GO_term());
      flag = true;
      i = i + 1;
    } else if (flag && (lineString.substr(0,4) == "id: ")) {
      go_term_vector[i].id = lineString.substr(4);
    } else if (flag && (lineString.substr(0,8) == "alt_id: ")) {
      go_term_vector[i].alt_ids.insert(lineString.substr(8));
    } else if (flag && (lineString.substr(0,6) == "name: ")) {
      go_term_vector[i].name = lineString.substr(6);
    } else if (flag && (lineString.substr(0,11) == "namespace: ")) {
      go_term_vector[i].name_space = lineString.substr(11);
    } else if (flag && (lineString.substr(0,6) == "is_a: ")) {
      go_term_vector[i].parents.insert(lineString.substr(6,10));
    } else if (flag && (lineString.substr(0,13) == "is_obsolete: ") && (lineString.substr(13) == "true")) {
      go_term_vector[i].is_obsolete = true;
    } else if (lineString.empty()) {
      flag = false;
    }
  }

  // vector to map
  for (GO_term go_term : go_term_vector) {
    go_term_map[go_term.id] = go_term;
  }
}


// go
void assoc_parser(std::string & assoc_file_name, std::unordered_map <std::string, std::unordered_set <std::string>> & assoc_map) {
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
    std::unordered_set <std::string> go_id_set;
    std::string go_id;
    while (getline(slineString, go_id, '\t')) {
      go_id_set.insert(go_id);
    }
    if (assoc_map.find(gene_name) != assoc_map.end()) {
      assoc_map[gene_name].insert(go_id_set.begin(), go_id_set.end());
    } else {
      assoc_map[gene_name] = go_id_set;
    }
  }
}


// kegg
void assoc_parser(std::string & assoc_file_name, 
    std::unordered_map <std::string, std::unordered_set<std::string>> & K_map,
    std::unordered_map <std::string, std::unordered_set <std::string>> & assoc_map) {
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
    std::string gene;
    getline(slineString, gene, '\t');
    std::unordered_set <std::string> ko_set;
    std::string K;
    while (getline(slineString, K, '\t')) {
      ko_set.insert(K_map[K].begin(), K_map[K].end());
    }
    if (assoc_map.find(gene) != assoc_map.end()) {
      assoc_map[gene].insert(ko_set.begin(), ko_set.end());
    } else {
      assoc_map[gene] = ko_set;
    }
  }
}


void propagate(std::unordered_map <std::string, std::unordered_set <std::string>> & assoc_map, std::unordered_map <std::string, GO_term> & go_term_map) {
  for (auto & assoc_item : assoc_map) {
    std::unordered_set <std::string> all_go_ids = assoc_item.second;
    for (auto go_id : assoc_item.second) {
      std::unordered_set <std::string> all_parents = go_term_map[go_id].get_all_parents(go_term_map);
      all_go_ids.insert(all_parents.begin(), all_parents.end());
    }
    assoc_item.second.insert(all_go_ids.begin(), all_go_ids.end());
  }
}


void load_network(std::string & network_file_name,
    std::unordered_map <std::string, std::unordered_set <std::string>> & assoc_map,
    std::unordered_map <std::string, std::unordered_set <std::string>> & network,
    std::unordered_set <std::string> & background_gene_set) {
  // open network file
  std::ifstream network_file(network_file_name, std::ios::in);
  if (!network_file.good()) {
    std::cerr << "Error while opening " << network_file_name << ".\n";
    exit(-1);
  }

  // read file
  std::string lineString;
  while (getline(network_file, lineString)) {
    strim(lineString);
    if (lineString[0] == '#') {
      continue;
    }
    std::vector <std::string> str_vec;
    split_string(lineString, str_vec, "\t");
    std::string geneA = str_vec[0];
    std::string geneB = str_vec[1];
    if ((assoc_map.find(geneA) == assoc_map.end()) && (assoc_map.find(geneB) != assoc_map.end())) {
      if (network.find(geneA) != network.end()) {
        network[geneA].insert(geneB);
      } else {
        network[geneA] = std::unordered_set < std::string > {};
        network[geneA].insert(geneB);
      }
      background_gene_set.insert(geneB);
    } else if ((assoc_map.find(geneA) != assoc_map.end()) && (assoc_map.find(geneB) == assoc_map.end())) {
      if (network.find(geneB) != network.end()) {
        network[geneB].insert(geneA);
      } else {
        network[geneB] = std::unordered_set <std::string> {};
        network[geneB].insert(geneA);
      }
      background_gene_set.insert(geneA);
    }
  }
}


void load_module(std::string & module_file_name, std::vector <std::unordered_set <std::string>> & module_vec, 
      std::unordered_set <std::string> & background_gene_set) {
  // open file
  std::ifstream module_file(module_file_name, std::ios::in);
  if (!module_file.good()) {
    std::cerr << "Error while opening " << module_file_name << ".\n";
    exit(-1);
  }

  // read file
  std::string lineString;
  while (getline(module_file, lineString)) {
    strim(lineString);
    if (lineString[0] == '#') {
      continue;
    }
    std::unordered_set <std::string> module;
    std::stringstream slineString;
    slineString << lineString;
    std::string singleString;
    while (getline(slineString, singleString, '\t')) {
      module.insert(singleString);
      background_gene_set.insert(singleString);
    }
    module_vec.push_back(module);
  }
}


void load_gene_list(std::string & gene_list_file_name,
        std::unordered_map <std::string, std::unordered_set <std::string>> & assoc_map,
        std::unordered_map <std::string, std::unordered_set <std::string>> & gene_map) {
  // open gene list file
  std::ifstream gene_list_file(gene_list_file_name, std::ios::in);
  if (!gene_list_file.good()) {
    std::cerr << "Error while opening " << gene_list_file_name << ".\n";
    exit(-1);
  }

  // read file
  std::string lineString;  
  while (getline(gene_list_file, lineString)) {
    strim(lineString);
    if (lineString[0] == '#') {
      continue;
    }
    if (assoc_map.find(lineString) != assoc_map.end()) {
      gene_map[lineString] = assoc_map[lineString];
    }
  }
}


void gene_set_to_map(std::unordered_set <std::string> & gene_set, 
  std::unordered_map <std::string, std::unordered_set <std::string>> & assoc_map,
  std::unordered_map <std::string, std::unordered_set <std::string>> & gene_map) {
  for (auto & gene : gene_set) {
    if (assoc_map.find(gene) != assoc_map.end()) {
      gene_map[gene] = assoc_map[gene];
    }
  }
}


int count(std::unordered_set <std::string> & gene_set,
    std::unordered_map <std::string, std::unordered_set <std::string>> & assoc_map,
    std::unordered_map <std::string, int> & go_count_map) {
  int n = 0;
  std::unordered_set <std::string> go_set;
  std::vector <std::string> go_vector;
  for (auto gene : gene_set) {
    if (assoc_map.find(gene) != assoc_map.end()) {
      for (auto go : assoc_map[gene]){
        go_set.insert(go);
        go_vector.push_back(go);
      }
    }
    n += 1;
  }

  for (auto go : go_set) {
    int counts = std::count(go_vector.begin(), go_vector.end(), go);
    go_count_map[go] = counts;
  }

  return n;
}


void count(std::unordered_map <std::string, std::unordered_set <std::string>> & gene_go_map, std::unordered_map <std::string, int> & go_count_map) {
  std::unordered_set <std::string> go_set;
  std::vector <std::string> go_vector;
  for (auto gene_go_item : gene_go_map) {
    for (auto go : gene_go_item.second) {
      go_set.insert(go);
      go_vector.push_back(go);
    }
  }

  for (auto go : go_set) {
    int counts = std::count(go_vector.begin(), go_vector.end(), go);
    go_count_map[go] = counts;
  }
}


void K2ko_parser(std::string & K2ko_file_name, 
    std::unordered_map <std::string, std::unordered_set<std::string>> & K_map, 
    std::unordered_map <std::string, std::string> & ko_map) {
  // open K2ko file
  std::ifstream K2ko_file(K2ko_file_name, std::ios::in);
  if (!K2ko_file.good()) {
    std::cerr << "Error while opening " << K2ko_file_name << ".\n";
    exit(-1);
  }

  // read K2ko file
  std::string lineString;
  while (getline(K2ko_file, lineString)) {
    strim(lineString);
    if (lineString[0] == '#') {
      continue;
    }
    std::vector <std::string> str_vec;
    split_string(lineString, str_vec, "\t");
    std::string K = str_vec[0];
    std::string ko = str_vec[1];
    std::string desc = str_vec[2];

    if (K_map.find(K) != K_map.end()) {
      K_map[K].insert(ko);
    } else {
      K_map[K] = std::unordered_set <std::string> {ko};
    }

    if (ko_map.find(ko) == ko_map.end()) {
      ko_map[ko] = desc;
    }
    
  }
}


#endif
