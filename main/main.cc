/**
 * @file main.cc
 *
 */

#include "backtrack.h"
#include "candidate_set.h"
#include "common.h"
#include "graph.h"



int main(int argc, char* argv[]) {
  if (argc < 1) {
    std::cerr << "Usage: ./program <file name base>\n";
    return EXIT_FAILURE;
  }

  std::string base_file_name = argv[1];

  std::string data_file_name = "", query_file_name = "", candidate_set_file_name = "";

  data_file_name.append(base_file_name);
  candidate_set_file_name.append(base_file_name);
  query_file_name.append(base_file_name);

  data_file_name.pop_back();
  data_file_name.pop_back();
  data_file_name.pop_back();

  candidate_set_file_name.append(".cs");
  query_file_name.append(".igraph");
  data_file_name.append(".igraph");

  Graph data(data_file_name);
  Graph query(query_file_name, true);
  CandidateSet candidate_set(candidate_set_file_name);

  Backtrack backtrack;

  backtrack.PrintAllMatches(query_file_name, data, query, candidate_set);

  return EXIT_SUCCESS;
}
