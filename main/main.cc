/**
 * @file main.cc
 *
 */

#include "backtrack.h"
#include "candidate_set.h"
#include "common.h"
#include "graph.h"

int main(int argc, char* argv[]) {
  if (argc < 4) {
    std::cerr << "Usage: ./program <data graph file> <query graph file> "
                 "<candidate set file>\n";
    return EXIT_FAILURE;
  }

  std::string data_file_name = argv[1];
  std::string query_file_name = argv[2];
  std::string candidate_set_file_name = argv[3];

  Graph data(data_file_name);
  Graph query(query_file_name, true);
  CandidateSet candidate_set(candidate_set_file_name);

  Backtrack backtrack;
  
   time_t start = time(NULL);
  backtrack.PrintAllMatches(query_file_name, data, query, candidate_set);
   time_t end = time(NULL);

  double time = (double)(end - start);
  std::cout << time << std::endl;
  return EXIT_SUCCESS;
}
