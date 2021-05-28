/**
 * @file backtrack.h
 *
 */

#ifndef BACKTRACK_H_
#define BACKTRACK_H_

#include "candidate_set.h"
#include "common.h"
#include "graph.h"
#include <queue>
#include <vector>

class Backtrack {
 public:
  Backtrack();
  ~Backtrack();    
  
  void PrintAllMatches(std::string filename, const Graph &data, const Graph &query, const CandidateSet &cs);
  void bfsTraversal(const Graph &query, Vertex root_vertex, Tree *&tree, Vertex *&bfs_order);
};

#endif  // BACKTRACK_H_
