/**
 * @file backtrack.cc
 *
 */
#include "backtrack.h"
#include <stdio.h>
using namespace std;

Backtrack::Backtrack() {}
Backtrack::~Backtrack() {}



void Backtrack::bfsTraversal(const Graph &query, Vertex root_vertex, Tree *&tree, Vertex *&bfs_order) {
    size_t vertex_num = query.GetNumVertices();

    queue<Vertex> bfs_queue;
    vector<bool> visited(vertex_num, false);

    tree = new Tree[vertex_num];
    for (size_t i = 0; i < vertex_num; ++i) {
        tree[i].initialize(vertex_num);
    }
    bfs_order = new Vertex[vertex_num];

    size_t visited_vertex_count = 0;
    bfs_queue.push(root_vertex);
    visited[root_vertex] = true;
    tree[root_vertex].level_ = 0;

    while(!bfs_queue.empty()) {
        Vertex u = bfs_queue.front();
        bfs_queue.pop();
        bfs_order[visited_vertex_count++] = u;

        size_t u_nbrs_count;

        size_t start = query.GetNeighborStartOffset(u);
        size_t end = query.GetNeighborEndOffset(u);
        size_t maximum = 0;
        std::vector<Vertex> neighborSize(vertex_num, 0);

        for (size_t t=start; t<end; t++){
              size_t u_nbr = query.GetNeighbor(t);
              if (!visited[u_nbr]) {
                bfs_queue.push(u_nbr);
                visited[u_nbr] = true;
                tree[u_nbr].level_ = tree[u] .level_ + 1;
                tree[u_nbr].parent_[tree[u_nbr].parent_count_++] = u;
                tree[u].children_[tree[u].children_count_++] = u_nbr;
            }
        }
    }
}

std::vector<Vertex> Backtrack::FirstCSMin(const Graph &data, const Graph &query,
                                const CandidateSet &cs) { 
  size_t query_vertex_num = query.GetNumVertices();

  //TODO: select first node
  //Vertex root = GetRootVertex();
  size_t min_degree = data.GetNumVertices();
  Vertex root = 0;

  for(Vertex i=0; i<query_vertex_num; i++){
    if(min_degree > query.GetDegree(i)){
      min_degree = query.GetDegree(i);
      root = i;
    }
  }
  //bfs
  Tree *tree;
  Vertex *bfs_order;
  bfsTraversal(query, root, tree, bfs_order);

    /**
   * finding out the order of matching
   * neighborSize : number of neighbors within matched vertices
   * matchingOrder : order of matching
   * extendableVertex : extendable vertices from matched ones
   * one : vertex with candidate size of one, if any
   */
  
  vector<bool> visited(query_vertex_num, false);
  queue<Vertex> extendable_queue;
  vector<size_t> extendable_vertex_order(query_vertex_num); //extendable_vertex_order[ith] = u means u is ith number in extendable list
  size_t ex_index = 0 ;
  Vertex u = root;

  std::vector<Vertex> matchingOrder(query_vertex_num);
  for(size_t i = 0; i < query_vertex_num; i++){
    matchingOrder[i] = u;
    for(size_t j = 0; j < tree[u].children_count_; j++){
      Vertex child = tree[u].children_[j];
      tree[child].visted_parent++;

      if(!visited[child] && tree[child].visted_parent == tree[child].parent_count_){
        extendable_vertex_order[ex_index++] =child;
        visited[child] = true;
      }
    }

    if(ex_index == 0) break;
    Vertex max_candidate = extendable_vertex_order[ex_index-1] ;
    size_t max_cs = cs.GetCandidateSize(max_candidate);
    size_t max_index =ex_index-1;
    for(int j = 0; j < ex_index-1; j++){
      Vertex can = extendable_vertex_order[j];
      if(max_cs < cs.GetCandidateSize(can) ){
        max_candidate = can;
        max_cs = cs.GetCandidateSize(can);
        max_index = j;
      }
    }

    extendable_vertex_order[max_index] = extendable_vertex_order[ex_index-1] ;
    ex_index --;
    u=max_candidate;
  }

  return matchingOrder;
}

std::vector<Vertex> Backtrack::SecondRIMin(const Graph &data, const Graph &query,
                                const CandidateSet &cs) { 
  int32_t query_size = query.GetNumVertices();
  size_t maxDegree = data.GetNumVertices(), degree, cs_size;
  Vertex maxVertex = 0;
  Vertex v = 0;
  bool is_answer = true;

  // find vertex of max degree -> inital vertex
  for(v=0; v<query_size; v++) {
    degree = query.GetDegree(v);
    if (maxDegree > degree) {
      maxDegree = degree;
      maxVertex = v;
    }
  } v = 0;

  /**
   * finding out the order of matching
   * neighborSize : number of neighbors within matched vertices
   * matchingOrder : order of matching
   * extendableVertex : extendable vertices from matched ones
   * one : vertex with candidate size of one, if any
   */
  std::vector<Vertex> neighborSize(query_size, 0), next, matchingOrder(query_size);
  std::set<Vertex>  extendableVertex;
  Vertex neighbor, one;
  int start, end, maximum = 0, deg = 0;


  for (int j=0; j<query_size; j++) {

    matchingOrder[j] = v;
    neighborSize[v] = -1;
    next.clear();
    extendableVertex.erase(v);
    one = query_size;
    
    start = query.GetNeighborStartOffset(v);
    end = query.GetNeighborEndOffset(v);
    maximum = 0;

    //recomputing neighbor size
    for (int t=start; t<end; t++){
          neighbor = query.GetNeighbor(t);
          if (neighborSize[neighbor] != -1) {
            extendableVertex.insert(neighbor);
            deg = ++neighborSize[neighbor];
          }
    }

    /**
     * searching for next vertex
     */
    for (std::set<Vertex>::iterator s=extendableVertex.begin(); s!=extendableVertex.end(); s++) {
      if (cs.GetCandidateSize(*s) == 1) {
        if (*s < one) one = *s;
      } else {
        deg = neighborSize[*s];
        if (maximum < deg) {
          maximum = deg;
          next.clear();
          next.push_back(*s);
        } else if (maximum == deg) {
          next.push_back(*s);
        }
      }
    }

    /**
     * renewing current vertex
     * if there is vertex with candidate size one, choose it
     * else if there is only one vertex with maximum neighbors, choose it
     * else choose one with minimum candidate size from those with maximum neighbors
     */
    if (one < query_size) {
      v = one;
    } else if (next.size() == 1) {
      v = next.front();
    } else {
      size_t cs_min = data.GetNumVertices();
      while (next.size() > 0) {
        Vertex val = next.back();
        size_t cs_size = cs.GetCandidateSize(val);
        if (cs_size < cs_min) {
          v = val;
          cs_min = cs_size;
        } else if (cs_min == cs_size && val < v) v = val;
        next.pop_back();
      }
    } 

  }

  return matchingOrder;
}


bool Backtracking(std::vector<Vertex> matchingOrder, const Graph &data, const Graph &query,
                                const CandidateSet &cs, bool first) {


  /**
   * backtracking
   */

  bool is_answer = true; 
  int32_t query_size = query.GetNumVertices();
  std::vector<Vertex> answer(query_size, 0);
  std::vector<size_t> triedCandidate(query_size, 0);
  std::vector<bool> occupiedCandidate(data.GetNumVertices(), false);

  int j = 0;
  Vertex current;
  Vertex currCandidate;
  int cnt = 0;
  int j_cnt = 0;
  bool print_started = true;

  while (j > -1 && cnt < 100000) {
    j_cnt++;
    if(first && j_cnt > 1000000 && cnt==0) {
      print_started = false;
      break;
    }
    current = matchingOrder[j];

    // check if there is untried candidate of current vertex in query
    if (triedCandidate[current] < cs.GetCandidateSize(current)) {
      
      currCandidate = cs.GetCandidate(current, triedCandidate[current]);
      triedCandidate[current]++;
      is_answer = true;

      for (int i=0; i<j; i++) {

        // condition 1: check if current candidate is already included in the answer
        if (answer[matchingOrder[i]] == currCandidate){
          is_answer = false;
          break;
        }

        // condition 2: check if every edge in query have corresponding edge in data graph
        if (query.IsNeighbor(matchingOrder[i], current)) {
          if (!data.IsNeighbor(answer[matchingOrder[i]], currCandidate)) {
            is_answer = false;
            break;
          }
        }
      }

      // if both condition 1, 2 are matched, current candidate is selected
      if (is_answer == true) {
        answer[matchingOrder[j]] = currCandidate;

        // if answer is completed, print the answer, else go to next round
        if (j == (query_size-1)) {
          cnt++;
          std::printf("a");
          for (int k=0; k<query_size; k++)
            std::printf(" %d", answer[k]);
          std::printf("\n");
        } else {
          j++;
        }
        
      }

    } else {
      triedCandidate[current] = 0;
      j--;
    }
  } 
  std::printf("%d\n", cnt);
  return print_started;
}

void Backtrack::PrintAllMatches(const Graph &data, const Graph &query,
                                const CandidateSet &cs) {
  std::cout << "t " << query.GetNumVertices() << "\n";

  std::vector<Vertex> matchingOrder = FirstCSMin(data, query, cs);
  
  if ( !Backtracking(matchingOrder, data, query, cs, true) ) {
    std::cout << "matching algorithm switched"<<std::endl;
    matchingOrder = SecondRIMin(data, query, cs);
    Backtracking(matchingOrder, data, query, cs, false);
  }
}
