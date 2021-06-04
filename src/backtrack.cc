/**
 * @file backtrack.cc
 *
 */
#include "backtrack.h"
#include <stdio.h>
using namespace std;

Backtrack::Backtrack() {}
Backtrack::~Backtrack() {}

/**
 * @brief Backtrack to make DAG
 * save child, parent in dag
 */
void Backtrack::bfsTraversal(const Graph &query, Vertex root, DAG *&dag) {
    size_t vertex_num = query.GetNumVertices();

    queue<Vertex> bfs_queue;
    vector<bool> visited(vertex_num, false);

    dag = new DAG[vertex_num];
    for (size_t i = 0; i < vertex_num; ++i) {
        dag[i].init(vertex_num);
    }

    bfs_queue.push(root);
    visited[root] = true;

    while(!bfs_queue.empty()) {
        Vertex u = bfs_queue.front();
        bfs_queue.pop();

        //check neighbors, if not visited, push to queue
        size_t start = query.GetNeighborStartOffset(u);
        size_t end = query.GetNeighborEndOffset(u);
        std::vector<Vertex> neighborSize(vertex_num, 0);

        for (size_t t=start; t<end; t++){
              size_t u_nbr = query.GetNeighbor(t);
              if (!visited[u_nbr]) {
                bfs_queue.push(u_nbr);
                visited[u_nbr] = true;
                dag[u_nbr].parent_[dag[u_nbr].parent_count_++] = u;
                dag[u].children_[dag[u].children_count_++] = u_nbr;
            }
        }
    }
}

  /**
 * finding out the order of matching by candidate size
 * neighborSize : number of neighbors within matched vertices
 * matchingOrder : order of matching
 * extendableVertex : extendable vertices from matched ones
 */
std::vector<Vertex> Backtrack::FirstCSMin(const Graph &data, const Graph &query,
                                const CandidateSet &cs) { 
  size_t query_vertex_num = query.GetNumVertices();

  //select first node
  size_t min_degree = data.GetNumVertices();
  Vertex root = 0;

  for(size_t i=0; i<query_vertex_num; i++){
    if(min_degree > query.GetDegree(i)){
      min_degree = query.GetDegree(i);
      root = i;
    }
  }

  //build dag by bfs
  DAG *dag;
  bfsTraversal(query, root, dag);

  vector<bool> visited(query_vertex_num, false);
  queue<Vertex> extendable_queue;
  vector<size_t> extendable_vertex_order(query_vertex_num); //extendable_vertex_order[ith] = u  means u is ith number in extendable list
  size_t ex_index = 0 ;
  Vertex u = root;

  //calculate matching order using candidate size and extendable_vertex
  std::vector<Vertex> matchingOrder(query_vertex_num);
  for(size_t i = 0; i < query_vertex_num; i++){
    matchingOrder[i] = u;

    //find child_vertices of u that parents are all visited(extendable)
    //if extendable, check visited and put in the extendable_vertex_order vector
    //TODO: two graphs can be seperated
    for(size_t j = 0; j < dag[u].children_count_; j++){
      Vertex child = dag[u].children_[j];
      dag[child].visted_parent++;

      if(!visited[child] && dag[child].visted_parent == dag[child].parent_count_){
        extendable_vertex_order[ex_index++] =child;
        visited[child] = true;
      }
    }

    //if there is no child that is extendable, break
    //TODO: delete if possible
    if(ex_index == 0 && i == query_vertex_num -1) break;

    //find maximum candidate vertex in extendable_vertex_order
    //it would become next chosen matching order
    Vertex max_candidate = extendable_vertex_order[ex_index-1] ;
    size_t max_cs = cs.GetCandidateSize(max_candidate);
    size_t max_index =ex_index-1;
    for(size_t j = 0; j < ex_index-1; j++){
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

/**
 * Second matching order algorithm used.
 * 
 * return matchingOrder array
 */

std::vector<Vertex> Backtrack::SecondRIMin(const Graph &data, const Graph &query,
                                const CandidateSet &cs) { 

  int32_t query_size = query.GetNumVertices();


  // find vertex of minimum degree -> inital vertex

  size_t minDegree = data.GetNumVertices(), degree;
  Vertex minVertex = 0;
  Vertex v;

  for (v = 0; v < query_size; v++) {
    degree = query.GetDegree(v);
    if (minDegree > degree) {
      minDegree = degree;
      minVertex = v;
    }
  } v = minVertex;



  /**
   * finding out the order of matching
   * 
   * neighborSize : number of neighbors within matched vertices
   * matchingOrder : order of matching
   * reachableVertex : reachable vertices (those that have at least one edge with one of matched vertices)
   * one : vertex with candidate size of one, if any
   */

  std::vector<Vertex> neighborSize(query_size, 0), next, matchingOrder(query_size);
  std::set<Vertex>  reachableVertex;
  Vertex neighbor, one;
  int start, end, maximum = 0, deg = 0;


  for (int j=0; j<query_size; j++) {

    matchingOrder[j] = v;
    neighborSize[v] = -1;
    next.clear();
    reachableVertex.erase(v);
    one = query_size;
    

    //recomputing neighbor size

    start = query.GetNeighborStartOffset(v);
    end = query.GetNeighborEndOffset(v);
    maximum = 0;

    for (int t=start; t<end; t++){
          neighbor = query.GetNeighbor(t);
          if (neighborSize[neighbor] != -1) {
            reachableVertex.insert(neighbor);
            deg = ++neighborSize[neighbor];
          }
    }


    /**
     * searching for next vertex
     * 
     * Case 1: find vertex with candidate size one, if any
     * Case 2: if there is only one vertex with maximum neighbors, choose it
     * Case 3: else, choose one with minimum candidate size from those with maximum neighbors
     */

    for (std::set<Vertex>::iterator s=reachableVertex.begin(); s!=reachableVertex.end(); s++) {
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
     * 
     * Case 0: no reachable vertex, then find vertex with minimum degree again, within unchosen ones
     * Case 1: vertex with candidate size one exist, choose it
     * Case 2: if there is only one vertex with maximum neighbors, choose it
     * Case 3: else, choose one with minimum candidate size from those with maximum neighbors
     */
    if (reachableVertex.size() == 0) {
      size_t minDegree = data.GetNumVertices(), degree;
      Vertex minVertex = 0;
      Vertex v;

      for (v = 0; v < query_size; v++) {
        degree = query.GetDegree(v);
        if (minDegree > degree && neighborSize[neighbor] != -1) {
          minDegree = degree;
          minVertex = v;
        }
      } v = minVertex;
    } else if (one < query_size) {
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


/**
  * backtracking
  * 
  * parameter 'first' indicates whether this backtracking uses 
  * matching order of first algorithm 'FirstCSMin'
 */
bool Backtracking(std::vector<Vertex> matchingOrder, const Graph &data, const Graph &query,
                                const CandidateSet &cs, bool first) {

  int32_t query_size = query.GetNumVertices();
  std::vector<Vertex> answer(query_size, 0);
  std::vector<size_t> triedCandidate(query_size, 0);

  int j = 0, cnt = 0, j_cnt = 0;
  Vertex current, currCandidate;
  bool print_started = true, is_answer = true;

  /**
   *  perform backtracking until
   *  Case 1: j < -1 -> no answer left
   *  Case 2: cnt > 100000 -> finished printing 100000 answers
   */

  while (j > -1 && cnt < 100000) {

    // if too much time takes to print first answer, switch algorithm
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
  return print_started;
}

void Backtrack::PrintAllMatches(const Graph &data, const Graph &query,
                                const CandidateSet &cs) {
  std::cout << "t " << query.GetNumVertices() << "\n";

  // perform first algorithm
  std::vector<Vertex> matchingOrder = FirstCSMin(data, query, cs);
  
  // if algorithm should be switched, perform second algorithm
  if ( !Backtracking(matchingOrder, data, query, cs, true) ) {
    matchingOrder = SecondRIMin(data, query, cs);
    Backtracking(matchingOrder, data, query, cs, false);
  }
  return;
}
