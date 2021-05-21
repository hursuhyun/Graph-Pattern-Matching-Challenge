/**
 * @file backtrack.cc
 *
 */

#include "backtrack.h"

Backtrack::Backtrack() {}
Backtrack::~Backtrack() {}

void Backtrack::PrintAllMatches(const Graph &data, const Graph &query,
                                const CandidateSet &cs) {
  std::ofstream out("output.txt");
  
  out<<"t "<<query.GetNumVertices()<<std::endl;


  int32_t query_size = query.GetNumVertices();
  size_t maxDegree = 0, degree;
  Vertex maxVertex = 0;
  Vertex v;
  bool is_answer = true;

  // find vertex of max degree -> inital vertex
  for(v=0; v<query_size; v++) {
    degree = query.GetDegree(v);
    if (degree == 0) {
      is_answer = false;
      break;
    }
    if (maxDegree < degree) {
      maxDegree = degree;
      maxVertex = v;
    }
  } v = maxVertex;
  
  // return if there can't be any answer
  if (is_answer == false) return;


  /**
   * finding out the order of matching
   * neighborSize : number of neighbors within matched vertices
   * matchingOrder : order of matching
   * extendableVertex : extendable vertices from matched ones
   * one : vertex with candidate size of one, if any
   */
  std::vector<Vertex> neighborSize, next, matchingOrder;
  std::set<Vertex>  extendableVertex;
  Vertex neighbor, one;
  int start, end, maximum = 0, deg = 0;

  neighborSize.resize(query_size, 0);
  matchingOrder.resize(query_size);

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


  /**
   * backtracking
   */

  std::vector<Vertex> answer;
  std::vector<size_t> triedCandidate;

  answer.resize(query_size, 0);
  triedCandidate.resize(query_size, 0);

  int j=0;
  Vertex current;
  Vertex currCandidate;

  while (j > -1) {

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
      //std::cout << "j: " << j << "\tis_answer: " << is_answer <<
      //"\t" << triedCandidate[current] << "\t" << cs.GetCandidateSize(current) << "\n";

      // if both condition 1, 2 are matched, current candidate is selected
      if (is_answer == true) {
        answer[matchingOrder[j]] = currCandidate;

        // if answer is completed, print the answer, else go to next round
        if (j == (query_size-1)) {
          out<<"a";
          for (int k=0; k<query_size; k++)
            out<<" "<<answer[k];
          out<<std::endl;
          triedCandidate[current] = 0;
          j--;
        } else {
          j++;
        }
        
      }

    } else {
      triedCandidate[current] = 0;
      j--;
    }
  } 
}
