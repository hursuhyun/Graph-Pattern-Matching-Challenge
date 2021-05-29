/**
 * @file backtrack.cc
 *
 */

#include "backtrack.h"
#include <stdio.h>
#define shell 0
using namespace std;
Backtrack::Backtrack() {}
Backtrack::~Backtrack() {}

//make tree using bfs
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

//how to get vertices by label..?
void PreprocessingLabel(const Graph &graph){
  size_t vertex_num = graph.GetNumVertices();
  for(size_t i=0; i<vertex_num; i++){

  }

}

void getVerticesByLabel(const Graph &graph, Label label){

}

/*
Vertex Backtrack::GetRootVertex(const Graph &data, const Graph &query) {
    size_t vertex_num = query.GetNumVertices();
    double min_score = query.GetNumVertices();
    Vertex start_vertex = 0;

    PreprocessingLabel(data);

    for (size_t i = 0; i < query.GetNumVertices(); ++i) {
        size_t degree = query.GetDegree(i);
        if (degree <= 1)
            continue;

        size_t count = 0;    
        
        Label label = query.GetLabel(i);
        size_t degree = query.GetDegree(i);
        count = 0;
        const std::vector<Vertex> data_vertices = data.getVerticesByLabel(label);

        for (Vertex v : data_vertices) {
            if (data.GetDegree(v) >= degree) {
                count += 1;
            }
        }

        double cur_score = count / (double)degree;
        if (cur_score < min_score) {
            min_score = cur_score;
            start_vertex = i;
        }
    }

    return start_vertex;
}
*/

std::vector<Vertex> FindingMatchingOrder(Vertex u, const Graph &data, const Graph &query,
                                const CandidateSet &c, Tree *&tree) {
  /**
   * finding out the order of matching
   * neighborSize : number of neighbors within matched vertices
   * matchingOrder : order of matching
   * extendableVertex : extendable vertices from matched ones
   * one : vertex with candidate size of one, if any
   */
  
  size_t query_vertex_num = query.GetNumVertices();

  vector<bool> visited(query_vertex_num, false);
  queue<Vertex> extendable_queue;
  vector<size_t> extendable_vertex_order(query_vertex_num); //extendable_vertex_order[ith] = u means u is ith number in extendable list
  size_t ex_index = 0 ;

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

    //TODO: when RI is same?
    if(ex_index == 0) break;
    Vertex max_ri_vertex = extendable_vertex_order[ex_index-1];
    size_t max_ri = tree[max_ri_vertex].parent_count_;
    size_t max_ri_index=ex_index-1;
    bool flag = true;
    /*
    if(c.GetCandidateSize(max_ri_vertex) ==0) {
      cout << "no candidates" << endl;
      exit(1);
    }
    if(c.GetCandidateSize(max_ri_vertex) ==1) {
      u = max_ri_vertex;
      flag = false;
    }*/

    for(int j=0; j<ex_index-1 && flag; j++){
      Vertex v = extendable_vertex_order[j];
      /*
      if(c.GetCandidateSize(v) ==0) {
        cout << "no candidates" << endl;
        exit(1);
      }
      if(c.GetCandidateSize(v) ==1) {
        u = max_ri_vertex;
        max_ri_vertex = v; 
        max_ri_index = j;
        flag = false;
      }
      */
      //if( max_ri < tree[v].parent_count_ || (max_ri == tree[v].parent_count_ && c.GetCandidateSize(max_ri_vertex) > c.GetCandidateSize(v)) ) {
      if(max_ri < tree[v].parent_count_ ) {
        max_ri = tree[v].parent_count_;
        max_ri_vertex = v; 
        max_ri_index = j;
      }
    }
    extendable_vertex_order[max_ri_index] = extendable_vertex_order[ex_index-1];
    ex_index--;

    u=max_ri_vertex;
  }


  return matchingOrder;

}


void Backtracking(std::string filename, std::vector<Vertex> matchingOrder, const Graph &data, const Graph &query,
                                const CandidateSet &cs) {
  std::ofstream out(filename);
  
  out<<"t "<<query.GetNumVertices()<<std::endl;

  bool is_answer = true; 
  std::vector<Vertex> answer;
  std::vector<size_t> triedCandidate;
  int32_t query_size = query.GetNumVertices();

  answer.resize(query_size, 0);
  triedCandidate.resize(query_size, 0);

  int j=0;
  Vertex current;
  Vertex currCandidate;
  int cnt=0;

  while (j > -1 && cnt < 100000) {

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
          cnt++;
          out<<"a";
          for (int k=0; k<query_size; k++)
            out<<" "<<answer[k];
          out<<std::endl;
          if(shell){
            std::cout << "a"; 
            for (int k=0; k<query_size; k++)
              std::cout<<" "<<answer[k];
            std::cout<<std::endl;
          }
        } else {
          j++;
        }
        
      }

    } else {
      triedCandidate[current] = 0;
      j--;
    }
  }
  std::cout << "count: " << cnt << endl;

  return;
}

void Backtrack::PrintAllMatches(std::string filename, const Graph &data, const Graph &query,
                                const CandidateSet &cs) {
  for(int i=0; i<7; i++){
    filename.pop_back();
  }
  filename.append("_output.txt");

  size_t query_vertex_num = query.GetNumVertices();

  //TODO: select first node
  //Vertex root = GetRootVertex();
  size_t max_degree = 0;
  Vertex root = 0;

  for(Vertex i=0; i<query_vertex_num; i++){
    if(max_degree < query.GetDegree(i)){
      max_degree = query.GetDegree(i);
      root = i;
    }
  }
  //bfs
  Tree *tree;
  Vertex *bfs_order;
  bfsTraversal(query, root, tree, bfs_order);
  //build DAG -> first in bfs_order is the parent

  //select extendable vertices
  //select with RI
  //finding out matching order
    std::vector<Vertex> matchingOrder = FindingMatchingOrder(root, data, query, cs, tree);
  if(shell){
    for(Vertex i=0; i<query_vertex_num; i++)
      std::cout << matchingOrder[i] << " ";
    std::cout << endl;
  }

  //backtracking and printing out
  Backtracking(filename, matchingOrder, data, query, cs);

}
