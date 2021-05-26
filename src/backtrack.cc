/**
 * @file backtrack.cc
 *
 */

#include "backtrack.h"
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

void Backtrack::PrintAllMatches(const Graph &data, const Graph &query,
                                const CandidateSet &cs) {
  std::cout << "t " << query.GetNumVertices() << "\n";

  size_t query_vertex_num = query.GetNumVertices();

  //TODO: select first node
  //Vertex root = GetRootVertex();
  size_t max_degree = 0;
  Vertex root = 0;

  for(Vertex i=0; i<query_vertex_num; i++){
    if(max_degree > query.GetDegree(i)){
      max_degree = query.GetDegree(i);
      root = i;
    }
  }
  //bfs
  Tree *tree;
  Vertex *bfs_order;
  bfsTraversal(query, root, tree, bfs_order);

  size_t cnt_vertex = 0;
  size_t ans_count = 0;
  size_t max_num = 100000;

  std::vector<Vertex> answer(query_vertex_num);
  //build DAG -> first in bfs_order is the parent

  //select extendable vertices
  //select with RI
  vector<bool> visited(query_vertex_num, false);
  queue<Vertex> extendable_queue;
  vector<Vertex> matchingOrder(query_vertex_num);
  vector<Vertex> extendable_vertex(query_vertex_num);       //extendable_vertex[u] = 1 if u is extendable
  vector<size_t> extendable_vertex_order(query_vertex_num); //extendable_vertex_order[ith] = u means u is ith number in extendable list
  size_t ex_index = 0 ;

  Vertex u = root;
  int b = 0;
  bool is_answer = true;

  while(ans_count < max_num &&  b > -1){
    //printf("s%d %d\n", b, u);
    for(size_t i = 0; i < cs.GetCandidateSize(u); i++){
      Vertex candidate = cs.GetCandidate(u, i);
      printf("%d %d %d\n", b, u, candidate);

      //if not included in the answer, and there is an edge(u-v, data u-data v)
      for (size_t k=0; k<b; k++) {
        // condition 1: check if current candidate is already included in the answer
        if (answer[matchingOrder[k]] == candidate){
          is_answer = false;
          break;
        }

        // condition 2: check if every edge in query have corresponding edge in data graph
        if (query.IsNeighbor(matchingOrder[k], u)) {
          if (!data.IsNeighbor(answer[matchingOrder[k]], candidate)) {
            is_answer = false;
            break;
          }
        }
      }
      printf("%d %d\n", b, u);

      if(!is_answer) continue;
      printf("%d %d\n", b, u);
      matchingOrder[b] = u;
      answer[u] = candidate;
      visited[u] = true;
      b++;

      if(b == query_vertex_num) {
        ans_count++;
        printf("a ");
        for(size_t k=0; k<query_vertex_num; k++)
          printf("%d ", answer[k]);
        printf("\n");
        b--;
        continue;
      }

      //find extendable vertices
      for(size_t j = 0; j < tree[u].children_count_; j++){
        Vertex child = tree[u].children_[j];
        printf("child %d\n", child);
        tree[child].visted_parent++;

        if(!visited[child] && tree[child].visted_parent == tree[child].parent_count_){
          //extendable_queue.push(child);
          extendable_vertex_order[ex_index++] =child;
          extendable_vertex[child] = 1;
          visited[child] = true;
        }
      }

      //TODO: when RI is same?
      Vertex max_ri_vertex = extendable_vertex_order[ex_index-1];
      size_t max_ri = tree[max_ri_vertex].parent_count_;
      size_t max_ri_index=ex_index-1;
      printf("max:%d max_v:%d\n", max_ri, max_ri_vertex);

      for(int j=0; j<ex_index-1; j++){
        Vertex v = extendable_vertex_order[j];
        printf("v:%d parent_count:%d\n",v, tree[v].parent_count_);
        if(max_ri < tree[v].parent_count_) {
          max_ri = tree[v].parent_count_;
          max_ri_vertex = v; 
          max_ri_index = j;
        }
      }
      u = max_ri_vertex;
      
      extendable_vertex_order[max_ri_index] = 


      /* using queue
      //no answer
      if(extendable_queue.empty()) exit(0);
      Vertex max_ri_vertex = extendable_queue.front();
      extendable_queue.pop();
      size_t max_ri = tree[max_ri_vertex].parent_count_;
      printf("max:%d max_v:%d\n", max_ri, max_ri_vertex);

      for(int j=0; j<extendable_queue.size(); j++){
        Vertex v = extendable_queue.front();
        printf("v:%d parent_count:%d\n",v, tree[v].parent_count_);
        extendable_queue.pop();
        if(max_ri < tree[v].parent_count_) {
          extendable_queue.push(max_ri_vertex);
          max_ri = tree[v].parent_count_;
          max_ri_vertex = v; 
        }else{
          extendable_queue.push(v);
        }
      }
      */

      printf("b%d u%d\n", b, u);
    }
    b--;
    for(size_t j = 0; j < tree[u].children_count_; j++){
      Vertex child = tree[u].children_[j];
      printf("child %d\n", child);
      tree[child].visted_parent--;

      if(!visited[child] && tree[child].visted_parent == tree[child].parent_count_){
        extendable_queue.push(child);
        //extendable_vertex[ex_index++] = child;
        visited[child] = false;
      }
    }
  }

}
