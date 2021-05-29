#include "common.h"


using namespace std;
string ans_filename, data_filename, cs_filename, query_filename;

bool find_file(const std::string file_name_base) // search for this name
{
    std::string file_name = "../build/";
    bool exist = true;
    file_name.append(file_name_base);
    file_name.append("_output.txt");
    ifstream file(file_name);
    if(!file.is_open()) exist = false;
    file.close();
    return exist;
}

void checkfile(ifstream &ansf, ifstream &dataf, ifstream &csf, ifstream &queryf){
  
  if (!ansf.is_open()) {
    std::cout << "Answer File " << ans_filename << " not found!\n";
    exit(EXIT_FAILURE);
  }
  
  if (!dataf.is_open()) {
    std::cout << "Answer File " << ans_filename << " not found!\n";
    exit(EXIT_FAILURE);
  }
  
  if (!csf.is_open()) {
    std::cout << "Answer File " << cs_filename << " not found!\n";
    exit(EXIT_FAILURE);
  }
  
  if (!queryf.is_open()) {
    std::cout << "Answer File " << query_filename << " not found!\n";
    exit(EXIT_FAILURE);
  }
}

void getcs(ifstream &csf, vector<vector<Vertex>> &cs, size_t &cs_total_vertex, vector<size_t> &cs_versize){
  char type;
  csf >> type >> cs_total_vertex;
  cs.resize(cs_total_vertex);

  for(size_t i=0; i<cs_total_vertex; i++){
    Vertex v;
    cs_versize.resize(cs_total_vertex);
    csf >> type >> v >> cs_versize[i];
    cs[i].resize(cs_versize[i]);
    for(size_t j=0; j < cs_versize[i]; j++)
      csf >> cs[v][j];        
  }

  csf.close();
}

void getans(ifstream &ansf, vector<vector<Vertex>> &ans, size_t &ans_total_vertex, int &ans_num){
  char type;
  ansf >> type >> ans_total_vertex;
  if(type != 't') exit(1);

  while(ansf >> type){
      ans.resize(ans_num+1);
      ans[ans_num].resize(ans_total_vertex);
      for(size_t i=0; i<ans_total_vertex; i++){
        ansf >> ans[ans_num][i];
      }
      ans_num++;
    }
    ansf.close();
}

void getquery(ifstream &queryf, size_t &num_vertices, vector<Label> &label_, set<Label> &label_set, vector<vector<Vertex>> &adj_list){
  char type;
  int graph_id_;
  queryf >> type >> graph_id_ >> num_vertices;
  adj_list.resize(num_vertices);
  while (queryf >> type) {
    if (type == 'v') {
      Vertex id;
      Label l;
      queryf >> id >> l;
      label_[id] = l;
      label_set.insert(l);
    } else if (type == 'e') {
      Vertex v1, v2;
      Label l;
      queryf >> v1 >> v2 >> l;
      adj_list[v1].push_back(v2);
      adj_list[v2].push_back(v1);
    }
  }
  queryf.close();
}

void getdata(ifstream &dataf, size_t &num_vertices, vector<Label> &label_, set<Label> &label_set, vector<vector<Vertex>> &adj_list){
  char type;
  int graph_id_;
  dataf >> type >> graph_id_ >> num_vertices;
  adj_list.resize(num_vertices+1);
  label_.resize(num_vertices);
  while (dataf >> type) {
    if (type == 'v') {
      Vertex id;
      Label l;
      dataf >> id >> l;
      label_[id] = l;
      label_set.insert(l);
    } else if (type == 'e') {
      Vertex v1, v2;
      Label l;
      dataf >> v1 >> v2 >> l;

      adj_list[v1].push_back(v2);
      adj_list[v2].push_back(v1);
    }
  }
  dataf.close();
}


int main(int argc, char* argv[]){

  std::string hprd_file_list[8] = {"lcc_hprd_n1", "lcc_hprd_n3", "lcc_hprd_n5", "lcc_hprd_n8", "lcc_hprd_s1", "lcc_hprd_s3", "lcc_hprd_s5", "lcc_hprd_s8"};
  std::string human_file_list[8] = {"lcc_human_n1", "lcc_human_n3", "lcc_human_n5", "lcc_human_n8", "lcc_human_s1", "lcc_human_s3", "lcc_human_s5", "lcc_human_s8"};
  std::string yeast_file_list[8] = {"lcc_yeast_n1", "lcc_yeast_n3", "lcc_yeast_n5", "lcc_yeast_n8", "lcc_yeast_s1", "lcc_yeast_s3", "lcc_yeast_s5", "lcc_yeast_s8"};
  std::vector<std::string> existing_file_list;

  for (int i=0; i<8; i++){
    if (find_file(hprd_file_list[i])) {
      existing_file_list.push_back(hprd_file_list[i]);
    }
  }
  for (int i=0; i<8; i++){
    if (find_file(human_file_list[i])) {
      existing_file_list.push_back(human_file_list[i]);
    }
  }
  for (int i=0; i<8; i++){
    if (find_file(yeast_file_list[i])) {
      existing_file_list.push_back(yeast_file_list[i]);
    }
  }

  int size = existing_file_list.size();
  for(int i=0; i<size; i++){
      std::cout << "checking for " << existing_file_list[i] << "\n" << "-----------------------------------------" << "\n";
  
  std::string file = existing_file_list[i];
  std::string ans_filename="../build/", data_filename="../build/", cs_filename="../build/", query_filename="../build/";
  
  ans_filename.append(file);
  cs_filename.append(file);
  query_filename.append(file);

  file.pop_back();
  file.pop_back();
  file.pop_back();
  data_filename.append(file);

  ans_filename.append("_output.txt");
  cs_filename.append(".cs");
  query_filename.append(".igraph");
  data_filename.append(".igraph");


  ifstream ansf(ans_filename), dataf(data_filename), csf(cs_filename), queryf(query_filename);
  checkfile(ansf, dataf, csf, queryf);

  cout << "file checked" << endl;

  char type;
  size_t ans_total_vertex, cs_total_vertex;
  vector<vector<Vertex>> ans;
  int ans_num=0;
  getans(ansf, ans, ans_total_vertex, ans_num);

  cout << "ans checked" << endl;

  vector<vector<Vertex>> cs;
  vector<size_t> cs_versize;
  getcs(csf, cs, cs_total_vertex, cs_versize);
  
  cout << "cs checked" << endl;
  
  size_t query_num_vertices;
  vector<Label> query_label_(ans_total_vertex);
  set<Label> query_label_set;
  vector<vector<Vertex>> query_adj_list(ans_total_vertex);
  getquery(queryf, query_num_vertices, query_label_, query_label_set, query_adj_list);
  
  cout << "query checked" << endl;
  
  size_t data_num_vertices;
  vector<Label> data_label_;
  set<Label> data_label_set;
  vector<vector<Vertex>> data_adj_list;
  getdata(dataf, data_num_vertices, data_label_, data_label_set, data_adj_list);  

  
  cout << "data checked" << endl;
  
  if(ans_total_vertex != cs_total_vertex || ans_total_vertex != query_num_vertices) {
    cout << "vertices number not the same ";
    exit(1);
  }

  
  cout << "vertex number checked" << endl;

  for(int i=0; i<ans_num; i++){
    //check one-to-one
    for(int j=0; j<ans_total_vertex; j++){
      for(int k=0; k<ans_total_vertex; k++){
        if(k!=j && ans[i][j] == ans[i][k]) {
          cout << "answer is not one-to-one";
          exit(2);
        }
      }
    }
    //check edges connected
    for(Vertex u=0; u<ans_total_vertex; u++){
      //query-edge(u, v) && cs-edge(u, v)
      //for each vertex, check the edge
      for(int k=0; k<query_adj_list[u].size(); k++){  //u-v
        Vertex v = query_adj_list[u][k];
        Vertex cs_u = ans[i][u];
        Vertex cs_v = ans[i][v];
        bool flag = false;
        for(int j=0; j<data_adj_list[cs_u].size(); j++){
          if(data_adj_list[cs_u][j] == cs_v) flag = true;
        }
        if(!flag) {
          cout << "there is no edge in the graph one-to-one " << endl;
          cout << "query:" << u <<" " << v << " " << "data" << cs_u << " " << cs_v << endl;
          exit(3);
        }
      }
    }

    //check label? -> cs
    for(int j=0; j<ans_total_vertex; j++){
      bool flag = false;
      for(int k=0; k<cs_versize[j]; k++){
        if(cs[j][k] == ans[i][j]) flag = true;
      }

      if(!flag) {
        cout << "answer " << i << "of vertex" << j << " is not in the cs ";
        exit(4);
      }
    }
  }
  cout << "well done!"<< endl;
  std::cout << "\n";
  }

  return 0;
}