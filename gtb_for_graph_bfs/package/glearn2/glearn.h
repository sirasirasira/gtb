#ifndef GLEARN_H_
#define GLEARN_H_

#include <set>
#include <map>
#include <list>
#include <vector>
#include <iostream>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>

#include "tree.h"

using std::set;
using std::map;
using std::list;
using std::string;
using std::vector;
using std::greater;
using boost::unordered_map;
using boost::unordered_set;

//#include <ext/hash_map>
//using __gnu_cxx::hash_map;

struct Triplet {
public:
  int x, y, z;
  Triplet reverse();
  explicit Triplet(){}
  explicit Triplet(int _x,int _y,int _z): x(_x),y(_y),z(_z){}
};

struct Pair { public: int a, b; void set(int, int); };
inline void Pair::set(int _a, int _b){ a = _a; b = _b; }

struct VertexPair: public Pair { public: int id; };
struct Edge   { public: int to, id; Triplet labels; };
struct EdgeCode { public: Triplet labels; Pair time; };

typedef vector< vector<Edge> > AdjacentList;

// graph as "an edge list"
// vertex_id1: (to_vertex2, edge_id1, edge_label), (to_vertex2, edge_id2, edge_label2), ...
// vertex_id2: ...
//   :
// vertex_idn: ...
class Graph: public AdjacentList {
 public:
  double value;      // associated value to learn
  string name;       // name label of this graph
  vector<int> label; // vertex labels
  int num_of_edges;
  void resize(unsigned int s){ AdjacentList::resize(s); label.resize(s); }
};

// list of vertex pairs
class VpairP {
public:
    VertexPair  vpair;
    VpairP* predec;
    explicit VpairP(){};
    inline void set(int a,int b,int id, VpairP* pr){
      vpair.a = a;
      vpair.b = b;
      vpair.id = id;
      predec = pr;
    }
};

//typedef vector<VpairP> VpairPList;
class VpairPList: public list<VpairP> {
 public:
  void push(int a,int b,int id, VpairP* pr){
    resize (size()+1);
    VpairP &vp = (*this).back();
    vp.set(a,b,id,pr);
  }
};

typedef map<int,VpairPList> Graph_to_VpairPList;

//void copy_maps(Graph_to_VpairPList&, Graph_to_VpairPList&);

inline Triplet Triplet::reverse(){ return Triplet(z,y,x); }

inline bool operator< (const Triplet& left, const Triplet& right){
    if (left.x!=-1 && right.x!=-1 && left.x != right.x) return (left.x < right.x);
    if (left.y!=-1 && right.y!=-1 && left.y != right.y) return (left.y < right.y);
    return (left.z < right.z);
}

inline bool operator<= (const Triplet& left, const Triplet& right){
    return !(right < left);
}

inline bool operator== (const Triplet& left, const Triplet& right){
    return (left.x==right.x && left.y==right.y && left.z==right.z);
}

inline bool operator< (const Pair& left, const Pair& right){
    if (left.a != right.a) return (left.a < right.a);
    return (left.b < right.b);
}

inline bool operator== (const EdgeCode& left, const EdgeCode& right){
    if(left.time.a != right.time.a) return false;
    if(left.time.b != right.time.b) return false;
    if(left.labels.x != right.labels.x) return false;
    if(left.labels.y != right.labels.y) return false;
    return (left.labels.z == right.labels.z);
}

inline bool operator!= (const EdgeCode& x, const EdgeCode& y){
    return !(x==y);
}

inline std::ostream& operator<< (std::ostream& os, const vector<EdgeCode> dsfcode){
    if(dsfcode.empty()) return os;
    os << "(" << dsfcode[0].labels.x << ") " << dsfcode[0].labels.y << " (0f" << dsfcode[0].labels.z << ")";
    for(unsigned int i=1; i<dsfcode.size(); ++i){
        if(dsfcode[i].time.a < dsfcode[i].time.b){
            os << " " << dsfcode[i].labels.y << " (" << dsfcode[i].time.a << "f" << dsfcode[i].labels.z << ")";
        }else{
            os << " " << dsfcode[i].labels.y << " (b" << dsfcode[i].time.b << ")";
        }
    }
    return os;
}

inline double max(double a, double b){
  if(a > b){
    return a;
  }else{
    return b;
  }
}

inline double min(double a, double b){
  if(a < b){
    return a;
  }else{
    return b;
  }
}

const vector<Graph> read_graphs(std::istream&);
Graph toGraph(vector<EdgeCode>&);
unsigned int support(Graph_to_VpairPList&);

class Parameter {
private:
    int num_of_vars;
    unordered_map<string,int> dfscode_id;
    unordered_map<int,string> dfscode_str;
public:
    // map
    map<int,set<int> > h1;
    map<int,set<int> > h2;
    map<int,set<int> > h_tmp;
    map<int,set<int> > coef_idx;
    unordered_map<int,double>  d;
    unordered_map<int,double>  Hd;
    unordered_map<int,double>  theta;
    unordered_map<int,double>  delta;
    unordered_map<int,double>  codesize;
    unordered_map<int,set<int> >  instances;
    // set
    set<int>     nonzero_idx;
    set<int>     nonzero_d_idx;
    // vector
    vector<double>   mu;
    // scalar
    double max_d;
    //double max_Hd;
    // methods
    explicit Parameter(): num_of_vars(0) {}
    int register_node(vector<EdgeCode> &dfscode, std::string &dfsstr){
      unordered_map<string,int>::const_iterator it = dfscode_id.find(dfsstr);
      if(it != dfscode_id.end()){
	return it->second;
      }else{
	num_of_vars++;
	dfscode_id[dfsstr] = num_of_vars;
	dfscode_str[num_of_vars] = dfsstr;
	codesize[num_of_vars] = dfscode.size();
	return num_of_vars;
      }
      return 0;
    }
    string retrieve(int idx){
      unordered_map<int,string>::const_iterator it = dfscode_str.find(idx);
      if(it != dfscode_str.end()){
	return it->second;
      }else{
	std::cout << "Not found in index: " << idx << std::endl;
      }
      return "";
    }
    int get_num_vars(){
      return num_of_vars;
    }
};

class tnode {
 public: 
  int pid;
  EdgeCode ecode;
  bool has_children;
  Graph_to_VpairPList g2tracers;
};

class Glearn {
private:
  tree<tnode> search_tree;
  void get_rightmost_path(vector<EdgeCode>&, vector<int>&);
  bool is_min_dfs();
  bool min_checker(vector<EdgeCode>&, Graph&, VpairPList&);
  int  explore_child_nodes(Graph_to_VpairPList&,
			   map<Pair,Graph_to_VpairPList>&, 
			   map<int,map<Pair,Graph_to_VpairPList>,greater<int> >&);
  bool node_evaluate(Graph_to_VpairPList&);
  void postorder_op(int);
 public:
  std::string ofname;
  Parameter param;
  unordered_set< set<int> > coc;
  unordered_set<int>  irregular;
  // parameters
  double lm1;
  double lm2;
  double scaling;
  double gamma;
  double sigma;
  double c;
  double alpha;
  double gs_cnst;
  double tol;
  int maxpat;
  int maxitr;
  int numvisited;
  // bool options
  bool prune;
  bool unique;
  bool inexact;
  bool continuation;
  bool out_instances;
  // data
  vector<Graph>    gdata;
  vector<EdgeCode> loaded_dfscode;
  std::string      loaded_dfsstr;
  void dfs_to_str();
  void set_data(std::istream&);
  bool node_quality_check();
  void recursive_search_from(int);
  void run();
};


inline void init(VpairP& vplist,int a,int b,int id){
  vplist.vpair.a = a;
  vplist.vpair.b = b;
  vplist.vpair.id = id;
  vplist.predec = 0;
}


#endif /*GLEARN_H_*/

