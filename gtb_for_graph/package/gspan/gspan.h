#ifndef GSPAN_H_
#define GSPAN_H_

#include <map>
#include <list>
#include <vector>
#include <iostream>
#include <ext/hash_map>

using std::map;
using std::list;
using std::string;
using std::vector;
using std::greater;
using __gnu_cxx::hash_map;

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
struct DFSCode { public: Triplet labels; Pair time; };

typedef vector< vector<Edge> > AdjacentList;

class Graph: public AdjacentList {
 public:
  vector<int> label;
  int num_of_edges;
  void resize(unsigned int s){ AdjacentList::resize(s); label.resize(s); }
};

class EdgeTracer {
 public:
  VertexPair  vpair;
  EdgeTracer* predec;
  explicit EdgeTracer(){};
  void set(int,int,int,EdgeTracer*);
};

inline void EdgeTracer::set(int a,int b,int id, EdgeTracer* pr){
  vpair.a = a;
  vpair.b = b;
  vpair.id = id;
  predec = pr;
}

typedef vector<EdgeTracer> Tracers;
typedef map<int,Tracers> GraphToTracers;
typedef map<Pair,GraphToTracers> PairSorter;

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

inline bool operator== (const DFSCode& left, const DFSCode& right){
  if(left.time.a != right.time.a) return false;
  if(left.time.b != right.time.b) return false;	
  if(left.labels.x != right.labels.x) return false;
  if(left.labels.y != right.labels.y) return false;
  return (left.labels.z == right.labels.z);	
}

inline bool operator!= (const DFSCode& x, const DFSCode& y){
  return !(x==y);
}

inline std::ostream& operator<< (std::ostream& os, const vector<DFSCode> pattern){
  if(pattern.empty()) return os;
  os << "(" << pattern[0].labels.x << ") " << pattern[0].labels.y << " (0f" << pattern[0].labels.z << ")";
  for(unsigned int i=1; i<pattern.size(); ++i){
    if(pattern[i].time.a < pattern[i].time.b){
      os << " " << pattern[i].labels.y << " (" << pattern[i].time.a << "f" << pattern[i].labels.z << ")";
    }else{
      os << " " << pattern[i].labels.y << " (b" << pattern[i].time.b << ")";
    }
  }
  return os;
}

const vector<Graph> readGraphs(std::istream&);

class Gspan {
 private:
  bool is_min();
  unsigned int support(GraphToTracers&);
  void scan_rm(vector<DFSCode>&, vector<int>&);
  bool min_checker(vector<DFSCode>&, Graph&, Tracers&);
  int  scan_gspan(GraphToTracers&, PairSorter&, map<int,PairSorter,greater<int> >&);
 public:
  bool out_instances;
  unsigned int minsup;
  unsigned int maxpat;
  vector<Graph>   gdata;
  vector<DFSCode> pattern;
  void set_data(std::istream& is) {
    gdata = readGraphs(is);
  };
  void edge_grow(GraphToTracers&);
  void report(GraphToTracers&);
  void run();
};

Graph toGraph(vector<DFSCode>&);

inline void init(EdgeTracer& et,int a,int b,int id){
  et.vpair.a = a;
  et.vpair.b = b;
  et.vpair.id = id;
  et.predec = 0;
}

#endif /*GSPAN_H_*/
