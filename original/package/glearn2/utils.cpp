// gspan.cpp
#include "glearn.h"
#include "lossfunc.h"
#include <algorithm>
#include <sstream>
#include <fstream>
#include <set>

using std::map;
using std::vector;
using std::set;

bool Glearn::node_quality_check(){
  bool flag = true;
  if ((int) loaded_dfscode.size() > maxpat){
    flag = false;
  }else  if (!is_min_dfs()){
    flag = false;
  }
  return flag;
}

void Glearn::dfs_to_str(){
  std::ostringstream os;
  os << loaded_dfscode;
  loaded_dfsstr = os.str();
}

int Glearn::explore_child_nodes(Graph_to_VpairPList& g2tracers,
				map<Pair,Graph_to_VpairPList>& b_heap,
				map<int,map<Pair,Graph_to_VpairPList>,greater<int> >& f_heap){
  // build right most path
  vector<int> rm_path_index;

  get_rightmost_path(loaded_dfscode,rm_path_index);
	
  int maxtoc = loaded_dfscode[rm_path_index[0]].time.b;
	
  vector<VertexPair> vpairs(loaded_dfscode.size());
  int minlabel = loaded_dfscode[0].labels.x;
  VpairP *tracer;
    
  Pair pkey;
    
  for(Graph_to_VpairPList::iterator it = g2tracers.begin(); it != g2tracers.end(); ++it){
    int gid = it->first;
    Graph& g = gdata[gid];
    VpairPList& vpl = it->second;
    for(VpairPList::iterator it2 = vpl.begin(); it2 != vpl.end(); ++it2){
      
      // an instance (a sequence of vertex pairs) as vector "vpair"
      tracer = &(*it2);
            
      vector<bool> discovered(g.size());
      vector<bool> tested(g.num_of_edges);

      int i = vpairs.size()-1;
      for(VpairP* p = tracer; i>=0; --i, p = p->predec){
	vpairs[i] = p->vpair;
	tested[vpairs[i].id] = true;
	discovered[vpairs[i].a] = discovered[vpairs[i].b] = true;
      }
			
      Pair& rm_vpair = vpairs[rm_path_index[0]];
            
      for(unsigned int i=0; i<g[rm_vpair.b].size(); ++i){
	Edge& added_edge = g[rm_vpair.b][i];
	//backward from the right most vertex
	for(unsigned int j=1; j<rm_path_index.size(); ++j){
	  int idx = rm_path_index[j];
	  if(tested[added_edge.id]) continue;
	  if(vpairs[idx].a != added_edge.to) continue;
	  if(loaded_dfscode[idx].labels <= added_edge.labels.reverse()){
	    pkey.set(loaded_dfscode[idx].time.a,added_edge.labels.y);
	    b_heap[pkey][gid].push(rm_vpair.b,added_edge.to,added_edge.id,tracer);
	  }
	}
	// forward from the right most vertex
	if (minlabel > added_edge.labels.z || discovered[added_edge.to]){ continue; }
	pkey.set(added_edge.labels.y,added_edge.labels.z);
	f_heap[maxtoc][pkey][gid].push(rm_vpair.b,added_edge.to,added_edge.id,tracer);
      }
			
      // forward from the other nodes on the right most path
      for(unsigned int j=0; j<rm_path_index.size(); ++j){
	int i = rm_path_index[j];
	Pair& from_vpair = vpairs[i];
	for(unsigned int k=0; k<g[from_vpair.a].size(); ++k){
	  Edge& added_edge = g[from_vpair.a][k];
	  if (minlabel > added_edge.labels.z || discovered[added_edge.to]){ continue; }
					
	  if(loaded_dfscode[i].labels <= added_edge.labels){
	    pkey.set(added_edge.labels.y,added_edge.labels.z);
	    f_heap[loaded_dfscode[i].time.a][pkey][gid].push(from_vpair.a,added_edge.to,added_edge.id,tracer);
	  }
	}
      }
    }
  }
  return maxtoc;
    
}

void Glearn::set_data(std::istream& is) {
  gdata = read_graphs(is);
};

void Glearn::get_rightmost_path(vector<EdgeCode>& pat, vector<int>& idx){
  int prev = -1;
    
  for(int i=pat.size()-1; i>=0; --i){
    if(pat[i].time.a < pat[i].time.b){ // forward edge
      if(prev == pat[i].time.b || prev == -1){
	idx.push_back(i);
	prev = pat[i].time.a;
      }
    }
  }
}

bool Glearn::is_min_dfs(){
  if(loaded_dfscode.size() == 1) return true;
    
  // do test
  Graph g = toGraph(loaded_dfscode);
  map<Triplet,VpairPList> heap;
  for(unsigned int v=0; v<g.size(); ++v){
    for(vector<Edge>::iterator e = g[v].begin(); e != g[v].end(); ++e){
      if (e->labels.x <= e->labels.z){
	heap[e->labels].push(v,e->to,e->id,0);
      }
    }
  }
	
  vector<EdgeCode> comp;
  comp.resize(1);
  map<Triplet,VpairPList>::iterator it = heap.begin();
    
  comp[0].labels = it->first;
  comp[0].time.set(0,1);
    
  bool res = min_checker(comp, g, it->second);

  return res;
}

bool Glearn::min_checker(vector<EdgeCode>& comp, Graph& g, VpairPList& _tracers){
  // build right most path
  vector<int> rm_path_index;
  get_rightmost_path(comp,rm_path_index);
    
  int minlabel = comp[0].labels.x;
  int maxtoc = comp[rm_path_index[0]].time.b;
  VpairPList& tracers =_tracers;
	
  vector<VertexPair> vpairs(comp.size());
  map<Pair,VpairPList>  b_heap, f_heap;
  VpairP *tracer;
  Pair pkey;
	
  for(VpairPList::iterator it = tracers.begin(); it != tracers.end(); ++it){
    // an instance (a sequence of vertex pairs) as vector "vpair"
    tracer = &(*it);
    vector<bool> discovered(g.size());
    vector<bool> tested(g.num_of_edges);

    int i = vpairs.size()-1;
    for(VpairP* p = tracer; i>=0; --i, p = p->predec){
      vpairs[i] = p->vpair;
      tested[vpairs[i].id] = true;
      discovered[vpairs[i].a] = discovered[vpairs[i].b] = true;
    }
    
    Pair& rm_vpair = vpairs[rm_path_index[0]];
		
    // grow from the right most vertex
    for(unsigned int i=0; i<g[rm_vpair.b].size(); ++i){
      Edge& added_edge = g[rm_vpair.b][i];
      //backward from the right most vertex
      for(unsigned int j=1; j<rm_path_index.size(); ++j){
	int idx = rm_path_index[j];
	if(tested[added_edge.id]) continue;
	if(vpairs[idx].a != added_edge.to) continue;
	if(comp[idx].labels <= added_edge.labels.reverse()){
	  pkey.set(comp[idx].time.a,added_edge.labels.y);
	  b_heap[pkey].push(rm_vpair.b,added_edge.to,added_edge.id,&(*it));
	}
      }
      // forward from the right most vertex
      if (minlabel > added_edge.labels.z || discovered[added_edge.to]){ continue; }
      pkey.set(added_edge.labels.y,added_edge.labels.z);
      f_heap[pkey].push(rm_vpair.b,added_edge.to,added_edge.id,&(*it));
    }
  }
	
  EdgeCode dcode;
  map<Pair,VpairPList>::iterator b_it=b_heap.begin();
  if(b_it != b_heap.end()){
    dcode.labels = Triplet(-1,b_it->first.b,-1);
    dcode.time.set(maxtoc,b_it->first.a);
    if(dcode != loaded_dfscode[comp.size()]) return false;
    comp.push_back(dcode);
    return min_checker(comp, g,b_it->second);
  }
    
  map<Pair,VpairPList>::iterator f_it=f_heap.begin();
  if(f_it != f_heap.end()){
    dcode.labels = Triplet(-1,f_it->first.a,f_it->first.b);
    dcode.time.set(maxtoc,maxtoc+1);
    if(dcode != loaded_dfscode[comp.size()]) return false;
    comp.push_back(dcode);
    return min_checker(comp, g,f_it->second);
  }
	
  map<Pair,VpairPList> ff_heap;
  int from = -1;
  // forward from the other nodes on the right most path
  for(unsigned int j = 0; j < rm_path_index.size(); ++j){
    int i = rm_path_index[j];
    for(VpairPList::iterator it =tracers.begin(); it != tracers.end(); ++it){
      tracer = &(*it);
      vector<bool> discovered(g.size());
      vector<bool> tested(g.num_of_edges);
            
      int k = vpairs.size()-1;
      for(VpairP* p = tracer; k>=0; --k, p = p->predec){
	vpairs[k] = p->vpair;
	tested[vpairs[k].id] = true;
	discovered[vpairs[k].a] = discovered[vpairs[k].b] = true;
      }
			
      Pair& from_vpair = vpairs[i];
      for(unsigned int k=0; k<g[from_vpair.a].size(); ++k){
	Edge& added_edge = g[from_vpair.a][k];
	if (minlabel > added_edge.labels.z || discovered[added_edge.to]) continue;
	if(comp[i].labels <= added_edge.labels){
	  from = comp[i].time.a;
	  pkey.set(added_edge.labels.y,added_edge.labels.z);
	  ff_heap[pkey].push(from_vpair.a,added_edge.to,added_edge.id,&(*it));
	}
      }
    }
    if(from != -1) break;
  }
	
  map<Pair,VpairPList>::iterator ff_it=ff_heap.begin();
  if(ff_it != ff_heap.end()){
    dcode.labels = Triplet(-1,ff_it->first.a,ff_it->first.b);
    dcode.time.set(from,maxtoc+1);
    if(dcode != loaded_dfscode[comp.size()]) return false;
    comp.push_back(dcode);
    return min_checker(comp, g,ff_it->second);
  }
  return true;
}


/////////////

const vector<Graph> read_graphs(std::istream &is) {
  vector<Graph> graphs;
  Graph g;
  Triplet labels;
  Edge edge;
    
  double val;
  char c; unsigned int i, j;
  std::string line, gname;
  int eid=0;
    
  while (getline(is, line)) {
    if (line.empty()) {
      g.num_of_edges = eid;
      graphs.push_back(g);
    }
    std::stringstream stream(line);
    if (line[0] == 't') {
      g.clear();
      eid = 0;
      stream >> c >> c >> c >> val >> gname;
      g.name = gname;
      g.value = val;
    } else if (line[0] == 'v') {
      stream >> c >> i >> j;
      g.resize(i+1);
      g.label[i] = j;
    } else if (line[0] == 'e') {
      stream >> c >> i >> j >> labels.y;
      labels.x = g.label[i];
      labels.z = g.label[j];
      edge.to = j; edge.labels = labels; edge.id = eid++;
      g[i].push_back(edge);
      edge.to = i; edge.labels = labels.reverse(); edge.id = eid++;
      g[j].push_back(edge);
    }
  }
    
  return graphs;
}

Graph toGraph(vector<EdgeCode>& dfscode){
  Graph g;
  Edge edge;
  int eid = 0;
  for(vector<EdgeCode>::iterator p=dfscode.begin(); p!=dfscode.end(); ++p){
    g.resize(std::max(p->time.a,p->time.b)+1);
    if(p->labels.x != -1) g.label[p->time.a] = p->labels.x;
    if(p->labels.z != -1) g.label[p->time.b] = p->labels.z;
    edge.to = p->time.b;
    edge.labels.x = g.label[p->time.a];
    edge.labels.y = p->labels.y;
    edge.labels.z = g.label[p->time.b];
    edge.id = eid++;
    g[p->time.a].push_back(edge);
    edge.to = p->time.a; edge.labels = edge.labels.reverse(); edge.id = eid++;
    g[p->time.b].push_back(edge);
  }
  g.num_of_edges = eid;
  return g;
}

unsigned int support(Graph_to_VpairPList& g2tracers){
  //return g2tracers.size();
  int sum=0;
  for(Graph_to_VpairPList::iterator x = g2tracers.begin(); x != g2tracers.end(); ++x){
    sum += 1;
  }
  return sum;
}

