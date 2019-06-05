#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <deque>
#include "finder.h"
#include "getopt.h"

using std::deque;
using std::vector;
using std::string;

using std::istringstream;

using std::unordered_map;

#define USAGE " features graphs"

void gtraverse(int node_id, vector<DFSCode>& query, unsigned pnt, int cur,
		GraphToTracers& g2tracers, vector<Graph>& gdata,
		const patricia& P, map<int,set<int> >& saver, map<int,set<int> >& bsaver,
		int treelevel, map<int,int>& dfslabel){

	//if (g2tracers.size()==0) return;

	//std::cout << "gtraverse @ node " << node_id << ": " << g2tracers.size() << " query.size=" << query.size() << std::endl;

	if (query.size() > pnt){
		DFSCode& dfscur = query[pnt];
		int s = dfscur.time.a;
		int t = dfscur.time.b;
		EdgeTracer* tracer;
		EdgeTracer cursor;
		GraphToTracers newtracer;

		for(GraphToTracers::iterator x = g2tracers.begin(); x != g2tracers.end(); ++x){
			int gid = x->first;
			Graph& g = gdata[gid];
			for(vector<EdgeTracer>::iterator it = x->second.begin(); it != x->second.end(); ++it){

				tracer = &(*it);

				vector<unsigned> discovered(g.size(),0);
				int growpoint = -999;
				int backpoint = -999;
				for(int j = pnt-1; j>=0; --j, tracer = tracer->predec){
					discovered[tracer->vpair.a] = 1;
					discovered[tracer->vpair.b] = 1;
					if(query[j].time.a == s) growpoint = tracer->vpair.a;
					if(query[j].time.b == s) growpoint = tracer->vpair.b;
					if(s > t){
						if(query[j].time.a == t) backpoint = tracer->vpair.a;
						if(query[j].time.b == t) backpoint = tracer->vpair.b;
					}
				}

				for(unsigned int i=0; i<g[growpoint].size(); ++i){
					const Edge& e = g[growpoint][i];
					if (e.labels.y == dfscur.labels.y && e.labels.z == dfscur.labels.z){
						if (s < t){
							if (discovered[e.to] != 1){
								cursor.vpair.a = growpoint;
								cursor.vpair.b = e.to;
								cursor.vpair.id = e.id;
								cursor.predec  = &(*it);
								newtracer[gid].push_back(cursor);
							}
						}else{
							if (backpoint == e.to){
								cursor.vpair.a = growpoint;
								cursor.vpair.b = e.to;
								cursor.vpair.id = e.id;
								cursor.predec  = &(*it);
								newtracer[gid].push_back(cursor);
							}
						}
					}
				}
			}
		}
		gtraverse(node_id,query,pnt+1,cur,newtracer,gdata,P,saver,bsaver,treelevel,dfslabel);
	}else{

		if(P.outQ[node_id]>0){
			//std::cout << "==>" << g2tracers.size() << " " << query << std::endl;

			for(GraphToTracers::iterator x = g2tracers.begin(); x != g2tracers.end(); ++x){
				//std::cout << "PUT saver[" << P.outQ[node_id] << "] = " << x->first << std::endl;
				//saver[x->first].insert(P.outQ[node_id]);
				saver[P.outQ[node_id]].insert(x->first);
				bsaver[x->first].insert(P.outQ[node_id]);
			}
		}

		const list<int>& child = P.children[node_id];
		for(list<int>::const_iterator it=child.begin(); it != child.end(); ++it){

			const vector<string>& vec = P.node[*it];
			int addcount = 0;
			int c2 = cur;
			vector<int> addkeys;
			for(vector<string>::const_iterator it2=vec.begin(); it2 != vec.end(); ++it2){
				c2 = buildDFSCode(*it2, query, c2, dfslabel, addkeys);
				addcount += 1;
			}
			gtraverse(*it,query,pnt,c2,g2tracers,gdata,P,saver,bsaver,treelevel+1,dfslabel);
			for(vector<int>::iterator lit=addkeys.begin(); lit!=addkeys.end(); ++lit){
				//std::cout << "erase" << *lit << std::endl;
				dfslabel.erase(*lit);
			}
			for(int j=0; j<addcount; ++j){
				query.pop_back();
			}
		}
	}
}


int main(int argc, char **argv) {

	std::string tran_feat_fname = "tran_by_feat.out";
	std::string feat_tran_fname = "feat_by_tran.out";

	static struct option long_options[] =
	{
		{"tran_by_feat", required_argument, 0, 0}, // index = 0
		{"feat_by_tran", required_argument, 0, 0}  // index = 1
	};

	int option_index = 0;

	int opt;
	string outdir=".";
	while ((opt = getopt_long (argc, argv, "d:o:", 
					long_options, &option_index)) != -1) {
		switch (opt) {
			case 0:
				if (long_options[option_index].flag != 0)
					break;
				if(option_index==0){
					tran_feat_fname = optarg;
				}else if(option_index==1){
					feat_tran_fname = optarg;
				}
				break;
			case 'o':
				outdir = string(optarg);
				break;
			default:
				std::cerr << "Usage: "<< argv[0] << USAGE << std::endl;
				return -1;
		}
	}

	// *********** i/o check ***********

	if(argc-optind != 2){
		std::cerr << "Usage: "<< argv[0] << USAGE << std::endl;
		return -1;
	}

	// *********** mode check ***********

	Finder f;

	string filefeat   (argv[optind]);
	string filegraph  (argv[optind+1]);

	// *********** file check ***********

	std::ifstream feat_file(filefeat.c_str());
	if(feat_file.fail()){
		std::cerr << "File not found: " << filefeat << std::endl;
		return -1;
	}
	std::ifstream graph_file(filegraph.c_str());
	if(graph_file.fail()){
		std::cerr << "File not found: " << filegraph << std::endl;
		return -1;
	}

	std::cout << "Finder version 0.4"         << std::endl;
	std::cout << "      by Ichigaku Takigawa" << std::endl;
	std::cout << ">settings are:"             << std::endl;
	std::cout << "  Feature => " << filefeat   << std::endl;
	std::cout << "  Graph   => " << filegraph  << std::endl;

	// *********** main ***********

	f.read_features(feat_file);
	f.read_graphs(graph_file);

	std::cout << "# graphs: " << f.gdata.size() << std::endl;

	map<int,set<int> > gsave,gbacksave;

	vector<DFSCode> gquery;
	map<int,int> dfslabel;

	list<int>& childg = f.gpat_tree.children[patricia::root];
	for(list<int>::iterator it=childg.begin(); it != childg.end(); ++it){
		//std::cout << "tree: # " << *it << std::endl;
		vector<string>& vec = f.gpat_tree.node[*it];
		int cur = 0;
		int addcount = 0;
		vector<int> addkeys;
		for(vector<string>::iterator it2=vec.begin(); it2 != vec.end(); ++it2){
			cur = buildDFSCode(*it2, gquery, cur, dfslabel, addkeys);
			addcount += 1;
		}
		GraphToTracers& newtracer = f.gheap[gquery[0].labels];

		gtraverse(*it,gquery,1,cur,newtracer,f.gdata,f.gpat_tree,gsave,gbacksave,1,dfslabel);
		for(vector<int>::iterator lit=addkeys.begin(); lit!=addkeys.end(); ++lit){
			//std::cout << "erase" << *lit << std::endl;
			dfslabel.erase(*lit);
		}
		for(int j=0; j<addcount; ++j){
			gquery.pop_back();
		}
	}

	std::cout << "gsave size: " << gsave.size() << std::endl;

	std::ofstream ofs1(feat_tran_fname.c_str());

	for( map<int,set<int> >::iterator itr = gsave.begin(); itr != gsave.end(); ++itr){
		//std::cout << itr->second.size() << " " << f.fg_name[itr->first] << std::endl;
		//std::cout << itr->second.size() << std::endl;
		ofs1 << itr->first << ":";
		for( set<int>::iterator itr2 =itr->second.begin(); itr2 != itr->second.end(); ++itr2){
			ofs1 << *itr2+1 << " ";
		}
		ofs1 << std::endl;
	}
	ofs1.close();

	std::ofstream ofs2(tran_feat_fname.c_str());

	for( map<int,set<int> >::iterator itr = gbacksave.begin(); itr != gbacksave.end(); ++itr){
		ofs2 << itr->first+1 << " " << f.gdata[itr->first].value << ":";
		for( set<int>::iterator itr2 =itr->second.begin(); itr2 != itr->second.end(); ++itr2){
			ofs2 << *itr2 << " ";
		}
		ofs2 << std::endl;
	}
	ofs2.close();

	return 0;
}


//  

