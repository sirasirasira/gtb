#include "Gspan.h"
#define CLASS Gspan

#include "Database.h"
extern Database db;

/*
void CLASS::run() {
	auto& gdata = db.gdata;
	map<Triplet, GraphToTracers> heap;
	for (ID gid = 0; gid < (ID) gdata.size(); gid++) {
		EdgeTracer cursor;
		Graph& g = gdata[gid];
		for (ID vid = 0; vid < (ID) g.size(); vid++) {
			for (auto e : g[vid]) {
				if (e.labels.x <= e.labels.z) {
					cursor.set(vid, e.to, e.id, nullptr);
					heap[e.labels][gid].push_back(cursor);
				}
			}
		}
	}
	pattern.resize(1);
	for (auto itr = heap.begin(); itr != heap.end(); itr++) {
		pattern[0].labels = itr->first;
		pattern[0].time.set(0, 1);
		edgeGrow(itr->second);
	}
}

void CLASS::run(Pattern _pattern) {
	// std::cout << "debug frontier | " << _pattern  << std::endl; // debug
	pattern = _pattern;
	edgeGrow(cache[_pattern].g2tracers, true);
}
*/

void CLASS::makeRoot() {
	DFSCode dcode;
	dcode.labels = Triplet(0, 0, 0);
	dcode.time.set(0, 0);
	Pattern _root{dcode};
	root = _root;
	vector<Pattern> childs; // empty vec
	cache.insert({root, CacheRecord(childs)});

	auto& gdata = db.gdata;
	map<Triplet, GraphToTracers> heap;
	for (ID gid = 0; gid < (ID) gdata.size(); gid++) {
		EdgeTracer cursor;
		Graph& g = gdata[gid];
		for (ID vid = 0; vid < (ID) g.size(); vid++) {
			for (auto e : g[vid]) {
				if (e.labels.x <= e.labels.z) {
					cursor.set(vid, e.to, e.id, nullptr);
					heap[e.labels][gid].push_back(cursor);
				}
			}
		}
	}
	pattern.resize(1);

	for (auto itr = heap.begin(); itr != heap.end(); itr++) {
		if (support(itr->second) < minsup) {
			continue;
		}
		pattern[0].labels = itr->first;
		pattern[0].time.set(0, 1);
		cache[root].childs.push_back(pattern);
		cache.insert({pattern, CacheRecord(itr->second, childs)});
	}
}

// not minDFS
tuple<Pattern, EdgeTracer, ID> CLASS::oneEdgeSimulation(tuple<Pattern, EdgeTracer, ID>& pat_et_gid){
	// if there is no representation of expansion, pattern is not change and edgetracer.predec is nullptr
	// std::cout << "debug oneEdgeSimulation" << std::endl; // debug
	Pattern pattern = std::get<0>(pat_et_gid);
	EdgeTracer* tracer = &(std::get<1>(pat_et_gid));
	ID gid = std::get<2>(pat_et_gid);

	vector<VertexPair> vpairs(pattern.size());
	EdgeTracer cursor;

	Graph& g = db.gdata[gid];
	size_t maxtoc = 0;
	vector<bool> tested(g.num_of_edges); // edge
	vector<bool> discovered(g.size()); // node
	vector<int> vid2time(g.size(), -1);

	for (int i = vpairs.size()-1; i >= 0; --i, tracer = tracer->predec) {
		vpairs[i] = tracer->vpair;
		ID vidbase = vpairs[i].id - (vpairs[i].id % 2); // hit to_edge and from_edge
		tested[vidbase + 0] = true;
		tested[vidbase + 1] = true;
		discovered[vpairs[i].a] = discovered[vpairs[i].b] = true;

		vid2time[vpairs[i].b] = pattern[i].time.b;
		if (i == 0) {
			vid2time[vpairs[i].a] = pattern[i].time.a;
		}
		if (maxtoc < pattern[i].time.b) {
			maxtoc = pattern[i].time.b;
		}
	}

	bool valid_flg = false;
	DFSCode dcode;
	for (auto& i : Dice::shuffle(g.size())) {
		if (!discovered[i]) {
			continue;
		}
		for (auto& j : Dice::shuffle(g[i].size())) {
			Edge& added_edge = g[i][j];
			if (discovered[added_edge.to]) {
				// backward
				if (!tested[added_edge.id] and vid2time[i] > vid2time[added_edge.to]) {
					dcode.labels = Triplet(-1, added_edge.labels.y, -1);
					dcode.time.set(vid2time[i], vid2time[added_edge.to]);
					pattern.push_back(dcode);
					cursor.set(i,added_edge.to,added_edge.id,&(*tracer));
					valid_flg = true;
					break;
				}
			} else {
				// forward
				dcode.labels = Triplet(-1, added_edge.labels.y, added_edge.labels.z);
				dcode.time.set(vid2time[i], maxtoc);
				pattern.push_back(dcode);
				cursor.set(i,added_edge.to,added_edge.id,&(*tracer));
				valid_flg = true;
				break;
			}
		}
		if (valid_flg) {
			break;
		}
	}
	if (!valid_flg) {
		tracer->predec = nullptr;
	}
	return make_tuple(pattern, cursor, gid);
}

size_t CLASS::support(GraphToTracers& g2tracers) {
	size_t support = 0;
	for (auto x : g2tracers) {
		auto& id = x.first;
		if (id > (ID) db.gdata.num_train - 1) {
			break;
		}
		support++;
	}
	return support;
}

int CLASS::scanGspan(GraphToTracers& g2tracers, PairSorter& b_heap, map<int, PairSorter, std::greater<int>>& f_heap) const {
	// build right most path
	vector<size_t> rm_path_index;
	scan_rm(pattern, rm_path_index);

	int maxtoc = pattern[rm_path_index[0]].time.b;

	vector<VertexPair> vpairs(pattern.size());
	int minlabel = pattern[0].labels.x;
	EdgeTracer* tracer;

	Pair pkey;
	EdgeTracer cursor;

	for (auto x = g2tracers.begin(); x != g2tracers.end(); x++) {
		ID gid = x->first;
		Graph& g = db.gdata[gid];
		for (auto itr = x->second.begin(); itr != x->second.end(); itr++ ) {
			// an instance (a sequence of vertex pairs) as vector "vpair"
			tracer = &(*itr);
			vector<char> discovered(g.size(), false); // as bool vector
			vector<char> tested(g.num_of_edges, false); // as bool vector

			for (int i = vpairs.size() - 1; i >= 0; i--, tracer = tracer->predec) {
				vpairs[i] = tracer->vpair;
				tested[vpairs[i].id] = true;
				discovered[vpairs[i].a] = discovered[vpairs[i].b] = true;
			}

			Pair& rm_vpair = vpairs[rm_path_index[0]];

			for (size_t i = 0; i < g[rm_vpair.b].size(); i++) {
				Edge& added_edge = g[rm_vpair.b][i];
				// backward from the right most vertex
				for (size_t j = 1; j < rm_path_index.size(); j++) {
					int idx = rm_path_index[j];
					if (tested[added_edge.id]) continue;
					if (vpairs[idx].a != added_edge.to) continue;
					if (pattern[idx].labels <= added_edge.labels.reverse()) {
						pkey.set(pattern[idx].time.a, added_edge.labels.y);
						cursor.set(rm_vpair.b, added_edge.to, added_edge.id, &(*itr));
						b_heap[pkey][gid].push_back(cursor);
					}
				}
				// forward from the right most vertex
				if (minlabel > added_edge.labels.z or discovered[added_edge.to]) continue;
				pkey.set(added_edge.labels.y, added_edge.labels.z);
				cursor.set(rm_vpair.b, added_edge.to, added_edge.id, &(*itr));
				f_heap[maxtoc][pkey][gid].push_back(cursor);
			}
			// forward from the other nodes on the right most path
			for (size_t j = 0; j < rm_path_index.size(); j++) {
				size_t i = rm_path_index[j];
				Pair& from_vpair = vpairs[i];
				for (size_t k = 0; k < g[from_vpair.a].size(); k++) {
					Edge& added_edge = g[from_vpair.a][k];
					if (minlabel > added_edge.labels.z or discovered[added_edge.to]) continue;

					if (pattern[i].labels <= added_edge.labels) {
						pkey.set(added_edge.labels.y, added_edge.labels.z);
						cursor.set(from_vpair.a, added_edge.to, added_edge.id, &(*itr));
						f_heap[pattern[i].time.a][pkey][gid].push_back(cursor);
					}
				}
			}
		}
	}
	return maxtoc;
}


/*
void CLASS::edgeGrow(GraphToTracers& g2tracers, bool in_cache_flg) {
	// std::cout << "debug edgeGrow" << std::endl; // debug
	if (pattern.size() > maxpat) {
		// std::cout << "debug cut size" << std::endl; // debug
		return;
	}
	if (support(g2tracers) < minsup) {
		// std::cout << "debug cut minsup" << std::endl; // debug
		return;
	}
	if (in_cache_flg) {
		goto G_skip;
	}
	if (is_min.run(pattern) == false) {
		// std::cout << "debug cut is_min" << std::endl; // debug
		return;
	}

	//report(g2tracers);
	cache.insert({pattern, CacheRecord(g2tracers)});
G_skip:

	vector<ID> posi = getPosiIds(g2tracers);
	spliter->update(pattern, posi); 
	if (spliter->isBounded(posi)) {
		// std::cout << "debug cut bounded" << std::endl; // debug
		return;
	}
	// std::cout << "debug new ptn" << std::endl; // debug


	PairSorter b_heap;
	map<int, PairSorter, std::greater<int>> f_heap;
	int maxtoc = scanGspan(cache[pattern].g2tracers, b_heap, f_heap);
	// std::cout << "debug maxtoc " << maxtoc << std::endl; // debug

	// projecting
	DFSCode dcode;
	for (auto itr = b_heap.begin(); itr != b_heap.end(); itr++) {
		 // std::cout << "debug edgeGrow b_heap" << std::endl; // debug
		dcode.labels = Triplet(-1, itr->first.b, -1);
		dcode.time.set(maxtoc, itr->first.a);
		pattern.push_back(dcode);
		edgeGrow(itr->second);
		pattern.pop_back();
	}

	for (auto itr = f_heap.begin(); itr != f_heap.end(); itr++) {
		 // std::cout << "debug edgeGrow f_heap" << std::endl; // debug
		for (auto itr2 = itr->second.begin(); itr2 != itr->second.end(); itr2++) {
			dcode.labels = Triplet(-1, itr2->first.a, itr2->first.b);
			dcode.time.set(itr->first, maxtoc + 1);
			pattern.push_back(dcode);
			edgeGrow(itr2->second);
			pattern.pop_back();
		}
	}
}

void CLASS::report(GraphToTracers& tracers) {
	// std::cout << "debug report" << std::endl; // debug
	extern Setting setting;
	if (setting.out_instances) {
		cout << tracers.size() << "||" << pattern << "||";
		for (auto x : tracers) {
			cout << x.first << " ";
		}
		cout << endl;
	} else {
		// cout << tracers.size() << " " << ":" << pattern << endl;
		cout << " " << ":" << pattern << endl;
	}
}
*/
