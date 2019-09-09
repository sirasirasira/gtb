#include "Gspan.h"
#define CLASS Gspan

#include "Database.h"
extern Database db;

#include "Calculator.h" 
#include "Spliter.h" 

void CLASS::run() {
	// std::cout << "debug" << std::endl; // debug
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
		if (!check_pattern(pattern, itr->second)) {
			continue;
		}
		e1patterns.push_back(pattern);
		vector<DFSCode> childs;
		cache.insert({pattern, CacheRecord(itr->second, childs)});
		vector<ID> posi = getPosiIds(itr->second);
		spliter->update(pattern, posi);
		if (!spliter->isBounded(posi)) {
			edgeGrow();
		}
	}
}

void CLASS::run(Pattern _pattern) {
	// std::cout << "debug frontier | " << _pattern  << std::endl; // debug
	pattern = _pattern;
	edgeGrow();
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

bool CLASS::check_pattern(Pattern pattern, GraphToTracers& g2tracers) {
		if (pattern.size() > maxpat) {
			return false;
		}
		if (support(g2tracers) < minsup) {
			return false;
		}
		if (is_min.run(pattern) == false) {
			return false;
		}
		return true;
}

void CLASS::edgeGrow() {
	// std::cout << "debug edgeGrow" << std::endl; // debug
	// report(g2tracers);

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
		if (!check_pattern(pattern, itr->second)) {
			pattern.pop_back();
			continue;
		}
		vector<DFSCode> childs;
		cache.insert({pattern, CacheRecord(itr->second, childs)});
		vector<ID> posi = getPosiIds(itr->second);
		spliter->update(pattern, posi);
		if (!spliter->isBounded(posi) and pattern.size() < maxpat) {
			edgeGrow();
		}
		pattern.pop_back();
		cache[pattern].childs.push_back(dcode);
	}

	for (auto itr = f_heap.begin(); itr != f_heap.end(); itr++) {
		// std::cout << "debug edgeGrow f_heap" << std::endl; // debug
		for (auto itr2 = itr->second.begin(); itr2 != itr->second.end(); itr2++) {
			dcode.labels = Triplet(-1, itr2->first.a, itr2->first.b);
			dcode.time.set(itr->first, maxtoc + 1);
			pattern.push_back(dcode);
			if (!check_pattern(pattern, itr2->second)) {
				pattern.pop_back();
				continue;
			}
			vector<DFSCode> childs;
			cache.insert({pattern, CacheRecord(itr2->second, childs)});
			vector<ID> posi = getPosiIds(itr2->second);
			spliter->update(pattern, posi);
			if (!spliter->isBounded(posi) and pattern.size() < maxpat) {
				edgeGrow();
			}
			pattern.pop_back();
			cache[pattern].childs.push_back(dcode);
		}
	}
}

int CLASS::scanGspan(GraphToTracers& g2tracers, PairSorter& b_heap, map<int, PairSorter, std::greater<int>>& f_heap) const {
	cache[pattern].scan = true;
	// build right most path
	vector<size_t> rm_path_index;
	scan_rm(pattern, rm_path_index);

	int maxtoc = pattern[rm_path_index[0]].time.b;

	vector<VertexPair> vpairs(pattern.size());
	int minlabel = pattern[0].labels.x;
	EdgeTracer* tracer;

	Pair pkey;
	EdgeTracer cursor;

	// std::cout << pattern << std::endl; // debug
	for (auto x = g2tracers.begin(); x != g2tracers.end(); x++) {
		ID gid = x->first;
		Graph& g = db.gdata[gid];
		for (auto itr = x->second.begin(); itr != x->second.end(); itr++ ) {
			// an instance (a sequence of vertex pairs) as vector "vpair"
			tracer = &(*itr);

			//vector<bool> discovered(g.size());
			//vector<bool> tested(g.num_of_edges);
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

