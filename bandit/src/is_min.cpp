#include "IsMin.h"
#define CLASS IsMin

#include "Calculator.h" 

Pattern CLASS::convert(Pattern& pattern) {
	// std::cout << "debug isMin: " << pattern << std::endl; // debug
	if (pattern.size() == 1) return pattern;
	Graph g = toGraph(pattern);
	map<Triplet, Tracers> heap;
	EdgeTracer cursor;
	for (ID vid = 0; vid < (ID) g.size(); vid++) {
		for (const auto& e : g[vid]) {
			if (e.labels.x <= e.labels.z) {
				cursor.set(vid, e.to, e.id, nullptr);
				heap[e.labels].push_back(cursor);
			}
		}
	}
	Pattern comp;
	comp.resize(1);
	auto itr = heap.begin();
	comp[0].labels = itr->first;
	comp[0].time.set(0,1);
	//pattern_ptr = &pattern;
	return minChecker(comp, g, itr->second);
}

Pattern CLASS::minChecker(Pattern& comp, Graph& g, Tracers& tracers) {
	// std::cout << "debug minChecker" << std::endl; // debug

	// build right most path
	vector<size_t> rm_path_index;
	scan_rm(comp, rm_path_index);

	//const Pattern& pattern = *pattern_ptr;
	int minlabel = comp[0].labels.x;
	int maxtoc = comp[rm_path_index[0]].time.b;

	vector<VertexPair> vpairs(comp.size());
	map<Pair, Tracers> b_heap;
	map<Pair, Tracers> f_heap;
	EdgeTracer* tracer;
	EdgeTracer cursor;
	Pair pkey;

	for (auto itr = tracers.begin(); itr != tracers.end(); itr++) {
		// an instance (a sequence of vertex pairs) as vector "vpair"
		tracer = &(*itr);
		//vector<bool> discovered(g.size());
		//vector<bool> tested(g.num_of_edges);
		vector<bool> discovered(g.size(), false); // as bool vector
		vector<bool> tested(g.num_of_edges, false); // as bool vector

		for (int i = vpairs.size() - 1; i >= 0; i--, tracer = tracer->predec) {
			vpairs[i] = tracer->vpair;
			tested[vpairs[i].id] = true;
			discovered[vpairs[i].a] = discovered[vpairs[i].b] = true;
		}

		Pair& rm_vpair = vpairs[rm_path_index[0]];

		// grow from the right most vertex
		for (size_t i = 0; i < g[rm_vpair.b].size(); i++) {
			Edge& added_edge = g[rm_vpair.b][i];
			// backward from the right most vertex
			for (size_t j = 1; j < rm_path_index.size(); j++) {
				size_t idx = rm_path_index[j];
				if (tested[added_edge.id]) continue;
				if (vpairs[idx].a != added_edge.to) continue;
				if (comp[idx].labels <= added_edge.labels.reverse()) {
					pkey.set(comp[idx].time.a, added_edge.labels.y);
					cursor.set(rm_vpair.b, added_edge.to, added_edge.id, &(*itr));
					b_heap[pkey].push_back(cursor);
				}
			}
			// forward from the right most vertex
			if (minlabel > added_edge.labels.z or discovered[added_edge.to]) continue;
			pkey.set(added_edge.labels.y, added_edge.labels.z);
			cursor.set(rm_vpair.b, added_edge.to, added_edge.id, &(*itr));
			f_heap[pkey].push_back(cursor);
		}
	}

	DFSCode dcode;
	if (b_heap.empty() == false) {
		auto b_itr = b_heap.begin();
		dcode.labels = Triplet(-1, b_itr->first.b, -1);
		dcode.time.set(maxtoc, b_itr->first.a);
		//if (dcode != pattern[comp.size()]) return false;
		comp.push_back(dcode);
		return minChecker(comp, g, b_itr->second);
	}
	if (f_heap.empty() == false) {
		auto f_itr = f_heap.begin();
		dcode.labels = Triplet(-1, f_itr->first.a, f_itr->first.b);
		dcode.time.set(maxtoc, maxtoc + 1);
		//if (dcode!= pattern[comp.size()]) return false;
		comp.push_back(dcode);
		return minChecker(comp, g, f_itr->second);
	}

	map<Pair, Tracers> ff_heap;
	int from = -1;
	//forward from the other nodes on the right most path
	for (size_t j = 0; j < rm_path_index.size(); j++) {
		int i = rm_path_index[j];
		for (auto itr = tracers.begin(); itr != tracers.end(); itr++) {
			tracer = &(*itr);
			vector<bool> discovered(g.size(), false); // as bool vector
			vector<bool> tested(g.num_of_edges, false); // as bool vector
			for (int k = vpairs.size() - 1; k >= 0; k--, tracer = tracer->predec) {
				vpairs[k] = tracer->vpair;
				tested[vpairs[k].id] = true;
				discovered[vpairs[k].a] = discovered[vpairs[k].b] = true;
			}
			Pair& from_vpair = vpairs[i];
			for (size_t k = 0; k < g[from_vpair.a].size(); k++) {
				Edge& added_edge = g[from_vpair.a][k];
				if (minlabel > added_edge.labels.z or discovered[added_edge.to]) continue;
				if (comp[i].labels <= added_edge.labels) {
					from = comp[i].time.a;
					pkey.set(added_edge.labels.y, added_edge.labels.z);
					cursor.set(from_vpair.a, added_edge.to, added_edge.id, &(*itr));
					ff_heap[pkey].push_back(cursor);
				}
			}
		}
		if (from != -1) break;
	}
	if (ff_heap.empty() == false) {
		auto ff_itr = ff_heap.begin();
		dcode.labels = Triplet(-1, ff_itr->first.a, ff_itr->first.b);
		dcode.time.set(from, maxtoc + 1);
		//if (dcode != pattern[comp.size()]) return false;
		comp.push_back(dcode);
		return minChecker(comp, g, ff_itr->second);
	}
	//return true;
	return comp;
}

