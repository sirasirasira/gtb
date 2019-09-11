#include "Gspan.h"
#define CLASS Gspan

#include "Database.h"
extern Database db;
extern Setting setting;

void CLASS::makeRoot(const vector<ID>& targets) {
	// std::cout << "makeRoot" << std::endl; // debug
	Pattern pattern;
	DFSCode dcode;
	dcode.labels = Triplet(-1, -1, -1);
	dcode.time.set(0, 0);
	root = {dcode};
	vector<Pattern> childs; // empty vec
	cache.insert({root, CacheRecord(childs)});
	cache[root].scan = true;

	auto& gdata = db.gdata;
	map<Triplet, GraphToTracers> heap;
	for (ID gid : targets) {
		EdgeTracer cursor(1);
		Graph& g = gdata[gid];
		for (ID vid = 0; vid < (ID) g.size(); vid++) {
			for (auto e : g[vid]) {
				if (e.labels.x <= e.labels.z) {
					cursor[0] = VertexPair(vid, e.to, e.id, e.labels.y);
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
	cout << "1 edge graphs size: " << cache.size()-1 << endl;
}

// not minDFS, not support
Pattern CLASS::EdgeSimulation(const Pattern& _pattern, const EdgeTracer& _tracer, ID gid) {
	// std::cout << "EdgeSimulation: " << _pattern << std::endl; // debug

	Pattern pattern = _pattern;
	EdgeTracer cursor = _tracer;

	Graph& g = db.gdata[gid];
	size_t maxtoc = 0;
	vector<bool> tested(g.num_of_edges); // edge
	vector<bool> discovered(g.size()); // node
	vector<int> vid2time(g.size(), -1);

	for (int i = cursor.size()-1; i >= 0; --i) {
		ID eidbase = cursor[i].id - (cursor[i].id % 2); // hit to_edge and from_edge
		tested[eidbase + 0] = true;
		tested[eidbase + 1] = true;
		discovered[cursor[i].a] = discovered[cursor[i].b] = true;

		vid2time[cursor[i].b] = pattern[i].time.b;
		if (i == 0) {
			vid2time[cursor[i].a] = pattern[i].time.a;
		}
		if (maxtoc < pattern[i].time.b) {
			maxtoc = pattern[i].time.b;
		}
	}

	bool valid_flg;
	DFSCode dcode;
	do {
		valid_flg = false;
		for (auto& i : Dice::shuffle(g.size())) {
			if (!discovered[i]) {
				continue;
			}
			for (auto& j : Dice::shuffle(g[i].size())) {
				Edge& added_edge = g[i][j];
				ID eidbase = added_edge.id - (added_edge.id % 2);
				if (discovered[added_edge.to]) {
					// backward
					// if (!tested[added_edge.id] and vid2time[i] > vid2time[added_edge.to]) {
					if (!tested[added_edge.id]) {
						dcode.labels = Triplet(-1, added_edge.labels.y, -1);
						dcode.time.set(vid2time[i], vid2time[added_edge.to]);
						pattern.push_back(dcode);
						cursor.push_back(VertexPair(i, added_edge.to, added_edge.id, added_edge.labels.y));
						valid_flg = true;
						// update tested
						tested[eidbase + 0] = true;
						tested[eidbase + 1] = true;
						break;
					}
				} else {
					// forward
					dcode.labels = Triplet(-1, added_edge.labels.y, added_edge.labels.z);
					dcode.time.set(vid2time[i], maxtoc+1);
					pattern.push_back(dcode);
					cursor.push_back(VertexPair(i, added_edge.to, added_edge.id, added_edge.labels.y));
					valid_flg = true;
					// update discovered & tested & vid2time & maxtoc
					discovered[added_edge.to] = true;
					tested[eidbase + 0] = true;
					tested[eidbase + 1] = true;
					vid2time[added_edge.to] = maxtoc+1;
					maxtoc++;
					break;
				}
			}
			if (valid_flg) {
				break;
			}
		}
	} while (!stop_condition(pattern, valid_flg));
	return pattern;
}

bool CLASS::stop_condition(const Pattern pattern, bool valid_flg) {
	// std::cout << "stop_condition: " << pattern << std::endl; // debug
	if (pattern.size() >= maxpat) {
		// std::cout << "maxpat" << std::endl;
		return true;
	}
	if (!valid_flg) {
		// std::cout << "no childs" << std::endl;
		return true;
	}
	if (Dice::p(1 - pow(setting.stopping_rate, pattern.size()))) {
		// std::cout << "probability" << std::endl;
		return true;
	}
	return false;
}

// only minDFS in DAG
// !!! minimum pattern not correspond EdgeTracer
bool Gspan::scanGspan(const Pattern& pattern) {
	// std::cout << "scanGspan: " << pattern << std::endl; // debug
	cache[pattern].scan = true;
	if (pattern.size() >= maxpat) {
		return false;
	}

	EdgeTracer cursor;
	DFSCode dcode;
	Pattern min_pat;

	map<Pattern, GraphToTracers> heap;
	set<set<ID>> ids_dic;
	for (auto x = cache[pattern].g2tracers.begin(); x != cache[pattern].g2tracers.end(); ++x) {
		int gid = x->first;
		Graph& g = db.gdata[gid];
		ids_dic = {};
		for (auto it = x->second.begin(); it != x->second.end(); ++it) {
			// an instance (a sequence of vertex pairs) as vector "vpair"
			cursor = *it;
			vector<bool> tested(g.num_of_edges);
			vector<bool> discovered(g.size());
			set<ID> ids = {};

			for (int i = cursor.size()-1; i >= 0; --i) {
				ID eidbase = cursor[i].id - (cursor[i].id % 2); // hit to_edge and from_edge
				tested[eidbase + 0] = true;
				tested[eidbase + 1] = true;
				discovered[cursor[i].a] = discovered[cursor[i].b] = true;
				ids.insert({eidbase, eidbase+1});
			}

			// make heap
			for (unsigned int i = 0; i < g.size(); i++) {
				if (!discovered[i]) {
					continue;
				}
				for (unsigned int j = 0; j < g[i].size(); j++) {
					const Edge& added_edge = g[i][j];
					ID eidbase = added_edge.id - (added_edge.id % 2);
					ids.insert({eidbase, eidbase+1});
					if (find(ids_dic.begin(), ids_dic.end(), ids) != ids_dic.end()) { // found
						if (!discovered[added_edge.to]) {
							ids.erase(eidbase);
							ids.erase(eidbase+1);
						}
					} else { // not found
						if (discovered[added_edge.to]) {
							// backward
							if (!tested[added_edge.id]) {
								ids_dic.insert(ids);
								cursor.push_back(VertexPair(i,added_edge.to,added_edge.id, added_edge.labels.y));
								auto pat_cur = is_min.convert(cursor, gid);
								heap[pat_cur.first][gid].push_back(pat_cur.second);
								cursor.pop_back();
							}
						} else {
							// forward
							ids_dic.insert(ids);
							cursor.push_back(VertexPair(i,added_edge.to,added_edge.id, added_edge.labels.y));
							auto pat_cur = is_min.convert(cursor, gid);
							heap[pat_cur.first][gid].push_back(pat_cur.second);
							cursor.pop_back();
						}
						ids.erase(eidbase);
						ids.erase(eidbase+1);
					}
				}
			}
		}
	}

	/*
	cout << "heap" << endl;
	for (auto x : heap) {
		cout << x.first << endl;
	}
	*/
	if (heap.empty()) {
		return false;
	}

	vector<Pattern> childs;
	bool child_flg = false;
	for (auto itr = heap.begin(); itr != heap.end(); itr++) {
		if (support(itr->second) < minsup) {
			continue;
		}
		cache[pattern].childs.push_back(itr->first);
		cache.insert({itr->first, CacheRecord(itr->second, childs)});
		child_flg = true;
	}
	return child_flg;
}

size_t CLASS::support(GraphToTracers& g2tracers) {
	// std::cout << "support" << std::endl; // debug
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

/*
   void CLASS::one_edge_report(GraphToTracers& g2tracers){
   std::cout << g2tracers.size() << " || " << endl;
   for(GraphToTracers::iterator x = g2tracers.begin(); x != g2tracers.end(); ++x){
   std::cout << x->first << ":" << endl;
   for(Tracers::iterator y = x->second.begin(); y != x->second.end(); ++y){
   std::cout << "(" << y->vpair.a << " " << y->vpair.b << ") ";
   }
   std::cout << std::endl;
   }
   std::cout << std::endl;
   }
   */
