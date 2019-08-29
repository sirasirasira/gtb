#include "Gspan.h"
#define CLASS Gspan

#include "Database.h"
extern Database db;

void CLASS::makeRoot() {
	Pattern pattern;
	DFSCode dcode;
	dcode.labels = Triplet(-1, -1, -1);
	dcode.time.set(0, 0);
	Pattern _root{dcode};
	root = _root;
	vector<Pattern> childs; // empty vec
	cache.insert({root, CacheRecord(childs)});
	cache[root].scan = true;

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
Pattern CLASS::EdgeSimulation(const Pattern& _pattern, EdgeTracer& _tracer, ID gid){
	 std::cout << "debug EdgeSimulation" << std::endl; // debug
	std::cout << _pattern << std::endl; // debug
	std::cout << &_tracer << std::endl; // debug

	Pattern pattern = _pattern;
	EdgeTracer *tracer = &(_tracer);
	vector<VertexPair> vpairs(pattern.size());
	EdgeTracer cursor;

	Graph& g = db.gdata[gid];
	size_t maxtoc = 0;
	vector<bool> tested(g.num_of_edges); // edge
	vector<bool> discovered(g.size()); // node
	vector<int> vid2time(g.size(), -1);

	for (int i = vpairs.size()-1; i >= 0; --i, tracer = tracer->predec) {
	 std::cout << "111111111111" << std::endl; // debug
	 std::cout <<  "a:" << tracer->vpair.a << std::endl;
	 std::cout <<  "b:" << tracer->vpair.b << std::endl;
	 std::cout << "id:" << tracer->vpair.id << std::endl; // debug
		vpairs[i] = tracer->vpair;
	 std::cout << "222222222222" << std::endl; // debug
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
				if (discovered[added_edge.to]) {
					// backward
					if (!tested[added_edge.id] and vid2time[i] > vid2time[added_edge.to]) {
						dcode.labels = Triplet(-1, added_edge.labels.y, -1);
						dcode.time.set(vid2time[i], vid2time[added_edge.to]);
						pattern.push_back(dcode);
						cursor.set(i,added_edge.to,added_edge.id,&(*tracer));
						valid_flg = true;
						// update discovered & tested & maxtoc
						ID vidbase = added_edge.id - (added_edge.id % 2);
						tested[vidbase + 0] = true;
						tested[vidbase + 1] = true;
						break;
					}
				} else {
					// forward
					dcode.labels = Triplet(-1, added_edge.labels.y, added_edge.labels.z);
					dcode.time.set(vid2time[i], maxtoc+1);
					pattern.push_back(dcode);
					cursor.set(i,added_edge.to,added_edge.id,&(*tracer));
					valid_flg = true;
					// update discovered & tested & maxtoc
					discovered[added_edge.to] = true;
					ID vidbase = added_edge.id - (added_edge.id % 2);
					tested[vidbase + 0] = true;
					tested[vidbase + 1] = true;
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
	// std::cout << "stop_condition" << std::endl; // debug
	 std::cout << pattern << std::endl;
	if (pattern.size() > maxpat) {
		// std::cout << "maxpat" << std::endl;
		return true;
	}
	if (!valid_flg) {
		// std::cout << "no childs" << std::endl;
		return true;
	}
	//TODO
	if (Dice::p(0.1)) {
		// std::cout << "probability" << std::endl;
		return true;
	}
	return false;
}

// minDFS only DAG
// !!! minimum pattern not correspond EdgeTracer
void Gspan::scanGspan(const Pattern& _pattern) {
	// std::cout << "debug scanGspan" << std::endl; // debug
	Pattern pattern = _pattern;
	cache[pattern].scan = true;
	if (pattern.size() >= maxpat) {
		return;
	}
	GraphToTracers g2tracers = cache[pattern].g2tracers;
	size_t maxtoc = 0;
	for (auto it = pattern.rbegin(); it != pattern.rend(); it++) {
		if (it->time.a < it->time.b) {
			maxtoc = it->time.b;
			break;
		}
	}
	vector<VertexPair> vpairs(pattern.size());
	EdgeTracer *tracer;
	Pair pkey;
	EdgeTracer cursor;
	DFSCode dcode;
	Pattern min_pat;

	map<Pattern, GraphToTracers> heap;
	for (auto x = g2tracers.begin(); x != g2tracers.end(); ++x) {
		int gid = x->first;
		for (auto it = x->second.begin(); it != x->second.end(); ++it) {
			// an instance (a sequence of vertex pairs) as vector "vpair"
			tracer = &(*it);

			Graph& g = db.gdata[gid];
			vector<bool> tested(g.num_of_edges);
			vector<bool> discovered(g.size());
			vector<int> vid2time(g.size(), -1);

			for (int i = vpairs.size()-1; i >= 0; --i, tracer = tracer->predec) {
				vpairs[i] = tracer->vpair;
				auto vidbase = vpairs[i].id - (vpairs[i].id % 2); // hit to_edge and from_edge
				tested[vidbase + 0] = true;
				tested[vidbase + 1] = true;
				discovered[vpairs[i].a] = discovered[vpairs[i].b] = true;

				vid2time[vpairs[i].b] = pattern[i].time.b;
				if (i == 0) {
					vid2time[vpairs[i].a] = pattern[i].time.a;
				}
			}

			// make heap
			for (unsigned int i = 0; i < g.size(); i++) {
				if (!discovered[i]) {
					continue;
				}
				for (unsigned int j = 0; j < g[i].size(); j++) {
					Edge& added_edge = g[i][j];
					if (discovered[added_edge.to]) {
						// backward
						if (!tested[added_edge.id] and vid2time[i] > vid2time[added_edge.to]) {
							dcode.labels = Triplet(-1, added_edge.labels.y, -1);
							dcode.time.set(vid2time[i], vid2time[added_edge.to]);
							pattern.push_back(dcode);
							min_pat = is_min.convert(pattern);
							pkey.set(added_edge.labels.y,vid2time[added_edge.to]);
							cursor.set(i,added_edge.to,added_edge.id,&(*it));
							heap[min_pat][gid].push_back(cursor);
							pattern.pop_back();
						}
					} else {
						// forward
						dcode.labels = Triplet(-1, added_edge.labels.y, added_edge.labels.z);
						dcode.time.set(vid2time[i], maxtoc);
						pattern.push_back(dcode);
						min_pat = is_min.convert(pattern);
						pkey.set(added_edge.labels.y,added_edge.labels.z);
						cursor.set(i,added_edge.to,added_edge.id,&(*it));
						heap[min_pat][gid].push_back(cursor);
						pattern.pop_back();
					}
				}
			}
		}
	}

	vector<Pattern> childs;
	for (auto itr = heap.begin(); itr != heap.end(); itr++) {
		if (support(itr->second) < minsup) {
			continue;
		}
		cache[pattern].childs.push_back(itr->first);
		cache.insert({itr->first, CacheRecord(itr->second, childs)});
	}
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
