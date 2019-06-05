#include "Finder.h"
#define CLASS Finder

#include "Database.h"
extern Database db;

vector<ID> CLASS::run(const Pattern& best_pattern, const vector<ID>& targets) {
	//cout << "finder run" << endl;
	pattern = best_pattern;

	vector<ID> posi;
	for (auto gid : targets) {
		if (gcheck(gid)) {
			posi.push_back(gid);
		}
	}
	return posi;
}

bool CLASS::gcheck(const ID gid) {
	//cout << "gcheck id: " << gid << endl;
	exist = false;
	g = db.gdata[gid];

	vector<ID> vheap;
	for (ID i = 0; i < g.size(); i++) {
		if (g.label[i] == pattern[0].labels.x) {
			vheap.push_back(i);
		}
	}
	for (auto root : vheap) {
		vids.resize(1);
		vids[0] = root;
		gtraverse(0);
	}
	if (exist) {
		return true;
	} else {
		return false;
	}
}

void CLASS::gtraverse(unsigned int cur) {
	if (pattern.size() <= cur) {
		exist = true;
	}
	if (exist) {
		return;
	}

	const int from = pattern[cur].time.a;
	const int to = pattern[cur].time.b;
	if (to < from) {
		// backward
		for (auto e : g[vids[from]]) {
			if (e.to == vids[to]) {
				//TODO
				gtraverse(cur+1);
			}
		}
	} else {
		// forward
		int elabel = pattern[cur].labels.y;
		int vlabel = pattern[cur].labels.z;
		for (auto e : g[vids[from]]) {
			if (e.labels.y == elabel and e.labels.z == vlabel and ndiscovered(e.to)) {
				vids.push_back(e.to);
				gtraverse(cur+1);
				vids.pop_back();
			}
		}
	}
}

bool CLASS::ndiscovered(const int id) {
	if (find(vids.begin(), vids.end(), id) == vids.end()) {
		return true;
	} else {
		return false;
	}
}
