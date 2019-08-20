#include "UCT.h"
#define CLASS UCT

extern Setting setting;
#include "Database.h"
extern Database db;

void CLASS::run(const vector<ID>& _targets) {
	target = _targets;
	cache = db.gspan.getCache();
	e1patterns = db.gspan.gete1Patterns();
	for (int i = 0; i < db.setting.iteration; i++) {
		if (i < e1patterns.size()) {
			pattern = e1patterns[i];
			vector<ID> posi = db.gspan.getPosiIds(cache[pattern].g2tracer);
			db.spliter.update(pattern, posi);
			path = {};
			simulation(pattern);
			backpropagation(path);
		} else {
			pattern = selection();
			expansion();
			simulation();
			backpropagation();
		}
	}
}
