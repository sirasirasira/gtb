#include "UCT.h"
#define CLASS UCT

extern Setting setting;
#include "Database.h"
#include "Calculator.h" 
extern Database db;

void CLASS::run(const vector<ID>& _targets) {
	targets = _targets;
	db.gspan.clearUCB();
	Pattern pattern;
	vector<ID> posi;
	double score;

	for (unsigned int i = 0; i < setting.iteration; i++) {
		path = {};
		selection(root);
		expansion();
		if (update()) {
			//TODO
			backpropagation(-DBL_MAX);
			continue;
		}
		pattern = path[path.size()-1];
		pattern = simulation(pattern);
		if (cache.find(pattern) == cache.end()) {
			posi = db.finder.run(pattern, targets);
		} else {
			posi = db.gspan.getPosiIds(cache[pattern].g2tracers);
		}
		score = Calculator::score(db.ys, targets, posi);
		backpropagation(score);
	}
}

void CLASS::selection(const Pattern& pattern) {
	// std::cout << "selection : " << pattern << std::endl; // debug
	path.push_back(pattern);
	if (!cache[pattern].terminal) {
		Pattern best_child;
		double ucb;
		double max_ucb = -DBL_MAX;
		for (auto& c : cache[pattern].childs) {
			if (cache[c].prune) {
				continue;
			} else {
				if (cache[c].count == 0) {
					ucb = DBL_MAX;
				} else {
					ucb = (cache[c].sum_score / cache[c].count) + (setting.c * sqrt(2 * log(cache[pattern].count / cache[c].count)));
				}
				if (ucb > max_ucb) {
					best_child = c;
					max_ucb = ucb;
				}
			}
		}
		selection(best_child);
	}
}

void CLASS::expansion() {
	const Pattern pattern = path[path.size()-1];
	// std::cout << "expansion: " << pattern << endl; // debug
	if (cache[pattern].count >= setting.threshold) {
		if (!cache[pattern].scan) {
			if (db.gspan.scanGspan(pattern)) {
				cache[pattern].terminal = false;
				path.push_back(cache[pattern].childs[0]);
			}
		} else {
			if (cache[pattern].childs.size()) {
				cache[pattern].terminal = false;
				path.push_back(cache[pattern].childs[0]);
			}
		}
	}
}

bool CLASS::update() {
	vector<ID> posi;
	double score;
	double bound;
	bool pruning = false;
	int pruning_index;
	for (unsigned int i = 0; i < path.size(); i++) {
		if (pruning) {
			cache[path[i]].prune = true;
		} else {
			if (cache[path[i]].count == 0) {
				posi = db.gspan.getPosiIds(cache[path[i]].g2tracers);
				score = Calculator::score(db.ys, targets, posi);
				db.spliter.update(path[i], score);
				bound = Calculator::bound(db.ys, targets, posi);
				cache[path[i]].bound = bound;
				if (db.spliter.isBounded(bound)){
					cache[path[i]].prune = true;
					pruning = true;
					pruning_index = i;
				}
			}
		}
	}
	if (pruning) {
		path.resize(pruning_index);
	}
	return pruning;
}

Pattern CLASS::simulation(const Pattern& pattern) {
	// cout << "simulation: " << pattern << endl;
	auto& g2tracers = cache[pattern].g2tracers;
	// gid is in g2tracers and in targets
	vector<ID> keys(g2tracers.size());
	size_t i = 0;
	for (auto& x : g2tracers) {
		keys[i] = x.first;
		i++;
	}
	auto gids = Calculator::setIntersec(targets, keys);
	ID gid = gids[Dice::id(gids.size())];
	auto& tracers = g2tracers[gid];
	auto& tracer = tracers[Dice::id(tracers.size())];
	
	return db.gspan.EdgeSimulation(pattern, tracer, gid);
}
 
void CLASS::backpropagation(double score) {
	// cout << "backpropagation: " << score << endl;
	for (int i = path.size() - 1; i > 0; i--) {
		cache[path[i]].count++;
		cache[path[i]].sum_score -= score; // score is mse (lower is better)
	}
	cache[root].count++;
}
