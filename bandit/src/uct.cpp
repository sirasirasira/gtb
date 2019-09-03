#include "UCT.h"
#define CLASS UCT

extern Setting setting;
#include "Database.h"
#include "Calculator.h" 
extern Database db;

void CLASS::run(const vector<ID>& _targets) {
	targets = _targets;
	db.gspan.clearUCB();
	cache = db.gspan.getCache();
	root = db.gspan.getRoot();
	Pattern pattern;
	vector<ID> posi;
	double score;

	for (const auto& c : cache[root].childs) { // minimum itaration = one edge graphs size
		for (unsigned int i = 0; i < setting.threshold; i++) {
			path = {root, c};
			pattern = simulation(c);
			posi = db.finder.run(pattern, targets);
			score = Calculator::score(db.ys, targets, posi);
			backpropagation(score);
		}
	}

	for (unsigned int i = 0; i < setting.iteration - (cache[root].childs.size() * setting.threshold); i++) {
		path = {};
		selection(root);
		cout << "path" << endl;
		for (auto x : path) {
			cout << x << endl;
		}
		expansion();
		if (update()) {
			//TODO
			backpropagation(-DBL_MAX);
			continue;
		}
		pattern = path[path.size()-1];
		if (cache[pattern].terminal) {
			pattern = simulation(pattern);
			posi = db.finder.run(pattern, targets);
			score = Calculator::score(db.ys, targets, posi);
			backpropagation(score);
		} else {
			posi = db.gspan.getPosiIds(cache[pattern].g2tracers);
			score = Calculator::score(db.ys, targets, posi);
			backpropagation(score);
		}
	}
}

void CLASS::selection(const Pattern& pattern) {
	// std::cout << "selection : " << pattern << std::endl; // debug
	path.push_back(pattern);
	if (!cache[pattern].terminal) {
		Pattern best_child;
		double max_ucb = -DBL_MAX;
		for (auto& c : cache[pattern].childs) {
			if (cache[c].ucb > max_ucb) {
				best_child = c;
				max_ucb = cache[c].ucb;
			}
		}
		selection(best_child);
	}
}

void CLASS::expansion() {
	const Pattern pattern = path[path.size()-1];
	 std::cout << "expansion : " << pattern << endl; // debug
	if (cache[pattern].count >= setting.threshold) {
		if (!cache[pattern].scan) {
			if (db.gspan.scanGspan(pattern)) {
				cache[pattern].terminal = false;
	cout << "debug" << endl;
	cout << "childs size: " <<cache[pattern].childs.size() << endl;
	for (auto x : cache[pattern].childs) {
		cout << x << endl;
	}
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
			cache[path[i]].ucb = -DBL_MAX;
		} else {
			if (cache[path[i]].count == 0) {
				posi = db.gspan.getPosiIds(cache[path[i]].g2tracers);
				score = Calculator::score(db.ys, targets, posi);
				db.spliter.update(path[i], score);
				bound = Calculator::bound(db.ys, targets, posi);
				cache[path[i]].bound = bound;
				if (db.spliter.isBounded(bound)){
					cache[path[i]].ucb = -DBL_MAX;
					pruning = true;
					pruning_index = i;
				}
			}
		}
	}
	if (pruning) {
		path.resize(pruning_index+1);
	}
	return pruning;
}

Pattern CLASS::simulation(const Pattern& pattern) {
	auto& g2tracers = cache[pattern].g2tracers;
	auto itr = g2tracers.begin();
	std::advance(itr, Dice::id(g2tracers.size()));
	auto& tracers = itr->second;
	auto& tracer = tracers[Dice::id(tracers.size())];
	
	return db.gspan.EdgeSimulation(pattern, tracer, itr->first);
}
 
void CLASS::backpropagation(double score) {
	for (int i = path.size() - 1; i > 0; i--) {
		cache[path[i]].count++;
		cache[path[i]].sum_score -= score; // score is mse (lower is better)
		double count = cache[path[i]].count;
		double ave = cache[path[i]].sum_score / count;
		double p_count = cache[path[i-1]].count + 1;
		cache[path[i]].ucb = ave + setting.c * sqrt(2 * log(p_count) / count);
	}
	cache[root].count++;
}
