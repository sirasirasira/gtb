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
			std::cout << "root child : " << c << std::endl; // debug
			path = {root, c};
			const auto sim_pat = simulation(c);
			std::cout << sim_pat << std::endl; // debug
			std::cout << "finder b" << std::endl; // debug
			posi = db.finder.run(sim_pat, targets);
			std::cout << "finder e" << std::endl; // debug
			score = Calculator::score(db.ys, targets, posi);
			std::cout << "score : " << score << std::endl; // debug
			backpropagation(score);
		}
	}

	for (unsigned int i = 0; i < setting.iteration - (cache[root].childs.size() * setting.threshold); i++) {
		path = {};
		selection(root);
		if (update()) {
			//TODO
			backpropagation(-DBL_MAX);
			continue;
		}
		pattern = path[path.size()-1];
		if (cache[pattern].terminal) {
			pattern = simulation(pattern);
			std::cout << "finder start" << std::endl; // debug
			posi = db.finder.run(pattern, targets);
			std::cout << "finder end" << std::endl; // debug
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
	std::cout << "selection : " << pattern << std::endl; // debug
	std::cout << "selection : " << cache[pattern].terminal << std::endl; // debug
	path.push_back(pattern);
	if (cache[pattern].terminal == true) {
	std::cout << "selection2 : " << pattern << std::endl; // debug
		if (cache[pattern].count >= setting.threshold-1) {
			expansion(pattern);
		}
	} else {
	std::cout << "selection3 : " << pattern << std::endl; // debug
		if (cache[pattern].childs.size() != 0) {
			Pattern best_child;
			double max_ucb = -DBL_MAX;
			for (auto& c : cache[pattern].childs) {
				if (cache[c].ucb > max_ucb) {
					best_child = c;
					max_ucb = cache[c].ucb;
				}
			}
			std::cout << "best child : " << best_child << std::endl; // debug
			selection(best_child);
		}
	}
}

void CLASS::expansion(const Pattern & pattern) {
	std::cout << "expansion : " << pattern << std::endl; // debug
	cache[pattern].terminal = false;
	if (!cache[pattern].scan) {
		std::cout << "scan : " << pattern << std::endl; // debug
		db.gspan.scanGspan(pattern);
	}
	if (cache[pattern].childs.size() != 0) {
		path.push_back(cache[pattern].childs[0]);
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
	ID gid = Dice::id(g2tracers.size());
	auto& tracers = g2tracers[gid];
	auto& tracer = tracers[Dice::id(tracers.size())];
	 std::cout << "gid:" << gid << std::endl; // debug
	 std::cout << &tracer << std::endl;
	 std::cout <<  "a:" << tracer.vpair.a << std::endl;
	 std::cout <<  "b:" << tracer.vpair.b << std::endl;
	 std::cout << "id:" << tracer.vpair.id << std::endl; // debug
	
	return db.gspan.EdgeSimulation(pattern, tracer, gid);
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
