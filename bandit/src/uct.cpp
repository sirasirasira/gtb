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
	path.push_back(pattern);
	if (cache[pattern].terminal) {
		if (cache[pattern].count >= setting.threshold) {
			expansion(pattern);
		}
	} else {
		if (cache[pattern].childs.size() != 0) {
			Pattern best_child;
			double max_ucb = 0;
			for (auto& c : cache[pattern].childs) {
				if (cache[c].ucb > max_ucb) {
					best_child = c;
					max_ucb = cache[c].ucb;
				}
			}
			selection(best_child);
		}
	}
}

void CLASS::expansion(const Pattern & pattern) {
	cache[pattern].terminal = false;
	if (db.gspan.scanGspan(pattern)) {
		path.push_back(cache[pattern].childs[0]);
	}
}

Pattern CLASS::simulation(const Pattern& pattern) {
	auto& g2tracers = cache[pattern].g2tracers;
	ID gid = Dice::id(g2tracers.size());
	auto& tracers = g2tracers[gid];
	auto& tracer = tracers[Dice::id(tracers.size())];
	
	tuple<Pattern, EdgeTracer, ID> pat_et_gid = make_tuple(pattern, tracer, gid);
	do {
		pat_et_gid = db.gspan.oneEdgeSimulation(pat_et_gid);
	} while (!stop_condition(pat_et_gid));
	return std::get<0>(pat_et_gid);
}

bool CLASS::stop_condition(const tuple<Pattern, EdgeTracer, ID>& pat_et_gid) {
	if (std::get<0>(pat_et_gid).size() >= setting.maxpat) {
		return true;
	}
	if (std::get<1>(pat_et_gid).predec == nullptr) {
		return true;
	}
	//TODO
	if (Dice::p(0.1)) {
		return true;
	}
	return false;
}

void CLASS::backpropagation(double score) {
	for (int i = path.size() - 1; i > 0; i--) {
		cache[path[i]].count++;
		cache[path[i]].sum_score += score;
		double count = cache[path[i]].count;
		double ave = cache[path[i]].sum_score / count;
		double p_count = cache[path[i-1]].count + 1;
		cache[path[i]].ucb = ave + setting.c * sqrt(2 * log(p_count) / count);
	}
	cache[root].count++;
}
