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
	root = db.gspan.getroot();
	Pattern pattern;

	for (const auto& c : cache[root].childs) { // minimum itaration = one edge graphs size
		path = {c};
		pattern = simulation(c);
		backpropagation(pattern);
	}

	for (unsigned int i = 0; i < setting.iteration - cache[root].childs.size(); i++) {
		path = {};
		/*
		pattern = selection();
		expansion();
		simulation();
		backpropagation();
		*/
	}
}

Pattern CLASS::selection() {
	return root;
}

void CLASS::expansion() {
}

void CLASS::update(){
}

Pattern CLASS::simulation(const Pattern& pattern) {
	// random select edgetracer
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

void CLASS::backpropagation(const Pattern& pattern) {
	vector<ID> posi = db.finder.run(pattern, targets);
	double score = Calculator::score(db.ys, targets, posi);
	// update
	for (int i = path.size(); i >= 0; i--) {
		cache[path[i]].count++;
		cache[path[i]].sum_score += score;
		if (cache[path[i]].childs.size() * setting.c_threshold <= cache[path[i]].count) {
			cache[path[i]].terminal = false;
		}
		double count = cache[path[i]].count;
		double ave = cache[path[i]].sum_score / count;
		double p_count;
		if (i == 0) {
			p_count = cache[root].count + 1;
		} else {
			p_count = cache[path[i-1]].count + 1;
		}
		cache[path[i]].ucb = ave + setting.c * sqrt(2 * log(p_count) / count);
	}
	cache[root].count++;
}
