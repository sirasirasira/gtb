#include "UCT.h"
#define CLASS UCT

extern Setting setting;
#include "Database.h"
extern Database db;

void CLASS::run(const vector<ID>& _targets) {
	targets = _targets;
	db.gspan.clearUCB();
	cache = db.gspan.getCache();
	root = db.gspan.getroot();

	for (const auto& c : cache[root].childs) { // minimum itaration = one edge graphs size
		path = {c};
		simulation(c);
		backpropagation();
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

Pattern CLASS::simulation(const Pattern& _pattern) {
	//TODO random select _edgetracer
	EdgeTracer _edgetracer;
	
	pair<Pattern, EdgeTracer> pat_et = make_pair(_pattern, _edgetracer);
	do {
		pat_et = db.gspan.oneEdgeSimulation(pat_et);
	} while (!stop_condition(pat_et.first));
	return pat_et.first;
}

void CLASS::backpropagation() {
}

bool CLASS::stop_condition(const Pattern& pattern) {
	if (pattern.size() >= setting.maxpat) {
		return true;
	}
	//TODO
	if (Dice::p(0.1)) {
		return true;
	}
	return false;
}
