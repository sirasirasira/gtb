#include "Spliter.h"
#define CLASS Spliter

extern Setting setting;
#include "Database.h"
extern Database db;

#include "Calculator.h" 

void CLASS::initMinScore() {
	parent_score = Calculator::imp(db.ys, Calculator::trainOnly(targets));
	min_score = parent_score - setting.needed_impurity_decrease - std::numeric_limits<double>::epsilon();
}

void CLASS::prepare(const vector<ID>& _targets) {
	// std::cout << "debug Spliter prepare" << std::endl; // debug
	targets = _targets;
	// Debug::IDs(targets); // debug
	db.gspan.setSpliterPtr(this);
	db.gspan.makeRoot(targets);
	//std::cout << "prepare cache size: " << db.gspan.getCache().size() << std::endl;
}

vector<ID> CLASS::run(const vector<ID>& _targets) {
	// std::cout << "spliter run" << std::endl; // debug
	targets = _targets;
	best_pattern = {};
	initMinScore();
	if (min_score < 0) {
		valid_flg = false;
		return {};
	}
	db.uct.run(targets);
	// std::cout << "debug best_pattern " << best_pattern << std::endl; // debug
	if (best_pattern.size() == 0) {
		valid_flg = false;
		return {};
	}
	// std::cout << "debug parent_score " << parent_score << " min_score " << min_score << std::endl; // debug
	valid_flg = true;
	vector<ID> posi = db.gspan.getPosiIds(db.gspan.getCache().at(best_pattern).g2tracers);
	return Calculator::setIntersec(targets, posi);
}

void CLASS::update(const Pattern& pattern, double score) {
	if (score < min_score ) { // old pattern may be used (this func is called from gspan)
		min_score = score;
		best_pattern = pattern;
	}
}

bool CLASS::isBounded(double min_bound) {
	// cout << "isBounded" << endl;
	// cout << "min_score: " << min_score << " bound: " << min_bound << endl;
	if (min_score <= min_bound) {
		return true;
	} else {
		return false;
	}
}

