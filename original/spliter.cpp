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
	initMinScore();
	// Debug::IDs(targets); // debug
	db.gspan.setSpliterPtr(this);
	db.gspan.minsup = setting.minsup;
	db.gspan.maxpat = 3;
	db.gspan.run();
	db.gspan.maxpat = setting.maxpat;
	std::cout << "prepare cache size: " << db.gspan.getCache().size() << std::endl;
}

vector<ID> CLASS::run(const vector<ID>& _targets) {
	// std::cout << "debug spliter run" << std::endl; // debug
	targets = _targets;
	best_pattern = {};
	initMinScore();
	if (min_score < 0) {
		goto G_INVALID;
	}
	searchCashe();
	searchEnum();
	// std::cout << "debug best_pattern " << best_pattern << std::endl; // debug
	if (best_pattern.size() == 0) {
G_INVALID:
		valid_flg = false;
		return {};
	}
	// std::cout << "debug parent_score " << parent_score << " min_score " << min_score << std::endl; // debug
	valid_flg = true;
	vector<ID> posi = db.gspan.getPosiIds(db.gspan.getCache().at(best_pattern).g2tracers);
	return Calculator::setIntersec(targets, posi);
}

void CLASS::searchCashe() {
	const auto& cache = db.gspan.getCache();
	for (auto itr = cache.begin(); itr != cache.end(); itr++) {
		const auto& pattern = itr->first;
		const auto& tracers = itr->second.g2tracers;
		double score = Calculator::score(db.ys, targets, db.gspan.getPosiIds(tracers));
		if (score < min_score) {
			min_score = score;
			best_pattern = pattern;
		}
	}
}

bool CLASS::isFrontier(const Pattern& pattern, const Pattern& bef_pattern) {
	if (bef_pattern.size() - 1 == pattern.size()) {
		return false;
	}
	return true;
}

void CLASS::searchEnum() {
	const auto& cache = db.gspan.getCache();
	vector<Pattern> frontier_set;
	auto bef_pattern = cache.rbegin()->first;
	frontier_set.push_back(bef_pattern);
	auto bitr = cache.rbegin();
	bitr++;
	for (; bitr != cache.rend(); bitr++) {
		auto pattern = bitr->first;
		if (isFrontier(pattern, bef_pattern)) {
			frontier_set.push_back(pattern);
		}
		bef_pattern = pattern;
	}
	
	//std::shuffle(frontier_set.begin(), frontier_set.end(), Dice::mt);
	std::reverse(frontier_set.begin(), frontier_set.end());
	for (auto ptr : frontier_set) {
		db.gspan.run(ptr);
	}
}

void CLASS::update(Pattern pattern, vector<ID> posi) {
	if (Dice::p(setting.feature_used) == false) {
		return;
	}
	double score = Calculator::score(db.ys, targets, posi);
	if (score < min_score ) { // old pattern may be used (this func is called from gspan)
		min_score = score;
		best_pattern = pattern;
	}
}

bool CLASS::isBounded(vector<ID> posi) {
	double min_bound = Calculator::bound(db.ys, targets, posi);
	if (min_score <= min_bound) {
		return true;
	} else {
		return false;
	}
}

