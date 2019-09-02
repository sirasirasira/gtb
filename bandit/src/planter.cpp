#include "Planter.h"
#define CLASS Planter

extern Setting setting;
#include "Database.h"
extern Database db;

#include "Calculator.h"

// each value of additive_ys must be overwritten only once
const vector<double>& CLASS::run(const vector<ID>& train_targets, const vector<ID>& test_targets) {
	// std::cout << "debug Planter run" << std::endl; // debug
	cout << "BEGIN_tree" << endl;
	grow(train_targets, test_targets, 0);
	cout << "END_tree" << endl;
	return additive_ys;
}

void CLASS::grow(const vector<ID>& train_targets, const vector<ID>& test_targets, size_t depth) {
	// std::cout << "debug grow, depth " << depth << endl; // debug
	// Debug::IDs(targets, "targets "); // debug
	if (checkLeafConditions(train_targets, depth)) {
		makeLeaf(train_targets, test_targets, depth);
		return;
	}
	vector<ID> posi_train_targets = db.spliter.run(train_targets);
	if (db.spliter.valid() == false) {
		makeLeaf(train_targets, test_targets, depth);
		return;
	}
	const auto best_pattern = db.spliter.getBestPattern();
	db.gspan.updateFeatureImportance(best_pattern, db.spliter.getImportance());
	vector<ID> posi_test_targets = db.finder.run(best_pattern, test_targets);
	grow(posi_train_targets, posi_test_targets, depth + 1);
	cout << string(depth, '-') << "* " << best_pattern << endl;
	vector<ID> nega_train_targets = Calculator::setDiff(train_targets, posi_train_targets);
	vector<ID> nega_test_targets = Calculator::setDiff(test_targets, posi_test_targets);
	grow(nega_train_targets, nega_test_targets, depth + 1);
}

bool CLASS::checkLeafConditions(const vector<ID>& targets, size_t depth) {
	// if (targets.size() <= db.gspan.minsup * 2) return true;
	if (targets.size() < 2) return true;
	if (depth >= setting.max_depth) return true;
	return false;
}

// @change additive_ys
void CLASS::makeLeaf(const vector<ID>& train_targets, const vector<ID>& test_targets, size_t depth) {
	// std::cout << "debug makeLeaf" << std::endl; // debug
	// Debug::IDs(targets); // debug
	double sum = 0;
	size_t count = train_targets.size();
	for (ID id : train_targets) {
		sum += db.ys[id];
	}
	double mean = sum / (double) count;
	updateAdditiveYs(train_targets, test_targets, mean);
	// Debug::ys(additive_ys); // debug
	cout << string(depth, '-') << "* " << "[" << count << "] " << mean << endl;
}

void CLASS::updateAdditiveYs(const vector<ID>& train_targets, const vector<ID>& test_targets, double value) {
	for (auto id : train_targets) {
		additive_ys[id] = value;
	}
	for (auto id : test_targets) {
		additive_ys[id] = value;
	}
}
