#include "Planter.h"
#define CLASS Planter

extern Setting setting;
#include "Database.h"
extern Database db;

#include "Calculator.h"

// each value of additive_ys must be overwritten only once
const vector<double>& CLASS::run(const vector<ID>& targets) {
	// std::cout << "debug Planter run" << std::endl; // debug
	cout << "BEGIN_tree" << endl;
	grow(targets, 0);
	cout << "END_tree" << endl;
	return additive_ys;
}

void CLASS::grow(vector<ID> targets, size_t depth) {
	// std::cout << "debug grow, depth " << depth << endl; // debug
	// Debug::IDs(targets, "targets "); // debug
	if (checkLeafConditions(targets, depth)) {
		makeLeaf(targets, depth);
		return;
	}
	vector<ID> posi_targets = db.spliter.run(targets);
	if (db.spliter.valid() == false) {
		makeLeaf(targets, depth);
		return;
	}
	const auto best_pattern = db.spliter.getBestPattern();
	db.gspan.updateFeatureImportance(best_pattern, db.spliter.getImportance());
	grow(posi_targets, depth + 1);
	cout << string(depth, '-') << "* " << best_pattern << endl;
	vector<ID> nega_targets = Calculator::setDiff(targets, posi_targets);
	grow(nega_targets, depth + 1);
}

bool CLASS::checkLeafConditions(const vector<ID>& targets, size_t depth) {
	// if (targets.size() <= db.gspan.minsup * 2) return true;
	if (targets.size() < 2) return true;
	if (depth >= setting.max_depth) return true;
	return false;
}

// @change additive_ys
void CLASS::makeLeaf(const vector<ID>& targets, size_t depth) {
	// std::cout << "debug makeLeaf" << std::endl; // debug
	// Debug::IDs(targets); // debug
	double sum = 0;
	size_t count = 0;
	for (ID id : targets) {
		if ((id <= db.gdata.getLastTrainID()) == false) {
			break;
		}
		count++;
		sum += db.ys[id];
	}
	double mean = sum / (double) count;
	updateAdditiveYs(targets, mean);
	// Debug::ys(additive_ys); // debug
	cout << string(depth, '-') << "* " << "[" << count << "] " << mean << endl;
}

void CLASS::updateAdditiveYs(const vector<ID>& targets, double mean) {
	for (auto id : targets) {
		additive_ys[id] = mean;
	}
}

