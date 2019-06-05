#include "TreeEnsemble.h"
#define CLASS TreeEnsemble

extern Setting setting;
#include "Database.h"
extern Database db;

#include "Calculator.h" 

// --- Random Forest ---

void CLASS::runRandomForest() {
	// std::cout << "debug runRandomForest" << std::endl; // debug
	db.y_predictions.resize(db.raw_ys.size(), 0);
	for (tree_count = 1; tree_count <= setting.num_of_trees; tree_count++) {
		makeTargets();
		const vector<double>& additive_ys = db.planter.run(targets);
		for (ID id = 0; id < db.gdata.size(); id++) {
			db.y_predictions[id] = (tree_count * db.y_predictions[id] + additive_ys[id]) / (double) tree_count;
		}
		report();
	}
}

// --- Gradient Boosting ---

void CLASS::runGradientBoosting() {
	// std::cout << "debug runGradientBoosting" << std::endl; // debug
	plantFirst();
	for (tree_count = 1; tree_count <= setting.num_of_trees; tree_count++) {
		calcResidualErrors();
		makeTargets();
		const vector<double>& additive_ys = db.planter.run(targets);
		for (ID id = 0; id < db.gdata.size(); id++) {
			db.y_predictions[id] += setting.shrinkage * additive_ys[id];
		}
		report();
	}
}

void CLASS::plantFirst() {
	// std::cout << "debug plantFirst" << std::endl; // debug
	double sum = 0;
	for (ID id = db.gdata.getFirstTrainID(); id <= db.gdata.getLastTrainID(); id++) {
		sum += db.raw_ys[id];
	}
	double mean = sum / (double) db.gdata.num_train;
	for (ID id = 0; id < db.gdata.size(); id++) {
		db.y_predictions[id] = mean;
	}
}

// @change db.ys
void CLASS::calcResidualErrors() {
	// std::cout << "debug calcResidualErrors" << std::endl; // debug
	for (ID id = 0; id < db.ys.size(); id++) {
		db.ys[id] = Calculator::calcResidualErr(db.raw_ys[id], db.y_predictions[id]);
	}
	// Debug::ys(db.ys, "residual_errors"); // debug
}

// --- Common ---

// @change targets
void CLASS::makeTargets() {
	// std::cout << "debug makeTargets" << std::endl; // debug
	ID id;
	targets.clear();
	for (id = db.gdata.getFirstTrainID(); id <= db.gdata.getLastTrainID(); id++) {
		if (Dice::p(setting.data_used)) {
			targets.push_back(id);
		}
	}
	for (id = db.gdata.getFirstTestID(); id <= db.gdata.getLastTestID(); id++) {
		targets.push_back(id);
	}
	// Debug::IDs(targets); // debug
}

void CLASS::report() {
	reportFeatureImportance();
	db.eva.run(tree_count); // acc, auc
	cout << "REPORT " << tree_count << " cache_size " << db.gspan.getCache().size() << endl;
	cout << "REPORT " << tree_count << " gain_count " << gain_count << endl;
	cout << "REPORT " << tree_count << " bound_count " << bound_count << endl;
}

void CLASS::reportFeatureImportance() {
	size_t count = 0;
	for (const auto& pair : db.gspan.getCache()) {
		const auto& pattern = pair.first;
		double importance = pair.second.feature_importance;
		if (importance > 0) {
			count++;
			cout << "REPORT " << tree_count << " FI " << pattern << " : " << importance << endl;
		}
	}
	cout << "REPORT " << tree_count << " used_features " << count << endl;
}

