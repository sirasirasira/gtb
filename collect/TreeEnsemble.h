#pragma once

#include "MyInclude.h"

class TreeEnsemble {
	public:
		void runRandomForest();
		void runGradientBoosting();
		void incGainCount() {
			gain_count++;
		}
		void incBoundCount() {
			bound_count++;
		}
		void printGainCount() {
			cout << "gain_count " << gain_count << endl;
		}

	private:
		vector<ID> targets;
		size_t tree_count; // (gradient boosting) from 0, (random forest) from 1
		int gain_count = 0;
		int bound_count = 0;

		void makeTargets();
		void plantFirst();
		void calcResidualErrors();
		void plant();
		inline void report();
		inline void reportFeatureImportance();

		const vector<ID>& getTargets() {
			return targets;
		}

};
