#pragma once

#include "MyInclude.h"

class GradientBoosting {
	public:
		void run();
		void incGainCount() {
			gain_count++;
		}
		void incBoundCount() {
			bound_count++;
		}

	private:
		vector<ID> train_targets;
		vector<ID> test_targets;
		size_t tree_count; // (gradient boosting) from 0, (random forest) from 1
		int gain_count = 0;
		int bound_count = 0;

		void makeTargets();
		void plantFirst();
		void calcResidualErrors();
		inline void report();
		inline void reportFeatureImportance();
};
