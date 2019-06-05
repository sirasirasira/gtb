#pragma once

#include "MyInclude.h"

class Planter {
	public:
		const vector<double>& run(const vector<ID>& train_targets, const vector<ID>& test_targets);

		void resizeAdditiveYs(size_t s) {
			additive_ys.resize(s);
		}

	private:
		vector<double> additive_ys;
		void grow(const vector<ID>& train_targets, const vector<ID>& test_targets, size_t depth);
		inline bool checkLeafConditions(const vector<ID>& targets, size_t depth);
		void makeLeaf(const vector<ID>& train_targets, const vector<ID>& test_targets, size_t depth);
		inline void updateAdditiveYs(const vector<ID>& train_targets, const vector<ID>& test_targets, double mean);
};
