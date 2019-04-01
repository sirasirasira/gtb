#pragma once

#include "MyInclude.h"

class Planter {
	public:
		const vector<double>& run(const vector<ID>& targets);

		void resizeAdditiveYs(size_t s) {
			additive_ys.resize(s);
		}

	private:
		vector<double> additive_ys;
		void grow(vector<ID> targets, size_t depth);
		inline void updataFeatureImportance();
		void makeLeaf(const vector<ID>& targets, size_t depth);
		inline void updateAdditiveYs(const vector<ID>& targets, double mean);
		inline bool checkLeafConditions(const vector<ID>& targets, size_t depth);
};
