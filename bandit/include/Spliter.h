#pragma once

#include "MyInclude.h"
#include "StructuresGspan.h"

class Spliter {
	public:
		void prepare(const vector<ID>& _targets);
		vector<ID> run(const vector<ID>& _targets);
		void update(Pattern& pattern, double score);
		bool isBounded(double min_bound);
		inline bool valid() {
			return valid_flg;
		}
		inline const Pattern& getBestPattern() {
			assert(valid_flg);
			return best_pattern;
		}
		inline double getImportance() {
			assert(valid_flg);
			return parent_score - min_score;
		}

	private:
		bool valid_flg;
		vector<ID> targets;
		double parent_score;
		double min_score;
		Pattern best_pattern;

		void initMinScore();
};

