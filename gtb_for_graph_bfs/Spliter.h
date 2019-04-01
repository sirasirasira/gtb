#pragma once

#include "MyInclude.h"
#include "StructuresSpliter.h"

class Spliter {

	public:
		void prepare(const vector<ID>& _targets);
		vector<ID> run(const vector<ID>& _targets);
		void update(Pattern pattern, vector<ID> posi);
		void push_pq_enum(double bound, PQRecord pqrecord);
		bool isBounded(vector<ID> posi);
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
		inline vector<ID> gettargets() {
			return  targets;
		}

	private:
		bool valid_flg;
		vector<ID> targets;
		double parent_score;
		double min_score;
		Pattern best_pattern;
		priority_queue<std::pair<double, PQRecord>, vector<std::pair<double, PQRecord>>, std::greater<std::pair<double, PQRecord>>> pq_cache;
		priority_queue<std::pair<double, PQRecord>, vector<std::pair<double, PQRecord>>, std::greater<std::pair<double, PQRecord>>> pq_enum;

		void initMinScore();
		void searchCache();
		inline bool isFrontier(const Pattern& pattern, const Pattern& bef_pattern);
		void searchEnum();
};

