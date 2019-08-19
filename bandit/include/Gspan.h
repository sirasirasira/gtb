#pragma once

#include "MyInclude.h"
#include "StructuresGspan.h"
#include "Spliter.h"
#include "IsMin.h"

class Gspan {

	struct CacheRecord {
		vector<ID> posi;
		vector<DFSCode> childs;
		double bound;
		int count;
		int sum_value;
		double ucb;
		double feature_importance;
		CacheRecord() {
			bound = 0;
			count = 0;
			sum_value = 0;
			ucb = 0;
			feature_importance = 0;
		}
		CacheRecord(vector<ID> posi, vector<DFSCode> childs, double bound, int count, int sum_value, double ucb)
			: posi(posi), childs(childs) , bound(bound), count(count), sum_value(sum_value), ucb(ucb){
			feature_importance = 0;
		}
	};

	public:
		size_t minsup;
		size_t maxpat;

		void run();
		void run(Pattern pattern);

		void searche1Patterns();

		inline const map<Pattern, CacheRecord>& getCache() {
			return cache;
		}

		inline const vector<Pattern>& gete1Patterns() {
			return e1patterns;
		}

		inline void updataFeatureImportance(const Pattern& pattern, double importance) {
			cache.at(pattern).feature_importance += importance;
		}

		inline vector<ID> getPosiIds(const GraphToTracers& tracers) {
			vector<ID> vec(tracers.size());
			size_t i = 0;
			for (auto x : tracers) {
				vec[i] = x.first;
				i++;
			}
			return vec;
		}

		void setSpliterPtr(Spliter* ptr) {
			spliter = ptr;
		}

		inline void setMinsup(size_t _minsup) {
			minsup = _minsup;
		}

		inline void setMaxpat(size_t _maxpat) {
			maxpat = _maxpat;
		}

		inline void clearUCB() {
			for (auto [key, value] : cache) {
				value.bound = 0
				value.count = 0;
				value.sum_value = 0;
				value.ucb = 0;
			}
		}

	private:
		Spliter* spliter;
		Pattern pattern;
		IsMin is_min;
		map<Pattern, CacheRecord> cache; // inserted data must keep pointer
		vector<Pattern> e1patterns;

		void edgeGrow(GraphToTracers& g2tracers, bool in_cache_flg = false);
		size_t support(GraphToTracers& g2tracers);
		void report(GraphToTracers& tracers);
		int scanGspan(GraphToTracers& g2tracers, PairSorter& b_heap, map<int, PairSorter, std::greater<int>>& f_heap) const ;

};
