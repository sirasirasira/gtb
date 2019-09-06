#pragma once

#include "MyInclude.h"
#include "StructuresGspan.h"
#include "Spliter.h"
#include "IsMin.h"

class Gspan {

	struct CacheRecord {
		GraphToTracers g2tracers;
		double feature_importance;
		CacheRecord() {
			feature_importance = 0;
		}
		CacheRecord(GraphToTracers g2tracers) : g2tracers(g2tracers) {
			feature_importance = 0;
		}
	};

	public:
		size_t minsup;
		size_t maxpat;

		void run();
		void run(Pattern pattern);

		inline const map<Pattern, CacheRecord>& getCache() {
			return cache;
		}

		inline void updateFeatureImportance(const Pattern& pattern, double importance) {
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

	private:
		Spliter* spliter;
		Pattern pattern;
		IsMin is_min;
		map<Pattern, CacheRecord> cache; // inserted data must keep pointer

		void edgeGrow(GraphToTracers& g2tracers, bool in_cache_flg = false);
		size_t support(GraphToTracers& g2tracers);
		void report(GraphToTracers& tracers);
		int scanGspan(GraphToTracers& g2tracers, PairSorter& b_heap, map<int, PairSorter, std::greater<int>>& f_heap) const ;

};
