#pragma once

#include "MyInclude.h"
#include "StructuresGspan.h"
#include "Spliter.h"
#include "IsMin.h"

class Gspan {

	struct CacheRecord {
		GraphToTracers g2tracers;
		vector<std::pair<DFSCode, vector<ID>>> childs;
		double feature_importance;
		CacheRecord() {
			feature_importance = 0;
		}
		CacheRecord(GraphToTracers g2tracers, vector<std::pair<DFSCode, vector<ID>>> childs)
			: g2tracers(g2tracers), childs(childs) {
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

		inline void reportCache() {
			std::cout << "Report Cache Infomation" << std::endl;
			for (auto itr = cache.begin(); itr != cache.end(); ++itr) {
				std::cout << "pattern : " << itr->first << std::endl;
			}
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

	private:
		Spliter* spliter;
		Pattern pattern;
		IsMin is_min;
		map<Pattern, CacheRecord> cache; // inserted data must keep pointer

		void edgeGrow();
		void oneedgeGrow();
		size_t support(GraphToTracers& g2tracers);
		bool check_pattern(Pattern pattern, GraphToTracers& g2tracers);
		void report(GraphToTracers& tracers);
		int scanGspan(GraphToTracers& g2tracers, PairSorter& b_heap, map<int, PairSorter, std::greater<int>>& f_heap) const ;

};
