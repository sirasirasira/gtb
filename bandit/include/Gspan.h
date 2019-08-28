#pragma once

#include "MyInclude.h"
#include "StructuresGspan.h"
#include "Spliter.h"
#include "IsMin.h"

struct CacheRecord {
	GraphToTracers g2tracers;
	vector<Pattern> childs;
	bool terminal;
	double bound;
	size_t count;
	double sum_score;
	double ucb;
	double feature_importance;
	CacheRecord() {
		terminal = true;
		bound = 0;
		count = 0;
		sum_score = 0;
		ucb = DBL_MAX;
		feature_importance = 0;
	}
	CacheRecord(vector<Pattern> childs)
		: childs(childs){
			terminal = true;
			bound = 0;
			count = 0;
			sum_score = 0;
			ucb = 0;
			feature_importance = 0;
		}
	CacheRecord(GraphToTracers g2tracers, vector<Pattern> childs)
		: g2tracers(g2tracers), childs(childs){
			terminal = true;
			bound = 0;
			count = 0;
			sum_score = 0;
			ucb = 0;
			feature_importance = 0;
		}
};

class Gspan {

	public:
		size_t minsup;
		size_t maxpat;

		void run();
		void run(Pattern pattern);

		void makeRoot();

		inline const map<Pattern, CacheRecord>& getCache() {
			return cache;
		}

		inline const Pattern& getRoot() {
			return root;
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
			for (auto& x : cache) {
				x.second.terminal = true;
				x.second.bound = 0;
				x.second.count = 0;
				x.second.sum_score = 0;
				x.second.ucb = 0;
			}
		}

		tuple<Pattern, EdgeTracer, ID> oneEdgeSimulation(tuple<Pattern, EdgeTracer, ID>&);
		bool scanGspan(const Pattern&);

	private:
		Spliter* spliter;
		Pattern pattern;
		Pattern root;
		IsMin is_min;
		map<Pattern, CacheRecord> cache; // inserted data must keep pointer

		void edgeGrow(GraphToTracers& g2tracers, bool in_cache_flg = false);
		size_t support(GraphToTracers& g2tracers);
		void report(GraphToTracers& tracers);

};
