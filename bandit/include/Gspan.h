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
	bool scan;
	CacheRecord() {
		terminal = true;
		bound = 0;
		count = 0;
		sum_score = 0;
		ucb = DBL_MAX;
		feature_importance = 0;
		scan = false;
	}
	CacheRecord(vector<Pattern> childs)
		: childs(childs){
			terminal = true;
			bound = 0;
			count = 0;
			sum_score = 0;
			ucb = DBL_MAX;
			feature_importance = 0;
			scan = false;
		}
	CacheRecord(GraphToTracers g2tracers, vector<Pattern> childs)
		: g2tracers(g2tracers), childs(childs){
			terminal = true;
			bound = 0;
			count = 0;
			sum_score = 0;
			ucb = DBL_MAX;
			feature_importance = 0;
			scan = false;
		}
};

class Gspan {

	public:
		size_t minsup;
		size_t maxpat;

		Pattern EdgeSimulation(const Pattern&, EdgeTracer&, ID);
		bool scanGspan(const Pattern&);
		void makeRoot(const vector<ID>&);
		void one_edge_report(GraphToTracers& g2tracers);

		inline map<Pattern, CacheRecord>& getCache() {
			return cache;
		}

		inline const Pattern& getRoot() {
			return root;
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

		inline void clearUCB() {
			for (auto& x : cache) {
				x.second.terminal = true;
				x.second.bound = 0;
				x.second.count = 0;
				x.second.sum_score = 0;
				x.second.ucb = DBL_MAX;
			}
			cache[root].terminal = false;
		}

	private:
		Spliter* spliter;
		Pattern root;
		IsMin is_min;
		map<Pattern, CacheRecord> cache; // inserted data must keep pointer

		bool stop_condition(const Pattern, bool);
		size_t support(GraphToTracers& g2tracers);

};
