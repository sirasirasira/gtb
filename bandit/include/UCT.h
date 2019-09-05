#pragma once

#include "MyInclude.h"
#include "Gspan.h"
#include "StructuresGspan.h"

class UCT {
	public:
		UCT(map<Pattern, CacheRecord>& _cache, const Pattern& _root) : cache(_cache), root(_root) {
		}
		void run(const vector<ID>& _targets);

	private:
		vector<ID> targets;
		map<Pattern, CacheRecord>& cache;
		const Pattern& root;
		vector<Pattern> path;

		void selection(const Pattern&);
		void expansion();
		bool update();
		Pattern simulation(const Pattern&);
		void backpropagation(double);
};
