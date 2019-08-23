#pragma once

#include "MyInclude.h"
#include "Gspan.h"
#include "StructuresGspan.h"

class UCT {
	public:
		void run(const vector<ID>& _targets);

	private:
		vector<ID> targets;
		map<Pattern, CacheRecord> cache;
		Pattern pattern;
		Pattern root;
		vector<Pattern> path;

		Pattern selection();
		void expansion();
		Pattern simulation(const Pattern&);
		void backpropagation();
		bool stop_condition(const Pattern&);
};
