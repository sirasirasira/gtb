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
		Pattern root;
		vector<Pattern> path;

		void selection(const Pattern&);
		void expansion();
		bool update();
		Pattern simulation(const Pattern&);
		void backpropagation(double);
};
