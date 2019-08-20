#pragma once

#include "MyInclude.h"
#include "Gspan.h"
#include "StructuresGspan.h"

class UCT {
	public:
		void run(const vector<ID>& _targets);

	private:
		vector<ID> targets;
		map<Pattern, CacheRecord>* cache;
		vector<Pattern>* e1patterns;
		Pattern pattern;
		vector<Pattern> path;

		void selection();
		void expansion();
		void simulation();
		void backpropagation();
};
