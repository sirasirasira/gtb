#pragma once

#include "MyInclude.h"
#include "StructuresGspan.h"

class Finder {
	public:
		vector<ID> run(const Pattern& best_pattern, const vector<ID>& targets);

	private:
		Pattern pattern;
		Graph g;
		bool exist;
		vector<ID> vids;

		bool gcheck(const ID gid);
		void gtraverse(unsigned int cur);
		bool ndiscovered(const int id);
};
