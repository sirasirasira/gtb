#pragma once

#include "MyInclude.h"
#include "StructuresGspan.h"

class IsMin {
	public:
		bool run(const Pattern& pattern);
	private:
		const Pattern* pattern_ptr;
		bool minChecker(Pattern& comp, Graph& g, Tracers& tracers);
};
