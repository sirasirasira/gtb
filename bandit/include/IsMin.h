#pragma once

#include "MyInclude.h"
#include "StructuresGspan.h"

class IsMin {
	public:
		Pattern convert(Pattern& pattern);
	private:
		const Pattern* pattern_ptr;
		Pattern minChecker(Pattern& comp, Graph& g, Tracers& tracers);
};
