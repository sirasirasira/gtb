#pragma once

#include "MyInclude.h"
#include "StructuresGspan.h"

class IsMin {
	public:
		pair<Pattern, EdgeTracer> convert(const EdgeTracer& tracer, ID gid);
	private:
		// const Pattern* pattern_ptr;
		pair<Pattern, EdgeTracer> minChecker(Pattern& comp, Graph& g, Tracers& tracers);
};
