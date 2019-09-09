#pragma once

#include "IncludeLib.h"
#include "Structures.h"
#include "StructuresGspan.h"

struct PQRecord {
	Pattern pattern;
	vector<ID> posi;

	PQRecord() {}
	PQRecord(Pattern pattern, vector<ID> posi)
		: pattern(pattern), posi(posi) {}
};

inline bool operator < (const PQRecord& l, const PQRecord& r) {
	return (l.pattern < r.pattern);
}
