#pragma once

#include "IncludeLib.h"

using ID = size_t; // starts from 0

struct Triplet {
	int x;
	int y;
	int z;
	inline Triplet reverse() {
		return Triplet(z, y, x);
	}
	explicit Triplet() {};
	explicit Triplet(int x, int y, int z) : x(x), y(y), z(z) {};
};

inline std::ostream& operator << (std::ostream& os, const Triplet t) {
	os << t.x << "," << t.y << "," << t.z;
	return os;
}

inline bool operator < (const Triplet& l, const Triplet& r) {
	if (l.x != -1 and r.x != -1 and l.x != r.x) return (l.x < r.x);
	if (l.y != -1 and r.y != -1 and l.y != r.y) return (l.y < r.y);
	return (l.z < r.z);
}

inline bool operator <= (const Triplet& l, const Triplet& r) {
	return !(r < l);
}

inline bool operator == (const Triplet& l, const Triplet& r) {
	return (l.x == r.x and l.y == r.y and l.z == r.z);
}

struct Pair {
	ID a;
	ID b;
	void set(ID _a, ID _b) {
		a = _a;
		b = _b;
	}
};

inline bool operator < (const Pair& l, const Pair& r) {
	if (l.a != r.a) return (l.a < r.a);
	return (l.b < r.b);
}

struct VertexPair : public Pair {
	ID id;
};

struct Edge {
	ID id; // edge id
	ID to; // to vertex id (from vertex id is AdjacentList key)
	Triplet labels;

	Edge() {};
	Edge(ID id, ID to, Triplet labels) : id(id), to(to), labels(labels) {};
};

using AdjacentList = vector<vector<Edge>>;

struct Graph : public AdjacentList {
	vector<int> label;
	size_t num_of_edges;

	inline void resize(size_t s) {
		AdjacentList::resize(s);
		label.resize(s);
	}
};

struct GraphData : vector<Graph> {
	size_t num_train;
	size_t num_test;

	inline ID getFirstTrainID() {
		return 0;
	}
	inline ID getLastTrainID() {
		return num_train - 1;
	}
	inline ID getFirstTestID() {
		return num_train;
	}
	inline ID getLastTestID() {
		return size() - 1;
	}
};
