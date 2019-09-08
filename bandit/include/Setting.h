#pragma once

#include "IncludeLib.h"

struct Setting {
	//size_t minsup = 1; // same as min_node_size
	double minsup = 1; // 1~: support size, 0~1: support rate
	size_t maxpat = -1;
	size_t num_of_trees = 300;
	double shrinkage = 1.0;
	double needed_impurity_decrease = 0.0;
	size_t max_depth = -1;
	size_t iteration = 100;
	size_t threshold = 1;
	double exploration_strength = 1;
	double stopping_rate = 1 - pow(10 , -1);

	void print() {
		cout
			<< " minsup (min_node_size) " << minsup
			<< " maxpat " << maxpat
			<< " num_of_trees " << num_of_trees
			<< " shrinkage " << shrinkage
			<< " needed_impurity_decrease " << needed_impurity_decrease
			<< " max_depth " << max_depth
			<< " iteration " << iteration
			<< " threshold " << threshold
			<< " exploration_strength " << exploration_strength
			<< " stopping_rate " << stopping_rate
			<< endl;
	}

};
