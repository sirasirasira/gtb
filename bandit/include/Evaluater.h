#pragma once

#include "MyInclude.h"

class Evaluater {
	public:
		void run(size_t t, vector<ID>& train_targets, vector<ID>& test_targets);

	private:
		void calcACCAUCLoss(size_t tree_count, string type, vector<ID>& targets);
		double calcAUC(const map<double, vector<double>>& pred_map);
};
