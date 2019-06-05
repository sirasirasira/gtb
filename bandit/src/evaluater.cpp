#include "Evaluater.h"
#define CLASS Evaluater

extern Setting setting;
#include "Database.h"
extern Database db;

#include "Calculator.h" 

void CLASS::run(size_t tree_count, vector<ID>& train_targets, vector<ID>& test_targets) {
	calcACCAUCLoss(tree_count, "train", train_targets);
	calcACCAUCLoss(tree_count, "test", test_targets);
}

// @change pred_map
void CLASS::calcACCAUCLoss(size_t tree_count, string type, vector<ID>& targets) {
	size_t num_correct = 0;
	size_t num_all = targets.size();
	double loss_sum = 0;
	map<double, vector<double>> pred_map;
	for (auto gid : targets) {
		double y = db.raw_ys[gid];
		double p = db.y_predictions[gid];
		pred_map[p].push_back(y);
		if (Calculator::isSameSign(y, p)) num_correct++;
		loss_sum += Calculator::calcDeviation(y, p);
	}
	double acc = num_correct / (double) num_all;
	double auc = calcAUC(pred_map);
	double loss_mean = loss_sum / (double) num_all;
	cout << "REPORT " << tree_count << " " << type << "acc " << acc << endl;
	cout << "REPORT " << tree_count << " " << type << "auc " << auc << endl;
	cout << "REPORT " << tree_count << " " << type << "loss_mean " << loss_mean << endl;
}

double CLASS::calcAUC(const map<double, vector<double>>& pred_map) {
	double area = 0;
	int height = 0;
	int width = 0;
	for (const auto& pair : pred_map) {
		const auto& vec = pair.second;
		int count_true_label = 0;
		int count_false_label = 0;
		for (double label : vec) {
			if (label > 0) {
				count_true_label++;
			} else {
				count_false_label++;
			}
		}
		area += height * count_false_label + count_true_label * count_false_label / 2.0;
		height += count_true_label;
		width += count_false_label;
	}
	return 1 - (area / (double) height / (double) width);
}

