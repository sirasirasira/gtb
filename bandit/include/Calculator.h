#pragma once

#include "MyInclude.h"
#include "StructuresGspan.h"
 
extern Setting setting;
#include "Database.h"
extern Database db;

namespace Calculator {

	inline bool isSameSign(double a, double b) {
		if (a * b > 0) {
			return true;
		} else {
			return false;
		}
	}

	inline double calcDeviation(double ans, double pred) {
		using namespace std;
		return log(1 + exp(-2 * ans * pred));
	}

	inline double calcResidualErr(double ans, double pred) {
		using namespace std;
		return -1 * (-2 * ans) / (exp(2 * ans * pred) + 1);
	}

	inline vector<ID> setDiff(const vector<ID>& a, const vector<ID>& b) {
		vector<ID> vec;
		std::set_difference(a.begin(), a.end(), b.begin(), b.end(), std::back_inserter(vec));
		return vec;
	}

	inline vector<ID> setIntersec(vector<ID> a, vector<ID> b) {
		vector<ID> vec;
		std::set_intersection(a.begin(), a.end(), b.begin(), b.end(), std::back_inserter(vec));
		return vec;
	}

	// not sorted targets is ok
	inline double calcYsMean(const vector<double>& ys, const vector<ID>& train_targets) {
		double sum = 0;
		for (ID id : train_targets) {
			sum += ys[id];
		}
		double mean = sum / (double) train_targets.size();
		return mean;
	}

	// not sorted targets is ok
	inline double TSD(const vector<double>& ys, const vector<ID>& train_targets) {
		double mean = calcYsMean(ys, train_targets);
		double tss = 0;
		for (ID id : train_targets) {
			tss += calcDeviation(ys[id], mean);
		}
		return tss;
	}

	// not sorted targets is ok
	inline double TSS(const vector<double>& ys, const vector<ID>& train_targets) {
		double mean = calcYsMean(ys, train_targets);
		double tss = 0;
		for (ID id : train_targets) {
			tss += pow(ys[id] - mean, 2);
		}
		return tss;
	}

	inline double imp(const vector<double>& ys, const vector<ID>& train_targets) {
		return TSS(ys, train_targets);
	}

	// assert raw is sorted
	inline vector<ID> trainOnly(const vector<ID>& raw) {
		vector<ID> only;
		for (ID id : raw) {
			if ((id <= db.gdata.getLastTrainID()) == false) {
				break;
			}
			only.push_back(id);
		}
		return only;
	}

	// assert raw_* is sorted
	inline double score(const vector<double>& ys, const vector<ID>& targets, const vector<ID>& raw_posi) {
		db.gradient_boosting.incGainCount();
		vector<ID> posi = setIntersec(targets, raw_posi); // TODO
		/*
		if (posi.size() < db.gspan.minsup) {
			return DBL_MAX;
		}
		*/
		vector<ID> nega = setDiff(targets, posi);
		/*
		if (nega.size() < db.gspan.minsup) {
			return DBL_MAX;
		}
		*/
		return imp(ys, posi) + imp(ys, nega);
	}

	inline vector<ID> myCast(multimap<double, ID> sorted_posi_ids) {
		vector<ID> posi_ids;
		for (auto p : sorted_posi_ids) {
			posi_ids.push_back(p.second);
		}
		return posi_ids;
	}

	// assert raw_* is sorted
	inline double bound(const vector<double>& ys, const vector<ID>& targets, const vector<ID>& raw_posi) {
		db.gradient_boosting.incBoundCount();
		double min_score = DBL_MAX;
		vector<ID> posi = setIntersec(targets, raw_posi); // TODO
		multimap<double, ID> sorted_posi_ids;
		vector<ID> posi_ids;
		vector<ID> nega_ids;
		for (auto id : posi) {
			sorted_posi_ids.insert({ys[id], id});
		}
		for (int a = 0; a < 2; a++) {
			// move from back
			nega_ids = setDiff(targets, posi); // not depend on posi_ids
			posi_ids = myCast(sorted_posi_ids);
			if (a == 0) {
				// move from front
				std::reverse(begin(posi_ids), end(posi_ids));
			}
			for (size_t i = posi_ids.size(); i >= db.gspan.minsup; i--) {
				ID move_id = *(posi_ids.rbegin());
				posi_ids.pop_back();
				nega_ids.push_back(move_id);
				double score = imp(ys, posi_ids) + imp(ys, nega_ids);
				if (score < min_score) min_score = score;
			}
		}
		return min_score;
	}

}
