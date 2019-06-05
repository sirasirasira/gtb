#include "Spliter.h"
#define CLASS Spliter

extern Setting setting;
#include "Database.h"
extern Database db;

#include "Calculator.h" 
#include "StructuresGspan.h"

void CLASS::initMinScore() {
	parent_score = Calculator::imp(db.ys, Calculator::trainOnly(targets));
	min_score = parent_score - setting.needed_impurity_decrease - std::numeric_limits<double>::epsilon();
}

void CLASS::prepare(const vector<ID>& _targets) {
	// std::cout << "debug Spliter prepare" << std::endl; // debug
	targets = _targets;
	initMinScore();
	// Debug::IDs(targets); // debug
	db.gspan.setSpliterPtr(this);
	db.gspan.minsup = setting.minsup;
	db.gspan.maxpat = 3;
	db.gspan.run();
	db.gspan.maxpat = setting.maxpat;
	std::cout << "prepare cache size: " << db.gspan.getCache().size() << std::endl;
}

vector<ID> CLASS::run(const vector<ID>& _targets) {
	// std::cout << "debug spliter run" << std::endl; // debug
	targets = _targets;
	best_pattern = {};
	initMinScore();
	if (min_score < 0) {
		goto G_INVALID;
	}
	pq_enum = {};
	searchCache();
	searchEnum();
	// std::cout << "debug best_pattern " << best_pattern << std::endl; // debug
	if (best_pattern.size() == 0) {
G_INVALID:
		valid_flg = false;
		return {};
	}
	// std::cout << "debug parent_score " << parent_score << " min_score " << min_score << std::endl; // debug
	valid_flg = true;
	vector<ID> posi = db.gspan.getPosiIds(db.gspan.getCache().at(best_pattern).g2tracers);
	return Calculator::setIntersec(targets, posi);
}

void CLASS::searchCache() {
	// std::cout << "serch Cache" << std::endl;
	const auto& cache = db.gspan.getCache();
	for (auto itr = cache.begin(); itr != cache.end(); itr++) {
		const auto& pattern = itr->first;
		const auto& g2tracers = itr->second.g2tracers;
		vector<ID> posi = db.gspan.getPosiIds(g2tracers);
		if (itr->second.childs.size() == 0) {
			double min_bound = Calculator::bound(db.ys, targets, posi);
			pq_enum.push(std::make_pair(min_bound, PQRecord(pattern, posi)));
		} else {
			update(pattern, posi);
		}
	}
}

void CLASS::searchEnum() {
	// std::cout << "search Enum" << std::endl;
	while (!pq_enum.empty()) {
		double min_bound = pq_enum.top().first;
		if (min_score <= min_bound) {
			break;
		}

		PQRecord pqrecord = pq_enum.top().second;
		update(pqrecord.pattern, pqrecord.posi);
		pq_enum.pop();
		db.gspan.run(pqrecord.pattern);
	}
}

void CLASS::update(Pattern pattern, vector<ID> posi) {
	if (Dice::p(setting.feature_used) == false) {
		return;
	}
	double score = Calculator::score(db.ys, targets, posi);
	if (score < min_score ) { // old pattern may be used (this func is called from gspan)
		min_score = score;
		best_pattern = pattern;
	}
}

void CLASS::push_pq_enum(double bound, PQRecord pqrecord) {
	pq_enum.push(std::make_pair(bound, pqrecord));
}

bool CLASS::isBounded(vector<ID> posi) {
	double min_bound = Calculator::bound(db.ys, targets, posi);
	if (min_score <= min_bound) {
		return true;
	} else {
		return false;
	}
}

