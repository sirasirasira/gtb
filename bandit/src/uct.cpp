#include "UCT.h"
#define CLASS UCT

extern Setting setting;
#include "Database.h"
#include "Calculator.h" 
extern Database db;

void CLASS::run(const vector<ID>& _targets) {
	// std::cout << "uct run" << std::endl; // debug
	targets = _targets;
	db.gspan.clearUCB();
	Pattern pattern;
	vector<ID> posi;
	double score;

	for (unsigned int i = 0; i < setting.iteration; i++) {
		path = {root};
		if(!selection(root)) { // all node is searched
			break;
		}
		expansion();
		pattern = path[path.size()-1];
		pattern = simulation(pattern);
		if (cache.find(pattern) == cache.end()) {
			posi = db.finder.run(pattern, targets);
		} else {
			posi = db.gspan.getPosiIds(cache[pattern].g2tracers);
		}
		score = Calculator::score(db.ys, targets, posi);
		backpropagation(score);
	}
}

bool CLASS::selection(const Pattern& pattern) {
	// std::cout << "selection: " << pattern << std::endl; // debug
	if (cache[pattern].terminal) return true;

	Pattern best_child;
	double max_ucb = -DBL_MAX;
	double ucb;
	for (auto& c : cache[pattern].childs) {
		if (cache[c].prune) {
			continue;
		} else {
			if (cache[c].count == 0) {
				if (update(c)) { // prune
					continue;
				} else { // not prune
					best_child = c;
					break;
				}
			} else {
				ucb = (cache[c].sum_score / cache[c].count)
					+ setting.exploration_strength
					* (sqrt(2 * log(cache[pattern].count) / cache[c].count));
				ucb -= setting.bound_rate * cache[c].bound; //TODO
			}

			if (ucb > max_ucb) {
				best_child = c;
				max_ucb = ucb;
			}
		}
	}

	if (best_child.size() != 0) {
		path.push_back(best_child);
		return selection(best_child);
	} else {
		if (path.size() == 1) { // all node is pruned
			return false;
		} else {
			cache[pattern].prune = true;
			path.pop_back();
			return selection(path[path.size()-1]);
		}
	}
}

bool CLASS::update(const Pattern& pattern) {
	vector<ID> posi;
	double score;
	double bound;
	posi = db.gspan.getPosiIds(cache[pattern].g2tracers);
	score = Calculator::score(db.ys, targets, posi);
	db.spliter.update(pattern, score);
	bound = Calculator::bound(db.ys, targets, posi);
	cache[pattern].bound = bound;
	if (db.spliter.isBounded(bound)){
		cache[pattern].prune = true;
		return true;
	}
	return false;
}

void CLASS::expansion() {
	const Pattern pattern = path[path.size()-1];
	// std::cout << "expansion: " << pattern << endl; // debug
	if (cache[pattern].count >= setting.threshold) {
		if (!cache[pattern].scan) {
			if (db.gspan.scanGspan(pattern)) {
				cache[pattern].terminal = false;
				expand_selection(pattern);
			}
		} else {
			if (cache[pattern].childs.size()) {
				cache[pattern].terminal = false;
				expand_selection(pattern);
			}
		}
	}
}

void CLASS::expand_selection(const Pattern& pattern) {
	Pattern selected_child;
	for (auto& c : cache[pattern].childs) {
		if (update(c)) { // prune
			continue;
		} else { // not prune
			selected_child = c;
			break;
		}
	}

	if (selected_child.size() != 0) {
		path.push_back(selected_child);
	} else {
		cache[pattern].prune = true;
	}
}

Pattern CLASS::simulation(const Pattern& pattern) {
	// cout << "simulation: " << pattern << endl;
	auto& g2tracers = cache[pattern].g2tracers;
	// gid is in g2tracers and in targets
	vector<ID> keys(g2tracers.size());
	size_t i = 0;
	for (auto& x : g2tracers) {
		keys[i] = x.first;
		i++;
	}
	auto gids = Calculator::setIntersec(targets, keys);
	ID gid = gids[Dice::id(gids.size())];
	auto& tracers = g2tracers[gid];
	auto& tracer = tracers[Dice::id(tracers.size())];
	
	return db.gspan.EdgeSimulation(pattern, tracer, gid);
}
 
void CLASS::backpropagation(double score) {
	// cout << "backpropagation: " << score << endl;
	for (int i = path.size() - 1; i > 0; i--) {
		cache[path[i]].count++;
		cache[path[i]].sum_score -= score; // score is mse (lower is better)
	}
	cache[root].count++;
}
