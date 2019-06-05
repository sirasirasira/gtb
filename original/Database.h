#pragma once

#include "Structures.h"
#include "Setting.h"

#include "TreeEnsemble.h"
#include "Planter.h"
#include "Spliter.h"
#include "Gspan.h"
#include "Evaluater.h"

struct Database {
	GraphData gdata;
	vector<double> raw_ys;
	vector<double> ys; // (gradient boosting) residual error, (random forest) raw ys 
	vector<double> y_predictions;

	TreeEnsemble tree_ensemble;
	Planter planter;
	Spliter spliter;
	Gspan gspan;
	Evaluater eva;
};
