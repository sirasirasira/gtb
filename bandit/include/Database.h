#pragma once

#include "Structures.h"
#include "Setting.h"

#include "GradientBoosting.h"
#include "Planter.h"
#include "Spliter.h"
#include "Finder.h"
#include "Gspan.h"
#include "Evaluater.h"

struct Database {
	GraphData gdata;
	vector<double> raw_ys; // true value
	vector<double> ys; // residual error
	vector<double> y_predictions; // predict value

	GradientBoosting gradient_boosting;
	Planter planter;
	Spliter spliter;
	Finder finder;
	Gspan gspan;
	Evaluater eva;
};
