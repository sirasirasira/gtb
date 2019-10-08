#include "MyInclude.h"
#include "Database.h"
#include "getopt.h"

Setting setting;
Database db;

void read(std::istream& is) {
	GraphData& graphs = db.gdata;
	Graph g;
	double y;
	Triplet labels;
	Edge edge;
	char c;
	string line;
	size_t line_num = 0;
	ID vid;
	ID eid = 0;
	while (getline(is, line)) {
		line_num++;
		if (line.empty()) {
			g.num_of_edges = eid;
			graphs.push_back(g);
		}
		std::stringstream stream(line);
		if (line[0] == 't') {
			g.clear();
			eid = 0;
			int i;
			stream >> c >> c >> i >> y;
			db.raw_ys.push_back(y);
		} else if (line[0] == 'v') {
			int l;
			stream >> c >> vid >> l;
			if (vid > (ID) g.size()) {
				std::cerr << "vertex ids are not sorted! line: " << line_num << " vid " << vid << " size " << g.size() << endl;
				std::exit(-1);
			}
			g.resize(vid + 1);
			g.label[vid] = l;
		} else if (line[0] == 'e') {
			size_t v1;
			size_t v2;
			stream >> c >> v1 >> v2 >> labels.y;
			labels.x = g.label[v1];
			labels.z = g.label[v2];
			g[v1].push_back(Edge(eid++, v2, labels));
			g[v2].push_back(Edge(eid++, v1, labels.reverse()));
		}
	}
}

void readData(std::istream& train_is, std::istream& test_is) {
	read(train_is);
	db.gdata.num_train = db.gdata.size();
	read(test_is);
	db.gdata.num_test = db.gdata.size() - db.gdata.num_train;
	db.y_predictions.resize(db.gdata.size(), 0);
	//Debug::ys(db.y_predictions); // debug
	db.planter.resizeAdditiveYs(db.gdata.size());
}

inline vector<ID> allTargets() {
	vector<ID> all_targets(db.gdata.size());
	ID id = 0;
	for (auto& val : all_targets) {
		val = id;
		id++;
	}
	return all_targets;
}

int main(int argc, char** argv) {
	int opt;
	string usage = "Usage: " + std::string(argv[0])
		+ " [-m minsup]"
		+ " [-x maxpat]"
		+ " [-t num_of_trees]"
		+ " [-s (double) shrinkage]"
		+ " [-n (double) needed_impurity_decrease]"
		+ " [-d max_depth]"
		+ " [-i iteration]"
		+ " [-h threshold]"
		+ " [-c exploration_strength]"
		+ " [-p stopping_rate]"
		+ " [-b UCB:bound_rate]"
		+ " traing_data_file"
		+ " test_data_file"
		+ "";
	while ((opt = getopt(argc, argv, "m:x:t:s:n:d:i:h:c:p:b:re")) != -1) {
		switch (opt) {
			case 'm':
				//setting.minsup = atoi(optarg);
				setting.minsup = atof(optarg);
				break;
			case 'x':
				setting.maxpat = atoi(optarg);
				break;
			case 't':
				setting.num_of_trees = atoi(optarg);
				break;
			case 's':
				setting.shrinkage = atof(optarg);
				break;
			case 'n':
				setting.needed_impurity_decrease = atof(optarg);
				break;
			case 'd':
				setting.max_depth = atoi(optarg);
				break;
			case 'i':
				setting.iteration = atoi(optarg);
				break;
			case 'h':
				setting.threshold = atoi(optarg);
				break;
			case 'c':
				setting.exploration_strength = atof(optarg);
				break;
			case 'p':
				setting.stopping_rate = 1 - pow(10, atoi(optarg));
				break;
			case 'b':
				setting.bound_rate = atof(optarg);
				break;
			default:
				std::cerr << usage << std::endl;
				return -1;
		}
	}

	if (argc - optind != 2) {
		std::cerr << usage << std::endl;
		return -1;
	}
	string train_filename(argv[optind++]);
	std::ifstream train_file(train_filename);
	if (train_file.fail()) {
		std::cerr << "File not found : " << train_filename << endl;
		return -1;
	}
	string test_filename(argv[optind++]);
	std::ifstream test_file(test_filename);
	if (test_file.fail()) {
		std::cerr << "File not found : " << test_filename << endl;
		return -1;
	}
	readData(train_file, test_file);
	cout << "train_file : " << train_filename << " , size " << db.gdata.num_train << endl;
	cout << "test_file : " << test_filename << " , size " << db.gdata.num_test << endl;

	db.ys.resize(db.gdata.num_train);

	if (setting.minsup == 0) {
		setting.minsup = 1;
	} else if (0 < setting.minsup and setting.minsup < 1) {
		setting.minsup = int(db.gdata.num_train * setting.minsup);
	}
	setting.print();

	db.gspan.minsup = setting.minsup;
	db.gspan.maxpat = setting.maxpat;
	db.gradient_boosting.run();

	std::cout << "\e[38;5;0m\e[48;5;40m --- end ---  \e[m" << std::endl; // debug
	return 0;
}
