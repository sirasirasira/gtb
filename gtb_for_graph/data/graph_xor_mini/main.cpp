#include <iostream>
#include <vector>
#include <string>
#include <cassert>

using std::cin;
using std::cout;
using std::endl;
using std::vector;
using std::string;

vector<string> parts = {
	"111",
	"222",
	//"333",
	//"112",
	//"122",
	"121",
	"212",
	//"223",
	//"233",
	//"232",
	//"323",
	//"113",
	//"133",
	//"131",
	//"313",
	//"123",
	//"132",
	//"213",
};

int cnc_ptn_size = 9;

void printEdge(int from, int to) {
	int label = 0;
	cout << "e " << from << " " << to << " " << label << endl;
}

void connect(int cnc) {
	// 0 <= cnc < cnc_ptn_size
	assert(cnc_ptn_size == 9);
	assert(0 <= cnc);
	assert(cnc < 9);
	vector<int> from = {0, 1, 2};
	vector<int> to   = {3, 4, 5};
	int from_key = cnc / 3;
	int to_key = cnc % 3;
	printEdge(from[from_key], 6);
	printEdge(to[to_key], 6);
}

bool symmetrical(int key) {
	return parts[key][0] == parts[key][2]; 
}

inline bool either(int a, int b1, int b2, int b3) {
	return (a == b1) or (a == b2) or (a == b3);

}

void printGraph(int kp1, int kp2, int cnc) { // kp1 = key parts 1, cnc = connect
	assert(cnc_ptn_size == 9);
	if ((symmetrical(kp1) and either(cnc, 6, 7, 8)) or (symmetrical(kp2) and either(cnc, 2, 5, 8))) {
		return;
	}
	if (kp1 == kp2 and either(cnc, 3, 6, 7)) {
		return;
	}
	int res = ((kp1 + kp2) % 2 == 1) ? +1 : -1; // as response variable
	if (res > 0) {
	} else {
	}
	cout << "t # 0 " << res << endl;
	int vid = 0;
	for (auto label : parts[kp1]) {
		cout << "v " << vid << " " << label << endl;
		vid++;
	}
	for (auto label : parts[kp2]) {
		cout << "v " << vid << " " << label << endl;
		vid++;
	}
	cout << "v " << vid << " " << 4 << endl;
	vid++;
	printEdge(0, 1);
	printEdge(1, 2);
	printEdge(3, 4);
	printEdge(4, 5);
	connect(cnc);
	cout << endl;
}

int main() {

	bool swt = false; // used for bucket 
	for (size_t i = 0; i < parts.size(); i++) {
		for (size_t j = i; j < parts.size(); j++) {
			for (int k = 0; k < cnc_ptn_size; k++) {
				if (true or swt) { // true for all
					printGraph(i, j, k);
				}
				swt = !swt;
			}
		}
	}

	//std::cout << "\e[38;5;0m\e[48;5;40m --- end ---  \e[m" << std::endl; // debug
	return 0;
}
