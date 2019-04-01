#pragma once

#include "MyInclude.h"

class Debug {
	public:
		static void IDs(vector<ID> vec, string str = "") {
			std::cout << "Debug " << str << " IDs : "; // debug
			for (auto id : vec) {
				cout << id << " ";
			}
			cout << endl;
		}

		static void ys(vector<double> ys, string str = "") {
			std::cout << "Debug " << str << " ys : "; // debug
			for (auto y : ys) {
				cout << y << " ";
			}
			cout << endl;
		}

};
