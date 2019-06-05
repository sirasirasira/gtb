#include <iostream>

using namespace std;

void graph(int id) {
	int label = -1;
	if (id <= 3 or 12 <= id) label = 1;
	cout << "t # " << id << " " << label << endl;

	cout << "v 0 1" << endl;
	cout << "v 1 1" << endl;
	cout << "v 2 0" << endl;
	label = 2;
	if ((4 <= id and id <= 7) or (12 <= id)) label = 3;
	cout << "v 3 " << label << endl;
	cout << "v 4 " << label << endl;
	label = id % 4;
	cout << "v 5 " << label << endl;

	label = 0;
	if (8 <= id) label = 1;
	cout << "e 0 1 " << label << endl;
	cout << "e 1 2 0" << endl;
	cout << "e 2 3 0" << endl;
	cout << "e 2 4 0" << endl;
	cout << "e 3 4 0" << endl;
	cout << "e 4 5 0" << endl;

	cout << endl;
}

int main() {

	for (int i = 0; i < 16; i++) {
		graph(i);
	}

	//std::cout << "\e[38;5;0m\e[48;5;40m --- end ---  \e[m" << std::endl; // debug
	return 0;
}
