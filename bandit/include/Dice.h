#pragma once

#include "MyInclude.h"
#include <random>

struct Dice {

	static std::random_device seed_gen;
	static std::mt19937 mt;
	static std::uniform_real_distribution<> dice;

	static bool p(double t) {
		assert(0 <= t and t <= 1);
		if (dice(mt) < t) {
			return true;
		} else {
			return false;
		}
	}

	static int id(int end) {
		static std::uniform_int_distribution<> dice_i(0, end);
		return dice_i(mt);
	}

};
