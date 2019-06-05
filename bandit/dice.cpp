#include "Dice.h"
#define CLASS Dice

std::random_device CLASS::seed_gen;
//std::mt19937 CLASS::mt(seed_gen());
std::mt19937 CLASS::mt(222);
std::uniform_real_distribution<> CLASS::dice(0.0, 1.0);
