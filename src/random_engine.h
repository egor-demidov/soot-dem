/****************************************/
/*** Copyright (c) 2024, Egor Demidov ***/
/****************************************/

#ifndef SOOT_AFM_RANDOM_ENGINE_H
#define SOOT_AFM_RANDOM_ENGINE_H

#include <random>

std::mt19937_64 & get_random_engine();
void seed_random_engine(long seed);

#endif //SOOT_AFM_RANDOM_ENGINE_H
