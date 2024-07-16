/****************************************/
/*** Copyright (c) 2024, Egor Demidov ***/
/****************************************/

#include "random_engine.h"

static std::mt19937_64 mt(0);

std::mt19937_64 & get_random_engine() {
    return mt;
}

void seed_random_engine(long seed) {
    mt.seed(seed);
}
