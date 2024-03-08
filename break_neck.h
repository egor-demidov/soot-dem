//
// Created by egor on 2/29/24.
//

#ifndef SOOT_AFM_BREAK_NECK_H
#define SOOT_AFM_BREAK_NECK_H

#include <vector>

void break_random_neck(std::vector<bool> & bonded_contacts, size_t n_part);
void seed_rng(long seed);

#endif //SOOT_AFM_BREAK_NECK_H
