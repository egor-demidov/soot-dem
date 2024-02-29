//
// Created by egor on 2/29/24.
//

#include <random>

#include "break_neck.h"

static std::mt19937_64 mt(0);

void break_random_neck(std::vector<bool> & bonded_contacts, size_t n_part) {
    std::vector<std::pair<size_t, size_t>> necks;
    for (size_t i = 0; i < n_part - 1; i ++) {
        for (size_t j = i + 1; j < n_part; j ++) {
            bool bond = bonded_contacts[i * n_part + j];
            if (bond) necks.emplace_back(i, j);
        }
    }

    if (necks.empty())
        return;

    std::uniform_int_distribution<size_t> dist(0, necks.size() - 1);
    auto neck_to_break = necks[dist(mt)];
    bonded_contacts [neck_to_break.first * n_part + neck_to_break.second] = false;
    bonded_contacts [neck_to_break.second * n_part + neck_to_break.first] = false;
}
