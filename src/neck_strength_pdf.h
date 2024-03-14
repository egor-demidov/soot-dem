//
// Created by egor on 3/13/24.
//

#ifndef SOOT_AFM_NECK_STRENGTH_PDF_H
#define SOOT_AFM_NECK_STRENGTH_PDF_H

#include <random>

template <typename real_t>
struct neck_strength_pdf_t {

    template<typename generator_t>
    real_t operator () (generator_t & generator) {
        std::uniform_real_distribution<real_t> uniform_dist(0.0, 1.0);
        real_t cdf = uniform_dist(generator);

    }
};

#endif //SOOT_AFM_NECK_STRENGTH_PDF_H
