//
// Created by egor on 10/27/26.
//

#include "core/Environment.h"
#include <iostream>

int main(int argc, char** argv)
{
    walberla::Environment environment(argc, argv);
    std::cout << "waLBerla initialized\n";
    return 0;
}
