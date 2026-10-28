//
// Created by egor on 10/27/26.
//

#ifndef SOOT_AFM_DEM_LB_PARAMETERS_H
#define SOOT_AFM_DEM_LB_PARAMETERS_H

#include <Eigen/Dense>
#include <cstddef>
#include <functional>
#include <vector>

struct DemLbParameters {
    using Vec = Eigen::Vector3d;
    using Field = std::vector<Vec>;

    // Identical spherical particles for now.
    double particleRadiusSI;
    double particleMassSI;
    double particleInertiaSI;
    double demTimeStepSI;

    double liquidDensitySI;
    double liquidDynamicViscositySI;
    double liquidSurfaceTensionSI;
    double liquidContactAngleDeg;

    unsigned cellsPerParticleDiameter;
    std::size_t particleCount = 1;

    // One force and torque per particle, in particle-index order.
    std::function<void(const Field &forcesSI, const Field &torquesSI)> submitFluidLoadSI;

    std::function<Vec(std::size_t)> particlePositionSI;
    std::function<Vec(std::size_t)> particleVelocitySI;
    std::function<Vec(std::size_t)> particleAngularVelocitySI;

    // Advances the entire DEM system once.
    std::function<void()> advanceDemStep;

    // Applies a remapping impulse to one particle.
    std::function<void(std::size_t particleIndex, const Vec &impulseSI, const Vec &angularImpulseSI)>
        applyRemapImpulseSI;
};

int runCapillaryLb(int argc, char **argv, const DemLbParameters &dem);

#endif // SOOT_AFM_DEM_LB_PARAMETERS_H
