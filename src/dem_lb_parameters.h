//
// Created by egor on 10/27/26.
//

#ifndef SOOT_AFM_DEM_LB_PARAMETERS_H
#define SOOT_AFM_DEM_LB_PARAMETERS_H

#include <Eigen/Dense>
#include <functional>

struct DemLbParameters {
    double particleRadiusSI;
    double particleMassSI;
    double particleInertiaSI;
    double demTimeStepSI;

    double liquidDensitySI;
    double liquidDynamicViscositySI;
    double liquidSurfaceTensionSI;
    double liquidContactAngleDeg;

    unsigned cellsPerParticleDiameter;

    std::function<void(
        const Eigen::Vector3d& forceSI,
        const Eigen::Vector3d& torqueSI
    )> submitFluidLoadSI;
    std::function<Eigen::Vector3d()> particlePositionSI;
    std::function<Eigen::Vector3d()> particleVelocitySI;
    std::function<Eigen::Vector3d()> particleAngularVelocitySI;

    std::function<void()> advanceDemStep;

    std::function<void(
        const Eigen::Vector3d&,
        const Eigen::Vector3d&
    )> applyRemapImpulseSI;
};

int runCapillaryLb(
    int argc,
    char** argv,
    const DemLbParameters& dem
);

#endif //SOOT_AFM_DEM_LB_PARAMETERS_H
