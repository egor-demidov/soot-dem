//
// Created by egor on 10/27/26.
//

#ifndef SOOT_AFM_LB_UNITS_H
#define SOOT_AFM_LB_UNITS_H

#include "dem_lb_parameters.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

struct LbUnits {
    double dxSI;
    double dtSI;

    unsigned lbStepsPerDemStep;

    double viscosityLB;
    double relaxationRate;
    double surfaceTensionLB;

    double forceScaleSI;
    double torqueScaleSI;
};

inline LbUnits makeLbUnits(const DemLbParameters& dem)
{
    for (double value : {
        dem.particleRadiusSI,
        dem.demTimeStepSI,
        dem.liquidDensitySI,
        dem.liquidDynamicViscositySI
    })
    {
        if (!std::isfinite(value) || value <= 0.0)
        {
            throw std::invalid_argument(
                "Radius, DEM timestep, liquid density and "
                "viscosity must be finite and positive"
            );
        }
    }

    if (dem.cellsPerParticleDiameter == 0 ||
        !std::isfinite(dem.liquidSurfaceTensionSI) ||
        dem.liquidSurfaceTensionSI < 0.0)
    {
        throw std::invalid_argument(
            "Invalid LB resolution or surface tension"
        );
    }

    LbUnits units{};

    units.dxSI =
        2.0 * dem.particleRadiusSI
        / dem.cellsPerParticleDiameter;

    const double dx2 = units.dxSI * units.dxSI;

    const double viscositySI =
        dem.liquidDynamicViscositySI
        / dem.liquidDensitySI;

    // Numerical target corresponding to SRT omega = 1.
    const double targetViscosityLB = 1.0 / 6.0;

    const double targetDtSI =
        targetViscosityLB * dx2 / viscositySI;

    // Divide each DEM step into an integer number of LB steps.
    const double substeps = std::max(
        1.0,
        std::ceil(dem.demTimeStepSI / targetDtSI)
    );

    if (!std::isfinite(substeps) ||
        substeps > std::numeric_limits<unsigned>::max())
    {
        throw std::invalid_argument(
            "LB substep count is too large"
        );
    }

    units.lbStepsPerDemStep =
        static_cast<unsigned>(substeps);

    units.dtSI =
        dem.demTimeStepSI / units.lbStepsPerDemStep;

    units.viscosityLB =
        viscositySI * units.dtSI / dx2;

    units.relaxationRate =
        1.0 / (0.5 + 3.0 * units.viscosityLB);

    // Reference liquid density in the LB solver is 1.
    const double massScaleSI =
        dem.liquidDensitySI * dx2 * units.dxSI;

    const double dt2 = units.dtSI * units.dtSI;

    units.surfaceTensionLB =
        dem.liquidSurfaceTensionSI * dt2 / massScaleSI;

    units.forceScaleSI =
        massScaleSI * units.dxSI / dt2;

    units.torqueScaleSI =
        units.forceScaleSI * units.dxSI;

    return units;
}

#endif //SOOT_AFM_LB_UNITS_H
