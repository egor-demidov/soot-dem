#pragma once
#include "capillary_setup.h"
inline std::filesystem::path writeCapillaryPrm(const CapillarySetup &c, unsigned long long steps) {
    std::string text = R"PRM(BlockForestParameters
{
   cellsPerBlock < @NX@, @NY@, @NZ@ >;
   periodicity < 0, 0, 0 >;
   loadBalancingFrequency        0;
   printLoadBalancingStatistics  false;
}

DomainParameters
{
   dropDiameter 1;
   domainSizeFactor < @NX@, @NY@, @NZ@ >;    // values multiplied with dropDiameter
}

PhysicsParameters
{
   relaxationRate 1.96;
   surfaceTension 8.37e-3;
   force <0, 0, 0>;
   enableWetting true;
   timesteps @STEPS@;
   diagnosticFrequency @LOG@;
}

ModelParameters
{
   pdfReconstructionModel        OnlyMissing;
   pdfRefillingModel             EquilibriumRefilling;
   excessMassDistributionModel   EvenlyAllInterfaceFallbackLiquid;
   curvatureModel                FiniteDifferenceMethod;
   useSimpleMassExchange         false;
   enableBubbleModel             false;
   enableBubbleSplits            false; // only used if enableBubbleModel=true
   cellConversionThreshold       1e-2;
   cellConversionForceThreshold  1e-1;
}

EvaluationParameters
{
   performanceLogFrequency 1000;
   evaluationFrequency 200;
   convergenceThreshold 0;
}

BoundaryParameters
{
   // X
   Border { direction W;  walldistance -1; NoSlip{} }
   Border { direction E;  walldistance -1; NoSlip{} }

   // Y
   Border { direction N;  walldistance -1; NoSlip{} }
   Border { direction S;  walldistance -1; NoSlip{} }

   // Z
   Border { direction T;  walldistance -1; NoSlip{} }
   Border { direction B;  walldistance -1; NoSlip{} }
}

MeshOutputParameters
{
   writeFrequency @MESH@;
   baseFolder mesh-out;
}

VTK
{
   fluid_field
   {
      writeFrequency @VTK@;
      ghostLayers 0;
      baseFolder vtk-out;
      samplingResolution 1;

      writers
      {
         fill_level;
         mapped_flag;
         velocity;
         density;
         curvature;
         normal;
         obstacle_normal;
         //pdf;
         //flag;
         //force_density;
      }

      inclusion_filters
      {
         // only include liquid and interface cells in VTK output
         //liquidInterfaceFilter;
      }

      before_functions
      {
         //ghost_layer_synchronization;   // only needed if writing the ghost layer
         gas_cell_zero_setter;            // sets velocity=0 and density=1 all gas cells before writing VTK output
      }
   }

   domain_decomposition
   {
      writeFrequency @VTK@;
      baseFolder vtk-out;
      outputDomainDecomposition true;
   }
})PRM";
    auto replace = [&](const std::string &from, const std::string &to) {
        std::size_t p = 0;
        while ((p = text.find(from, p)) != std::string::npos) {
            text.replace(p, from.size(), to);
            p += to.size();
        }
    };
    replace("@NX@", std::to_string(c.cells[0]));
    replace("@NY@", std::to_string(c.cells[1]));
    replace("@NZ@", std::to_string(c.cells[2]));
    replace("@STEPS@", std::to_string(steps));
    replace("@LOG@", std::to_string(c.logFrequency));
    replace("@MESH@", std::to_string(c.meshFrequency));
    replace("@VTK@", std::to_string(c.vtkFrequency));
    std::filesystem::create_directories("run");
    const std::filesystem::path path = "run/resolved_lb.prm";
    std::ofstream out(path);
    out << text;
    out.close();
    if (!out)
        throw std::runtime_error("Cannot write resolved LB configuration");
    return path;
}
