//
// Created by egor on 10/27/26.
//

//======================================================================================================================
//
//  This file is part of waLBerla. waLBerla is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  waLBerla is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with waLBerla (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \file DropWetting.cpp
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//
// This showcase simulates a droplet on a plane with specified target contact angle. The resulting drop height can be
// compared to an analytical model, as done in e.g.:
// dissertation of S. Bogner (2017), section 4.4.3.2 , urn:nbn:de:bvb:29-opus4-87191
//======================================================================================================================

#include <utility>
#include <array>
#include <queue>
#include <set>
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>

#include "blockforest/Initialization.h"

#include "core/Environment.h"

#include "lbm/PerformanceLogger.h"
#include "lbm/field/AddToStorage.h"
#include "lbm/free_surface/LoadBalancing.h"
#include "lbm/free_surface/SurfaceMeshWriter.h"
#include "lbm/free_surface/TotalMassComputer.h"
#include "lbm/free_surface/VtkWriter.h"
#include "lbm/free_surface/bubble_model/Geometry.h"
#include "lbm/free_surface/dynamics/SurfaceDynamicsHandler.h"
#include "lbm/free_surface/surface_geometry/SurfaceGeometryHandler.h"
#include "lbm/free_surface/surface_geometry/Utility.h"
#include "lbm/lattice_model/D3Q19.h"

#include "contact_line_force.h"
#include "sphere_surface.h"
#include "interface_plane_field.h"
#include "surface_contact_line.h"
#include "contact_line_output.h"
#include <sstream>

#include "dem_lb_parameters.h"
#include "lb_units.h"

#include "audited_surface_dynamics_handler.h"

namespace walberla::free_surface
{
   template< typename ScalarField_T, typename FlagField_T >
class CapillarySurfaceMeshWriter
{
 public:
   CapillarySurfaceMeshWriter(const std::weak_ptr< StructuredBlockForest >& blockForest, const ConstBlockDataID& fillFieldID,
                     const ConstBlockDataID& flagFieldID, const Set< FlagUID >& liquidInterfaceGasFlagIDSet,
                     real_t obstacleFillLevel, uint_t writeFrequency, const std::string& baseFolder)
      : blockForest_(blockForest), fillFieldID_(fillFieldID), flagFieldID_(flagFieldID),
        liquidInterfaceGasFlagIDSet_(liquidInterfaceGasFlagIDSet), obstacleFillLevel_(obstacleFillLevel),
        writeFrequency_(writeFrequency), baseFolder_(baseFolder), executionCounter_(uint_c(0)),
        fillFieldCloneID_(
           (abortIfNullptr::abortIfNullptr(blockForest),
            field::addCloneToStorage< ScalarField_T >(blockForest_.lock(), fillFieldID_, "Fill level field clone")))
   {}

   // config block must be named "MeshOutputParameters"
   CapillarySurfaceMeshWriter(const std::weak_ptr< StructuredBlockForest >& blockForest, const ConstBlockDataID& fillFieldID,
                     const ConstBlockDataID& flagFieldID, const Set< FlagUID >& liquidInterfaceGasFlagIDSet,
                     real_t obstacleFillLevel, const std::weak_ptr< Config >& config)
      : CapillarySurfaceMeshWriter(
           blockForest, fillFieldID, flagFieldID, liquidInterfaceGasFlagIDSet, obstacleFillLevel,
           (abortIfNullptr::abortIfNullptr(config),
            config.lock()->getOneBlock("MeshOutputParameters").getParameter< uint_t >("writeFrequency")),
           (abortIfNullptr::abortIfNullptr(config),
            config.lock()->getOneBlock("MeshOutputParameters").getParameter< std::string >("baseFolder")))
   {}

   void operator()()
   {
      if (writeFrequency_ == uint_c(0)) { return; }

      if (executionCounter_ == uint_c(0))
      {
         createBaseFolder();
         writeMesh();
      }
      else { writeMesh(); }

      ++executionCounter_;
   }

 private:
   void createBaseFolder() const
   {
      WALBERLA_ROOT_SECTION()
      {
         const std::filesystem::path basePath(baseFolder_);
         if (std::filesystem::exists(basePath)) { std::filesystem::remove_all(basePath); }
         std::filesystem::create_directories(basePath);
      }
      WALBERLA_MPI_BARRIER();
   }

   void writeMesh()
   {
      // only write mesh in given frequency
      if (executionCounter_ % writeFrequency_ != uint_c(0)) { return; }

      // rank=0 is just an arbitrary choice here
      const int targetRank = 0;

      auto blockForest = blockForest_.lock();
      WALBERLA_CHECK_NOT_NULLPTR(blockForest);

      // update clone of fill level field and set fill level of all non-liquid, -interface, or -gas cells to zero
      updateFillFieldClone(blockForest);

      const auto surfaceMesh = postprocessing::realFieldToSurfaceMesh< ScalarField_T >(
         blockForest, fillFieldCloneID_, real_c(0.5), uint_c(0), true, targetRank, MPI_COMM_WORLD);

      surfaceMesh->translate(
    Vector3<real_t>(
        real_c(0.5) * blockForest->dx(),
        real_c(0.5) * blockForest->dy(),
        real_c(0.5) * blockForest->dz()
    )
);

      WALBERLA_EXCLUSIVE_WORLD_SECTION(targetRank)
      {
         std::ostringstream filename;
         filename << "simulation_step_"
                  << std::setfill('0')
                  << std::setw(12)
                  << executionCounter_
                  << ".vtp";

         const auto path =
             std::filesystem::path(baseFolder_) / filename.str();

         std::ofstream output(path);
         if (!output.is_open())
         {
            WALBERLA_ABORT("Cannot open surface mesh output file.");
         }

         geometry::writeMeshVtp(output, *surfaceMesh);
      }
   }

   // update clone of fill level field and set fill level of all non-liquid, -interface, or -gas cells to zero;
   // explicitly use shared_ptr instead of weak_ptr to avoid checking the latter's validity (is done in writeMesh()
   // already)
   void updateFillFieldClone(
    const std::shared_ptr<StructuredBlockForest>& blockForest)
{
    for (auto block = blockForest->begin();
         block != blockForest->end(); ++block)
    {
        auto* output =
            block->template getData<ScalarField_T>(fillFieldCloneID_);

        const auto* input =
            block->template getData<ScalarField_T>(fillFieldID_);

        const auto* flags =
            block->template getData<FlagField_T>(flagFieldID_);

        const auto fluidMask =
            flags->getMask(liquidInterfaceGasFlagIDSet_);

        const cell_idx_t nx = cell_idx_c(input->xSize());
        const cell_idx_t ny = cell_idx_c(input->ySize());
        const cell_idx_t nz = cell_idx_c(input->zSize());

        const cell_idx_t inputGhost =
            cell_idx_c(input->nrOfGhostLayers());
        const cell_idx_t flagGhost =
            cell_idx_c(flags->nrOfGhostLayers());

        const cell_idx_t ng =
            inputGhost < flagGhost ? inputGhost : flagGhost;

        for (const auto& cell : output->xyzSizeWithGhostLayer())
        {
            real_t sum = real_c(0);
            uint_t validCount = uint_c(0);

            // Value at a grid corner: average its eight adjacent cells.
            for (cell_idx_t dz = -1; dz <= 0; ++dz)
            for (cell_idx_t dy = -1; dy <= 0; ++dy)
            for (cell_idx_t dx = -1; dx <= 0; ++dx)
            {
                const cell_idx_t x = cell.x() + dx;
                const cell_idx_t y = cell.y() + dy;
                const cell_idx_t z = cell.z() + dz;

                if (x < -ng || x >= nx + ng ||
                    y < -ng || y >= ny + ng ||
                    z < -ng || z >= nz + ng)
                {
                    continue;
                }

                // Exclude solid bookkeeping fill values from the mesh.
                if ((flags->get(x, y, z) & fluidMask) != 0)
                {
                    sum += input->get(x, y, z);
                    ++validCount;
                }
            }

            output->get(cell.x(), cell.y(), cell.z()) =
                validCount > uint_c(0)
                    ? sum / real_c(validCount)
                    : real_c(0);
        }
    }
}

   std::weak_ptr< StructuredBlockForest > blockForest_;
   ConstBlockDataID fillFieldID_;
   ConstBlockDataID flagFieldID_;
   Set< FlagUID > liquidInterfaceGasFlagIDSet_;
   real_t obstacleFillLevel_;
   uint_t writeFrequency_;
   std::string baseFolder_;
   uint_t executionCounter_;

   BlockDataID fillFieldCloneID_;
}; // class SurfaceMeshWriter
}

namespace walberla
{
namespace free_surface
{
namespace DropWetting
{
using ScalarField_T = GhostLayerField< real_t, 1 >;
using VectorField_T = GhostLayerField< Vector3< real_t >, 1 >;

using CollisionModel_T      = lbm::collision_model::SRT;
using ForceModel_T          = lbm::force_model::GuoField< VectorField_T >;
using LatticeModel_T        = lbm::D3Q19< CollisionModel_T, true, ForceModel_T, 2 >;
using LatticeModelStencil_T = LatticeModel_T::Stencil;
using PdfField_T            = lbm::PdfField< LatticeModel_T >;
using PdfCommunication_T    = blockforest::SimpleCommunication< LatticeModelStencil_T >;

// the geometry computations in SurfaceGeometryHandler require meaningful values in the ghost layers in corner
// directions (flag field and fill level field); this holds, even if the lattice model uses a D3Q19 stencil
using CommunicationStencil_T =
   std::conditional_t< LatticeModel_T::Stencil::D == uint_t{2}, stencil::D2Q9, stencil::D3Q27 >;
using Communication_T = blockforest::SimpleCommunication< CommunicationStencil_T >;

using flag_t                        = uint32_t;
using FlagField_T                   = FlagField< flag_t >;
using FreeSurfaceBoundaryHandling_T = FreeSurfaceBoundaryHandling< LatticeModel_T, FlagField_T, ScalarField_T >;

// compute the L2 norm of the relative difference between the current and the previous time-averaged quantities
template< typename FreeSurfaceBoundaryHandling_T >
class DropHeightEvaluator
{
 public:
   DropHeightEvaluator(const std::weak_ptr< const StructuredBlockForest >& blockForest,
                       const std::weak_ptr< const FreeSurfaceBoundaryHandling_T >& boundaryHandlingID,
                       const ConstBlockDataID& fillFieldID, const std::shared_ptr< real_t >& dropHeight,
                       uint_t interval)
      : blockForest_(blockForest), freeSurfaceBoundaryHandling_(boundaryHandlingID), fillFieldID_(fillFieldID),
        dropHeight_(dropHeight), interval_(interval), executionCounter_(uint_c(0))
   {}

   void operator()()
   {
      ++executionCounter_;

      // only evaluate in given intervals
      if (executionCounter_ % interval_ != uint_c(0)) { return; }

      // compute the height of the spherical drop cap
      computeDropHeight();
   }

 private:
   std::weak_ptr< const StructuredBlockForest > blockForest_;
   std::weak_ptr< const FreeSurfaceBoundaryHandling_T > freeSurfaceBoundaryHandling_;
   ConstBlockDataID fillFieldID_;

   std::shared_ptr< real_t > dropHeight_;

   uint_t interval_;

   uint_t executionCounter_; // number of times operator() has been called

   void computeDropHeight()
   {
      auto blockForest = blockForest_.lock();
      WALBERLA_CHECK_NOT_NULLPTR(blockForest);

      auto freeSurfaceBoundaryHandling = freeSurfaceBoundaryHandling_.lock();
      WALBERLA_CHECK_NOT_NULLPTR(freeSurfaceBoundaryHandling);

      const BlockDataID flagFieldID = freeSurfaceBoundaryHandling->getFlagFieldID();
      const typename FreeSurfaceBoundaryHandling_T::FlagInfo_T& flagInfo = freeSurfaceBoundaryHandling->getFlagInfo();

      *dropHeight_ = real_c(0);

      for (auto blockIt = blockForest->begin(); blockIt != blockForest->end(); ++blockIt)
      {
         real_t maxDropHeight = real_c(0);

         const FlagField_T* const flagField   = blockIt->template getData< const FlagField_T >(flagFieldID);
         const ScalarField_T* const fillField = blockIt->template getData< const ScalarField_T >(fillFieldID_);

         WALBERLA_FOR_ALL_CELLS(
            flagFieldIt, flagField, fillFieldIt, fillField,

            if (flagInfo.isInterface(flagFieldIt)) {
               Cell globalCell = flagFieldIt.cell();
               blockForest->transformBlockLocalToGlobalCell(globalCell, *blockIt);

               real_t currentDropHeight = real_c(globalCell.z()) + *fillFieldIt;
               maxDropHeight            = currentDropHeight > maxDropHeight ? currentDropHeight : maxDropHeight;
            }) // WALBERLA_FOR_ALL_CELLS

         if (maxDropHeight > *dropHeight_) { *dropHeight_ = maxDropHeight; }
      }
      // get maximum dropHeight among all processes
      mpi::allReduceInplace< real_t >(*dropHeight_, mpi::MAX);

      // WALBERLA_LOG_DEVEL_VAR_ON_ROOT(*dropHeight_);
   }
}; // class DropHeightEvaluator

int runWithDem(int argc, char** argv, const DemLbParameters& dem)
{
   Environment walberlaEnv(argc, argv);

   if (argc < 2) { WALBERLA_ABORT("Please specify a parameter file as input argument.") }

    const LbUnits units = makeLbUnits(dem);

    if (!dem.particlePositionSI)
    {
        throw std::runtime_error(
            "DEM particle-position callback is missing"
        );
    }

    const Eigen::Vector3d demReferencePositionSI =
        dem.particlePositionSI();

    if (!dem.particleVelocitySI ||
        !dem.particleAngularVelocitySI)
    {
        throw std::runtime_error(
            "DEM velocity callbacks are missing"
        );
    }

    const Eigen::Vector3d initialVelocityLB =
        dem.particleVelocitySI()
        * (units.dtSI / units.dxSI);

    const Eigen::Vector3d initialAngularVelocityLB =
        dem.particleAngularVelocitySI()
        * units.dtSI;

    WALBERLA_LOG_INFO_ON_ROOT(
        "DEM initial velocity [m/s] = "
        << dem.particleVelocitySI().transpose()
        << ", velocity LB = "
        << initialVelocityLB.transpose()
        << ", angular velocity LB = "
        << initialAngularVelocityLB.transpose()
    );

    WALBERLA_LOG_INFO_ON_ROOT(
        "Derived LB scales:"
        << " dxSI=" << units.dxSI
        << " dtSI=" << units.dtSI
        << " LBstepsPerDEMstep=" << units.lbStepsPerDemStep
        << " viscosityLB=" << units.viscosityLB
        << " omega=" << units.relaxationRate
        << " surfaceTensionLB=" << units.surfaceTensionLB
        << " forceScaleSI=" << units.forceScaleSI
        << " torqueScaleSI=" << units.torqueScaleSI
    );

   // get block forest parameters from parameter file
   auto blockForestParameters              = walberlaEnv.config()->getOneBlock("BlockForestParameters");
   const Vector3< uint_t > cellsPerBlock   = blockForestParameters.getParameter< Vector3< uint_t > >("cellsPerBlock");
   const Vector3< bool > periodicity       = blockForestParameters.getParameter< Vector3< bool > >("periodicity");
   const uint_t loadBalancingFrequency     = blockForestParameters.getParameter< uint_t >("loadBalancingFrequency");
   const bool printLoadBalancingStatistics = blockForestParameters.getParameter< bool >("printLoadBalancingStatistics");

   // read domain parameters from parameter file
   const auto domainParameters              = walberlaEnv.config()->getOneBlock("DomainParameters");
   const real_t dropDiameter                = domainParameters.getParameter< real_t >("dropDiameter");
   const Vector3< real_t > domainSizeFactor = domainParameters.getParameter< Vector3< real_t > >("domainSizeFactor");

   // define domain size
   Vector3< uint_t > domainSize = domainSizeFactor * dropDiameter;
   domainSize[0]                = uint_c(domainSizeFactor[0] * dropDiameter);
   domainSize[1]                = uint_c(domainSizeFactor[1] * dropDiameter);
   domainSize[2]                = uint_c(domainSizeFactor[2] * dropDiameter);

   // compute number of blocks as defined by domainSize and cellsPerBlock
   Vector3< uint_t > numBlocks;
   numBlocks[0] = uint_c(std::ceil(real_c(domainSize[0]) / real_c(cellsPerBlock[0])));
   numBlocks[1] = uint_c(std::ceil(real_c(domainSize[1]) / real_c(cellsPerBlock[1])));
   numBlocks[2] = uint_c(std::ceil(real_c(domainSize[2]) / real_c(cellsPerBlock[2])));

   // get number of (MPI) processes
   uint_t numProcesses = uint_c(MPIManager::instance()->numProcesses());
   WALBERLA_CHECK_LESS_EQUAL(numProcesses, numBlocks[0] * numBlocks[1] * numBlocks[2],
                             "The number of MPI processes is greater than the number of blocks as defined by "
                             "\"domainSize/cellsPerBlock\". This would result in unused MPI processes. Either decrease "
                             "the number of MPI processes or increase \"cellsPerBlock\".")

   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(numProcesses);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(cellsPerBlock);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(numBlocks);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(dropDiameter);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(domainSizeFactor);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(periodicity);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(domainSize);

   Vector3< uint_t > realDomainSize;
   realDomainSize[0] = cellsPerBlock[0] * numBlocks[0];
   realDomainSize[1] = cellsPerBlock[1] * numBlocks[1];
   realDomainSize[2] = cellsPerBlock[2] * numBlocks[2];

   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(realDomainSize);

   if (domainSize[0] != realDomainSize[0] && periodicity[0])
   {
      WALBERLA_ABORT(
         "The specified domain size in x-direction can not be obtained with the number of blocks you specified.")
   }
   if (domainSize[1] != realDomainSize[1] && periodicity[1])
   {
      WALBERLA_ABORT(
         "The specified domain size in y-direction can not be obtained with the number of blocks you specified.")
   }
   // the z-direction must be no slip in this setup and is simply extended (see below) to obtain the specified size
   int boundaryThicknessZ = int_c(realDomainSize[2]) - int_c(domainSize[2]);
   if (boundaryThicknessZ < 0)
   {
      WALBERLA_ABORT("Something went wrong: the resulting domain size in z-direction is less than specified.")
   }

   // read physics parameters from parameter file
   const auto physicsParameters  = walberlaEnv.config()->getOneBlock("PhysicsParameters");
    const real_t relaxationRate =
     real_c(units.relaxationRate);
    const real_t surfaceTension =
     real_c(units.surfaceTensionLB);
   const Vector3< real_t > force = physicsParameters.getParameter< Vector3< real_t > >("force");

   const bool enableWetting  = physicsParameters.getParameter< bool >("enableWetting");
    const real_t contactAngle =
     real_c(dem.liquidContactAngleDeg);

    WALBERLA_LOG_INFO_ON_ROOT(
    "Effective contact angle [deg] = " << contactAngle
    << ", enableWetting = " << enableWetting);

   const CollisionModel_T collisionModel = CollisionModel_T(relaxationRate);
   const real_t viscosity                = collisionModel.viscosity();

   const uint_t timesteps = physicsParameters.getParameter< uint_t >("timesteps");

   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(relaxationRate);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(surfaceTension);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(force);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(enableWetting);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(contactAngle);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(timesteps);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(viscosity);

   // read model parameters from parameter file
   const auto modelParameters               = walberlaEnv.config()->getOneBlock("ModelParameters");
   const std::string pdfReconstructionModel = modelParameters.getParameter< std::string >("pdfReconstructionModel");
   const std::string pdfRefillingModel      = modelParameters.getParameter< std::string >("pdfRefillingModel");
   const std::string excessMassDistributionModel =
      modelParameters.getParameter< std::string >("excessMassDistributionModel");
   const std::string curvatureModel          = modelParameters.getParameter< std::string >("curvatureModel");
   const bool useSimpleMassExchange          = modelParameters.getParameter< bool >("useSimpleMassExchange");
   const bool enableBubbleModel              = modelParameters.getParameter< bool >("enableBubbleModel");
   const bool enableBubbleSplits             = modelParameters.getParameter< bool >("enableBubbleSplits");
   const real_t cellConversionThreshold      = modelParameters.getParameter< real_t >("cellConversionThreshold");
   const real_t cellConversionForceThreshold = modelParameters.getParameter< real_t >("cellConversionForceThreshold");

   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(pdfReconstructionModel);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(pdfRefillingModel);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(excessMassDistributionModel);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(curvatureModel);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(useSimpleMassExchange);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(enableBubbleModel);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(enableBubbleSplits);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(cellConversionThreshold);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(cellConversionForceThreshold);

   // read evaluation parameters from parameter file
   const auto evaluationParameters      = walberlaEnv.config()->getOneBlock("EvaluationParameters");
   const uint_t performanceLogFrequency = evaluationParameters.getParameter< uint_t >("performanceLogFrequency");
   const uint_t evaluationFrequency     = evaluationParameters.getParameter< uint_t >("evaluationFrequency");
   const real_t convergenceThreshold    = evaluationParameters.getParameter< real_t >("convergenceThreshold");

   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(performanceLogFrequency);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(evaluationFrequency);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(convergenceThreshold);

   // create non-uniform block forest (non-uniformity required for load balancing)
   const std::shared_ptr< StructuredBlockForest > blockForest =
      createNonUniformBlockForest(domainSize, cellsPerBlock, numBlocks, periodicity);

   // add force field
   const BlockDataID forceDensityFieldID =
      field::addToStorage< VectorField_T >(blockForest, "Force field", force, field::fzyx, uint_c(1));

   // create lattice model
   const LatticeModel_T latticeModel = LatticeModel_T(collisionModel, ForceModel_T(forceDensityFieldID));

   // add pdf field
   const BlockDataID pdfFieldID = lbm::addPdfFieldToStorage(blockForest, "PDF field", latticeModel, field::fzyx);

   // add fill level field (initialized with 0, i.e., gas everywhere)
   const BlockDataID fillFieldID =
      field::addToStorage< ScalarField_T >(blockForest, "Fill level field", real_c(0.0), field::fzyx, uint_c(2));

    // Fixed solid sphere, resting on the bottom plane.
    const real_t solidRadius =
     real_c(0.5) * real_c(dem.cellsPerParticleDiameter);

    // Map the initial DEM position onto the existing LB placement.
    Eigen::Vector3d lbAnchor(
        0.5 * double(domainSize[0]),
        0.5 * double(domainSize[1]),
        0.5 * double(domainSize[2]));

    auto demPositionToLb =
        [&](const Eigen::Vector3d& positionSI)
            -> Eigen::Vector3d
        {
            return lbAnchor
                 + (positionSI - demReferencePositionSI)
                   / units.dxSI;
        };

    Eigen::Vector3d solidCenterEigen =
        demPositionToLb(dem.particlePositionSI());

    Vector3<real_t> solidCenter(
        real_c(solidCenterEigen.x()),
        real_c(solidCenterEigen.y()),
        real_c(solidCenterEigen.z())
    );

   // add boundary handling
   const std::shared_ptr< FreeSurfaceBoundaryHandling_T > freeSurfaceBoundaryHandling =
      std::make_shared< FreeSurfaceBoundaryHandling_T >(blockForest, pdfFieldID, fillFieldID);
   const BlockDataID flagFieldID                                      = freeSurfaceBoundaryHandling->getFlagFieldID();
   const typename FreeSurfaceBoundaryHandling_T::FlagInfo_T& flagInfo = freeSurfaceBoundaryHandling->getFlagInfo();

   // add spherical drop to fill level field
    // Preserve the original droplet-to-particle displacement.
    const real_t dropOffsetZ =
        real_c(1.1) * real_c(dropDiameter) - solidRadius;

    const geometry::Sphere sphereDrop(
        Vector3<real_t>(
            solidCenter[0],
            solidCenter[1],
            solidCenter[2] + dropOffsetZ),
        real_c(0.5) * real_c(dropDiameter));
   bubble_model::addBodyToFillLevelField< geometry::Sphere >(*blockForest, fillFieldID, sphereDrop, false);

   // initialize boundary conditions from config file
   const auto boundaryParameters = walberlaEnv.config()->getOneBlock("BoundaryParameters");
   freeSurfaceBoundaryHandling->initFromConfig(boundaryParameters);

   //// ADDED /////
    WALBERLA_LOG_INFO_ON_ROOT(
        "DEM reference position [m] = "
        << demReferencePositionSI.transpose()
        << ", LB particle center = "
        << solidCenter
    );

    auto solidSurface = makeSphereSurface(
        solidCenterEigen,
        double(solidRadius),
        4
    );

    struct InterfacePlaneSample {
        Eigen::Vector3d point;
        Eigen::Vector3d normal; // Points from liquid toward gas.
    };

    std::vector<InterfacePlaneSample> interfacePlanes;

    {
    using Vec = Eigen::Vector3d;

    const double radius = double(solidRadius);
    const double sigma = double(surfaceTension);
    const double pi = std::acos(-1.0);

    // Synthetic droplet used only to test reconstruction.
    // It does not change the simulated droplet.
    const Vec testDropletCenter =
        solidCenterEigen + Vec{0, 0, radius};

    ContactLineLoad total;

    double lineLength = 0.0;
    int segmentCount = 0;
    int degenerateCount = 0;

    for (const auto& triangle : solidSurface)
    {
        std::array<double, 3> values;

        for (int i = 0; i < 3; ++i)
        {
            values[i] =
                (triangle[i] - testDropletCenter).norm()
                - radius;
        }

        Vec a, b;

        const int result = cutSurfaceTriangle(
            triangle,
            values,
            1e-10 * radius,
            a,
            b
        );

        if (result < 0)
        {
            ++degenerateCount;
            continue;
        }

        if (result == 0)
            continue;

        const Vec solidNormal =
            (triangle[1] - triangle[0])
            .cross(triangle[2] - triangle[0])
            .normalized();

        const Vec midpoint = 0.5 * (a + b);
        const Vec tangent = (b - a).normalized();

        Vec interfaceNormal =
            (midpoint - testDropletCenter).normalized();

        // Align the interface normal with the discrete segment.
        interfaceNormal =
            (interfaceNormal
             - interfaceNormal.dot(tangent) * tangent)
            .normalized();

        const auto load = contactLineSegmentLoad(
            a,
            b,
            solidNormal,
            interfaceNormal,
            solidCenterEigen,
            sigma
        );

        total.force += load.force;
        total.torque += load.torque;

        lineLength += (b - a).norm();
        ++segmentCount;
    }

    const double expectedLength =
        std::sqrt(3.0) * pi * radius;

    const double expectedForceZ =
        1.5 * pi * radius * sigma;

    WALBERLA_LOG_INFO_ON_ROOT(
        "Smooth-surface check:"
        << " segments=" << segmentCount
        << " degenerate=" << degenerateCount
        << " length=" << lineLength
        << " expectedLength=" << expectedLength
        << " force=" << total.force.transpose()
        << " expectedFz=" << expectedForceZ
        << " torqueNorm=" << total.torque.norm()
    );
}

   using BoundaryHandling_T =
       FreeSurfaceBoundaryHandling_T::BoundaryHandling_T;
    using UBB_T =
        FreeSurfaceBoundaryHandling_T::UBB_T;

   for (auto blockIt = blockForest->begin();
        blockIt != blockForest->end();
        ++blockIt)
   {
      auto* handling =
          blockIt->getData<BoundaryHandling_T>(
              freeSurfaceBoundaryHandling->getHandlingID()
          );

      auto* flags = handling->getFlagField();

      // Include ghost cells so neighboring blocks see consistent geometry.
      for (const auto& cell : flags->xyzSizeWithGhostLayer())
      {
         const auto position =
             blockForest->getBlockLocalCellCenter(*blockIt, cell);

         if ((position - solidCenter).length() <= solidRadius)
         {
             handling->forceBoundary(
                 FreeSurfaceBoundaryHandling_T::ubbFlagID,
                 cell.x(), cell.y(), cell.z(),
                 UBB_T::Velocity(
                     real_c(0),
                     real_c(0),
                     real_c(0)
                 )
             );
         }
      }
   }
   ///////////////////////////////////////////////

   // IMPORTANT REMARK: this must be only called after every solid flag has been set; otherwise, the boundary handling
   // might not detect solid flags correctly
   freeSurfaceBoundaryHandling->initFlagsFromFillLevel();

    ConstBlockDataID remapExcessMassFieldID = ConstBlockDataID();

    struct ParticleChangedCell
    {
        IBlock* block;
        Cell cell;
    };

    struct CoveredLiquidCell
    {
        IBlock* block;
        Cell cell;

        double massLB;
        Eigen::Vector3d momentumLB;
    };

    struct ParticleCellChanges
    {
        uint_t covered = 0;
        uint_t uncovered = 0;
        uint_t retained = 0;

        std::vector<ParticleChangedCell> coveredCells;
        std::vector<ParticleChangedCell> uncoveredCells;

        double coveredLiquidVolumeLB = 0.0;
        double coveredLiquidMassLB = 0.0;

        std::vector<CoveredLiquidCell> coveredLiquidCells;

        Eigen::Vector3d coveredLiquidMomentumLB =
            Eigen::Vector3d::Zero();
    };

    auto detectParticleCellChanges =
        [&](const Eigen::Vector3d& proposedCenterLB)
        {
            ParticleCellChanges changes;

            const real_t radiusSquared = solidRadius * solidRadius;

            const Vector3<real_t> proposedCenter(
                real_c(proposedCenterLB.x()),
                real_c(proposedCenterLB.y()),
                real_c(proposedCenterLB.z()));

            for (auto& block : *blockForest)
            {
                const auto* flags =
                    block.getData<FlagField_T>(flagFieldID);

                const auto* fill =
                    block.getData<ScalarField_T>(fillFieldID);

                const auto* pdfs =
                    block.getData<PdfField_T>(pdfFieldID);

                const ScalarField_T* storedExcess = nullptr;

                if (remapExcessMassFieldID != ConstBlockDataID())
                {
                    storedExcess = block.getData<ScalarField_T>(
                        remapExcessMassFieldID);
                }

                const auto particleMask = flags->getFlag(
                    FreeSurfaceBoundaryHandling_T::ubbFlagID);

                // Count owned cells only; ghost cells would duplicate counts.
                for (const auto& cell : flags->xyzSize())
                {
                    const bool wasParticle =
                        (flags->get(cell.x(), cell.y(), cell.z())
                         & particleMask) != 0;

                    const auto position =
                        blockForest->getBlockLocalCellCenter(block, cell);

                    const auto displacement = position - proposedCenter;

                    const real_t distanceSquared =
                        displacement[0] * displacement[0] +
                        displacement[1] * displacement[1] +
                        displacement[2] * displacement[2];

                    const bool becomesParticle =
                        distanceSquared <= radiusSquared;

                    if (!wasParticle && becomesParticle)
                    {
                        ++changes.covered;
                        changes.coveredCells.push_back({&block, cell});

                        const auto cellFlag =
                            flags->get(cell.x(), cell.y(), cell.z());

                        if (flagInfo.isLiquid(cellFlag) ||
                            flagInfo.isInterface(cellFlag))
                        {
                            const double fraction = double(
                                fill->get(cell.x(), cell.y(), cell.z()));

                            const double density = double(
                                pdfs->getDensity(cell.x(), cell.y(), cell.z()));

                            const double extraMassLB =
                                storedExcess != nullptr
                                    ? double(storedExcess->get(
                                          cell.x(), cell.y(), cell.z()))
                                    : 0.0;

                            const double liquidMassLB =
                                fraction * density + extraMassLB;

                            // Momentum density of the underlying LB liquid.
                            const auto momentumDensity =
                                pdfs->getMomentumDensity(
                                    cell.x(), cell.y(), cell.z());

                            // Only the liquid-filled part contributes liquid momentum.
                            const Eigen::Vector3d liquidMomentumLB(
                                fraction * double(momentumDensity[0]),
                                fraction * double(momentumDensity[1]),
                                fraction * double(momentumDensity[2]));

                            if (!std::isfinite(liquidMassLB) ||
                                !liquidMomentumLB.allFinite())
                            {
                                throw std::runtime_error(
                                    "Non-finite liquid state in a covered cell");
                            }

                            changes.coveredLiquidCells.push_back({
                                &block,
                                cell,
                                liquidMassLB,
                                liquidMomentumLB
                            });

                            changes.coveredLiquidVolumeLB += fraction;
                            changes.coveredLiquidMassLB += liquidMassLB;
                            changes.coveredLiquidMomentumLB += liquidMomentumLB;
                        }
                    }
                    else if (wasParticle && !becomesParticle)
                    {
                        ++changes.uncovered;
                        changes.uncoveredCells.push_back({&block, cell});
                    }
                    else if (wasParticle && becomesParticle)
                    {
                        ++changes.retained;
                    }
                }
            }

            return changes;
        };

    struct RemapResult
    {
        Eigen::Vector3d linearImpulseLB;
        Eigen::Vector3d angularImpulseAboutOldCenterLB;
        Eigen::Vector3d oldCenterLB;
    };

    auto remapParticle =
        [&](const Eigen::Vector3d& shiftedCenter) -> RemapResult
        {
            if (!shiftedCenter.allFinite())
            {
                throw std::runtime_error(
                    "Non-finite proposed particle center");
            }

    const auto stationary =
    detectParticleCellChanges(solidCenterEigen);

    const auto shifted =
        detectParticleCellChanges(shiftedCenter);

    WALBERLA_LOG_INFO_ON_ROOT(
        "Stationary cell check:"
        << " covered=" << stationary.covered
        << " uncovered=" << stationary.uncovered
        << " retained=" << stationary.retained
    );

    WALBERLA_LOG_INFO_ON_ROOT(
        "One-cell translation check:"
        << " covered=" << shifted.covered
        << " uncovered=" << shifted.uncovered
        << " retained=" << shifted.retained
    );

    const double cellVolumeSI =
        units.dxSI * units.dxSI * units.dxSI;

    const double massUnitSI =
        dem.liquidDensitySI * cellVolumeSI;

    WALBERLA_LOG_INFO_ON_ROOT(
        "Proposed translation liquid:"
        << " volumeLB=" << shifted.coveredLiquidVolumeLB
        << " massLB=" << shifted.coveredLiquidMassLB
        << " volumeSI_m3="
        << shifted.coveredLiquidVolumeLB * cellVolumeSI
        << " massSI_kg="
        << shifted.coveredLiquidMassLB * massUnitSI
    );

    double recordedMassLB = 0.0;
    Eigen::Vector3d recordedMomentumLB =
        Eigen::Vector3d::Zero();

    for (const auto& entry : shifted.coveredLiquidCells)
    {
        recordedMassLB += entry.massLB;
        recordedMomentumLB += entry.momentumLB;
    }

    WALBERLA_LOG_INFO_ON_ROOT(
        "Covered-liquid records:"
        << " cells=" << shifted.coveredLiquidCells.size()
        << " massDifferenceLB="
        << recordedMassLB - shifted.coveredLiquidMassLB
        << " momentumDifferenceLB="
        << (recordedMomentumLB -
            shifted.coveredLiquidMomentumLB).norm()
        << " momentumLB="
        << recordedMomentumLB.transpose()
    );

    struct InterfaceRecipient
    {
        IBlock* block;
        Cell cell;
        Eigen::Vector3d positionLB;
    };

    std::vector<InterfaceRecipient> interfaceRecipients;

    for (auto& block : *blockForest)
    {
        const auto* flags =
            block.getData<FlagField_T>(flagFieldID);

        for (const auto& cell : flags->xyzSize())
        {
            const auto cellFlag =
                flags->get(cell.x(), cell.y(), cell.z());

            if (!flagInfo.isInterface(cellFlag))
                continue;

            const auto p =
                blockForest->getBlockLocalCellCenter(block, cell);

            const Eigen::Vector3d positionLB{
                static_cast<double>(p[0]),
                static_cast<double>(p[1]),
                static_cast<double>(p[2])
            };

            // A recipient must remain outside the proposed particle.
            if ((positionLB - shiftedCenter).squaredNorm()
                <= double(solidRadius) * double(solidRadius))
            {
                continue;
            }

            interfaceRecipients.push_back({
                &block, cell, positionLB
            });
        }
    }

            auto findInterfaceThroughLiquid =
                [&](const auto& source,
                    double minimumDonorMassLB = 0.0)
                    -> std::vector<std::size_t>
{
    std::vector<std::size_t> found;

    double reachedInterfaceMassLB = 0.0;

    const auto* flags =
        source.block->template getData<FlagField_T>(flagFieldID);

                    const auto* fill =
    source.block->template getData<ScalarField_T>(fillFieldID);

                    const auto* pdfs =
                        source.block->template getData<PdfField_T>(pdfFieldID);

    const auto ownedCells = flags->xyzSize();

    using CellKey = std::array<cell_idx_t, 3>;

    std::set<CellKey> visited;
    std::queue<Cell> frontier;

    visited.insert({
        source.cell.x(),
        source.cell.y(),
        source.cell.z()
    });

    // Start at the covered cell, then enter surviving fluid.
    frontier.push(source.cell);

    while (!frontier.empty())
    {
        const std::size_t layerSize = frontier.size();

        for (std::size_t n = 0; n < layerSize; ++n)
        {
            const Cell current = frontier.front();
            frontier.pop();

            for (std::size_t j = 0;
                 j < interfaceRecipients.size();
                 ++j)
            {
                const auto& recipient = interfaceRecipients[j];

                if (recipient.block == source.block &&
                    recipient.cell == current)
                {
                    found.push_back(j);

                    reachedInterfaceMassLB +=
                        double(fill->get(
                            current.x(), current.y(), current.z())) *
                        double(pdfs->getDensity(
                            current.x(), current.y(), current.z()));

                    break;
                }
            }

            for (auto d = LatticeModelStencil_T::beginNoCenter();
                 d != LatticeModelStencil_T::end();
                 ++d)
            {
                const Cell next(
                    current.x() + d.cx(),
                    current.y() + d.cy(),
                    current.z() + d.cz());

                if (!ownedCells.contains(next))
                    continue;

                const CellKey key{
                    next.x(), next.y(), next.z()
                };

                if (!visited.insert(key).second)
                    continue;

                const auto flag =
                    flags->get(next.x(), next.y(), next.z());

                if (!flagInfo.isLiquid(flag) &&
                    !flagInfo.isInterface(flag))
                {
                    continue;
                }

                const auto p =
                    blockForest->getBlockLocalCellCenter(
                        *source.block, next);

                const Eigen::Vector3d positionLB{
                    static_cast<double>(p[0]),
                    static_cast<double>(p[1]),
                    static_cast<double>(p[2])
                };

                // Never route through the proposed solid.
                if ((positionLB - shiftedCenter).squaredNorm()
                    <= double(solidRadius) * double(solidRadius))
                {
                    continue;
                }

                frontier.push(next);
            }
        }

        // Finish the whole layer to collect equally near recipients.
        if (!found.empty() &&
    reachedInterfaceMassLB >= minimumDonorMassLB)
        {
            return found;
        }
    }

    return found;
};

    struct CoveredCellRoute
    {
        std::size_t sourceIndex;
        std::vector<std::size_t> recipientIndices;
    };

    std::vector<CoveredCellRoute> coveredCellRoutes;

    double locallyRoutableMassLB = 0.0;
    double unroutedMassLB = 0.0;

    for (std::size_t i = 0;
         i < shifted.coveredLiquidCells.size();
         ++i)
    {
        const auto& source = shifted.coveredLiquidCells[i];

        const auto p =
            blockForest->getBlockLocalCellCenter(
                *source.block, source.cell);

        const Eigen::Vector3d sourcePositionLB{
            static_cast<double>(p[0]),
            static_cast<double>(p[1]),
            static_cast<double>(p[2])
        };

        CoveredCellRoute route;
        route.sourceIndex = i;

        for (std::size_t j = 0;
             j < interfaceRecipients.size();
             ++j)
        {
            const double distanceSquared =
                (interfaceRecipients[j].positionLB -
                 sourcePositionLB).squaredNorm();

            // D3Q19 neighbors: face or edge neighbors.
            if (distanceSquared > 0.5 &&
                distanceSquared <= 2.0 + 1e-12)
            {
                route.recipientIndices.push_back(j);
            }
        }

        if (route.recipientIndices.empty())
        {
            route.recipientIndices =
                findInterfaceThroughLiquid(source);
        }

        if (route.recipientIndices.empty())
            unroutedMassLB += source.massLB;
        else
            locallyRoutableMassLB += source.massLB;

        coveredCellRoutes.push_back(std::move(route));
    }

    WALBERLA_LOG_INFO_ON_ROOT(
        "Displacement routing:"
        << " survivingInterfaceCells=" << interfaceRecipients.size()
        << " routableMassLB=" << locallyRoutableMassLB
        << " unroutedMassLB=" << unroutedMassLB
        << " balanceErrorLB="
        << locallyRoutableMassLB + unroutedMassLB
           - shifted.coveredLiquidMassLB
    );

    struct LiquidTransfer
    {
        double massLB = 0.0;
        Eigen::Vector3d momentumLB =
            Eigen::Vector3d::Zero();
    };

    std::vector<LiquidTransfer> recipientTransfers(
        interfaceRecipients.size());

    LiquidTransfer pendingTransfer;

    for (const auto& route : coveredCellRoutes)
    {
        const auto& source =
            shifted.coveredLiquidCells.at(route.sourceIndex);

        if (route.recipientIndices.empty())
        {
            // Keep an explicit record of everything not yet routed.
            pendingTransfer.massLB += source.massLB;
            pendingTransfer.momentumLB += source.momentumLB;
            continue;
        }

        const double weight =
            1.0 / double(route.recipientIndices.size());

        for (const std::size_t index : route.recipientIndices)
        {
            auto& transfer = recipientTransfers.at(index);

            transfer.massLB += weight * source.massLB;
            transfer.momentumLB += weight * source.momentumLB;
        }
    }

    LiquidTransfer allocatedTransfer;

    for (const auto& transfer : recipientTransfers)
    {
        allocatedTransfer.massLB += transfer.massLB;
        allocatedTransfer.momentumLB += transfer.momentumLB;
    }

    const double massBalanceErrorLB =
        allocatedTransfer.massLB +
        pendingTransfer.massLB -
        shifted.coveredLiquidMassLB;

    const Eigen::Vector3d momentumBalanceErrorLB =
        allocatedTransfer.momentumLB +
        pendingTransfer.momentumLB -
        shifted.coveredLiquidMomentumLB;

    WALBERLA_LOG_INFO_ON_ROOT(
        "Displacement transfer plan:"
        << " allocatedMassLB=" << allocatedTransfer.massLB
        << " pendingMassLB=" << pendingTransfer.massLB
        << " massErrorLB=" << massBalanceErrorLB
        << " momentumErrorLB="
        << momentumBalanceErrorLB.norm()
    );

    struct UncoveredCellState
{
    IBlock* block;
    Cell cell;

    double fill;
    double density;
    Eigen::Vector3d velocityLB;
};

std::vector<UncoveredCellState> uncoveredStates;

double requestedRefillMassLB = 0.0;
uint_t unsupportedUncoveredCells = 0;

for (const auto& entry : shifted.uncoveredCells)
{
    const auto* flags =
        entry.block->getData<FlagField_T>(flagFieldID);

    const auto* fill =
        entry.block->getData<ScalarField_T>(fillFieldID);

    const auto* pdfs =
        entry.block->getData<PdfField_T>(pdfFieldID);

    double fractionSum = 0.0;
    double liquidMassSum = 0.0;

    Eigen::Vector3d liquidMomentumSum =
        Eigen::Vector3d::Zero();

    uint_t neighborCount = 0;

    for (auto d = LatticeModelStencil_T::beginNoCenter();
         d != LatticeModelStencil_T::end();
         ++d)
    {
        const Cell neighbor(
            entry.cell.x() + d.cx(),
            entry.cell.y() + d.cy(),
            entry.cell.z() + d.cz());

        // Single-block development implementation.
        if (!flags->xyzSize().contains(neighbor))
            continue;

        const auto flag = flags->get(
            neighbor.x(), neighbor.y(), neighbor.z());

        const bool isFluid =
            flagInfo.isLiquid(flag) || flagInfo.isInterface(flag);

        if (!isFluid && !flagInfo.isGas(flag))
            continue;

        const auto p =
            blockForest->getBlockLocalCellCenter(
                *entry.block, neighbor);

        const Eigen::Vector3d positionLB{
            static_cast<double>(p[0]),
            static_cast<double>(p[1]),
            static_cast<double>(p[2])
        };

        // Only use neighbors surviving the proposed move.
        if ((positionLB - shiftedCenter).squaredNorm()
            <= double(solidRadius) * double(solidRadius))
        {
            continue;
        }

        ++neighborCount;

        // Gas contributes zero liquid fraction.
        if (!isFluid)
            continue;

        const double fraction = double(fill->get(
            neighbor.x(), neighbor.y(), neighbor.z()));

        const double density = double(pdfs->getDensity(
            neighbor.x(), neighbor.y(), neighbor.z()));

        const auto momentum = pdfs->getMomentumDensity(
            neighbor.x(), neighbor.y(), neighbor.z());

        fractionSum += fraction;
        liquidMassSum += fraction * density;

        liquidMomentumSum += fraction * Eigen::Vector3d{
            static_cast<double>(momentum[0]),
            static_cast<double>(momentum[1]),
            static_cast<double>(momentum[2])
        };
    }

    if (neighborCount == 0)
    {
        ++unsupportedUncoveredCells;
        continue;
    }

    const double proposedFill = std::clamp(
        fractionSum / double(neighborCount), 0.0, 1.0);

    double proposedDensity = 1.0;
    Eigen::Vector3d proposedVelocity =
        Eigen::Vector3d::Zero();

    if (fractionSum > 1e-12 && liquidMassSum > 1e-12)
    {
        proposedDensity = liquidMassSum / fractionSum;
        proposedVelocity = liquidMomentumSum / liquidMassSum;
    }

    uncoveredStates.push_back({
        entry.block,
        entry.cell,
        proposedFill,
        proposedDensity,
        proposedVelocity
    });

    requestedRefillMassLB += proposedFill * proposedDensity;
}

WALBERLA_LOG_INFO_ON_ROOT(
    "Uncovered-cell refill estimate:"
    << " plannedCells=" << uncoveredStates.size()
    << " unsupported=" << unsupportedUncoveredCells
    << " requestedMassLB=" << requestedRefillMassLB
    << " coveredMinusRefillMassLB="
    << shifted.coveredLiquidMassLB - requestedRefillMassLB
);

    // Use interface cells reached by the displacement-routing search.
std::vector<bool> eligibleRecipient(
    interfaceRecipients.size(), false);

for (const auto& route : coveredCellRoutes)
{
    for (const std::size_t index : route.recipientIndices)
        eligibleRecipient.at(index) = true;
}

            // Newly uncovered wet cells also need access to liquid donors.
            for (const auto& state : uncoveredStates)
            {
                if (state.fill <= 0.0)
                    continue;

                const double requiredDonorMassLB = std::max(
                    0.0,
                    requestedRefillMassLB - shifted.coveredLiquidMassLB);

                const auto donorIndices =
                    findInterfaceThroughLiquid(state, requiredDonorMassLB);

                for (const std::size_t index : donorIndices)
                    eligibleRecipient.at(index) = true;
            }

std::vector<double> recipientMassDeltaLB(
    interfaceRecipients.size(), 0.0);

std::vector<double> availableMassLB(
    interfaceRecipients.size(), 0.0);

double totalAvailableMassLB = 0.0;
std::size_t eligibleCount = 0;

for (std::size_t j = 0; j < interfaceRecipients.size(); ++j)
{
    if (!eligibleRecipient[j])
        continue;

    ++eligibleCount;

    const auto& recipient = interfaceRecipients[j];

    const auto* fill =
        recipient.block->getData<ScalarField_T>(fillFieldID);

    const auto* pdfs =
        recipient.block->getData<PdfField_T>(pdfFieldID);

    const auto& cell = recipient.cell;

    const double fraction = double(
        fill->get(cell.x(), cell.y(), cell.z()));

    const double density = double(
        pdfs->getDensity(cell.x(), cell.y(), cell.z()));

    const double mass = fraction * density;

    if (!std::isfinite(mass) || mass < 0.0)
    {
        throw std::runtime_error(
            "Invalid interface mass in redistribution plan");
    }

    availableMassLB[j] = mass;
    totalAvailableMassLB += mass;
}

const double remainingMassLB =
    shifted.coveredLiquidMassLB - requestedRefillMassLB;

bool massPlanSupported =
    unsupportedUncoveredCells == 0 &&
    pendingTransfer.massLB == 0.0;

if (remainingMassLB < 0.0)
{
    const double withdrawalLB = -remainingMassLB;

    if (totalAvailableMassLB <= 0.0 ||
        withdrawalLB > totalAvailableMassLB)
    {
        massPlanSupported = false;
    }
    else
    {
        // Withdraw proportionally to each donor's liquid mass.
        // This cannot withdraw more than the donor contains.
        for (std::size_t j = 0;
             j < interfaceRecipients.size();
             ++j)
        {
            recipientMassDeltaLB[j] =
                -withdrawalLB *
                availableMassLB[j] / totalAvailableMassLB;
        }
    }
}
else if (remainingMassLB > 0.0)
{
    if (eligibleCount == 0)
    {
        massPlanSupported = false;
    }
    else
    {
        for (std::size_t j = 0;
             j < interfaceRecipients.size();
             ++j)
        {
            if (eligibleRecipient[j])
            {
                recipientMassDeltaLB[j] =
                    remainingMassLB / double(eligibleCount);
            }
        }
    }
}

double summedRecipientDeltaLB = 0.0;

for (const double delta : recipientMassDeltaLB)
    summedRecipientDeltaLB += delta;

WALBERLA_LOG_INFO_ON_ROOT(
    "Combined mass plan:"
    << " supported=" << massPlanSupported
    << " eligibleInterfaceCells=" << eligibleCount
    << " availableMassLB=" << totalAvailableMassLB
    << " recipientDeltaLB=" << summedRecipientDeltaLB
    << " balanceErrorLB="
    << requestedRefillMassLB + summedRecipientDeltaLB
       - shifted.coveredLiquidMassLB
);

    auto cellCenterEigen =
    [&](IBlock& block, const Cell& cell) -> Eigen::Vector3d
{
    const auto p =
        blockForest->getBlockLocalCellCenter(block, cell);

    return Eigen::Vector3d{
        static_cast<double>(p[0]),
        static_cast<double>(p[1]),
        static_cast<double>(p[2])
    };
};

Eigen::Vector3d fluidMomentumChangeLB =
    Eigen::Vector3d::Zero();

Eigen::Vector3d fluidAngularMomentumChangeLB =
    Eigen::Vector3d::Zero();

// Use one fixed reference point for every angular-momentum term.
const Eigen::Vector3d angularReferenceLB = solidCenterEigen;

auto accumulateMomentumChange =
    [&](const Eigen::Vector3d& positionLB,
        const Eigen::Vector3d& momentumChangeLB)
{
    fluidMomentumChangeLB += momentumChangeLB;

    fluidAngularMomentumChangeLB +=
        (positionLB - angularReferenceLB)
            .cross(momentumChangeLB);
};

// Covered cells lose their liquid and its momentum.
for (const auto& source : shifted.coveredLiquidCells)
{
    accumulateMomentumChange(
        cellCenterEigen(*source.block, source.cell),
        -source.momentumLB);
}

// Newly uncovered cells acquire their proposed liquid momentum.
for (const auto& state : uncoveredStates)
{
    const double massLB = state.fill * state.density;

    accumulateMomentumChange(
        cellCenterEigen(*state.block, state.cell),
        massLB * state.velocityLB);
}

// Surviving interface cells gain or lose liquid at their
// existing velocity. Their PDFs are not changed by this plan.
for (std::size_t j = 0; j < interfaceRecipients.size(); ++j)
{
    const double massDeltaLB = recipientMassDeltaLB[j];

    if (massDeltaLB == 0.0)
        continue;

    const auto& recipient = interfaceRecipients[j];

    const auto* pdfs =
        recipient.block->getData<PdfField_T>(pdfFieldID);

    const auto& cell = recipient.cell;

    const double density = double(
        pdfs->getDensity(cell.x(), cell.y(), cell.z()));

    if (!std::isfinite(density) || density <= 0.0)
    {
        throw std::runtime_error(
            "Invalid recipient density in momentum plan");
    }

    const auto momentumDensity = pdfs->getMomentumDensity(
        cell.x(), cell.y(), cell.z());

    const Eigen::Vector3d velocityLB{
        double(momentumDensity[0]) / density,
        double(momentumDensity[1]) / density,
        double(momentumDensity[2]) / density
    };

    accumulateMomentumChange(
        recipient.positionLB,
        massDeltaLB * velocityLB);
}

const Eigen::Vector3d particleRemapImpulseLB =
    -fluidMomentumChangeLB;

const Eigen::Vector3d particleRemapAngularImpulseLB =
    -fluidAngularMomentumChangeLB;

WALBERLA_LOG_INFO_ON_ROOT(
    "Remapping impulse plan:"
    << " particleImpulseLB="
    << particleRemapImpulseLB.transpose()
    << " angularImpulseAboutOldCenterLB="
    << particleRemapAngularImpulseLB.transpose()
);

    struct PlannedFluidCell
    {
        IBlock* block;
        Cell cell;

        double fill;
        double density;
        Eigen::Vector3d velocityLB;

        bool initializePdfs;
    };

    std::vector<PlannedFluidCell> plannedFluidCells;

    // Newly uncovered cells need initialized PDFs.
    for (const auto& state : uncoveredStates)
    {
        plannedFluidCells.push_back({
            state.block,
            state.cell,
            state.fill,
            state.density,
            state.velocityLB,
            true
        });
    }

    // Existing interface cells retain their PDFs.
    // Their liquid mass changes through their fill level.
    for (std::size_t j = 0; j < interfaceRecipients.size(); ++j)
    {
        const double massDeltaLB = recipientMassDeltaLB[j];

        if (massDeltaLB == 0.0)
            continue;

        const auto& recipient = interfaceRecipients[j];
        const auto& cell = recipient.cell;

        const auto* fill =
            recipient.block->getData<ScalarField_T>(fillFieldID);

        const auto* pdfs =
            recipient.block->getData<PdfField_T>(pdfFieldID);

        const double oldFill = double(
            fill->get(cell.x(), cell.y(), cell.z()));

        const double density = double(
            pdfs->getDensity(cell.x(), cell.y(), cell.z()));

        if (!std::isfinite(density) || density <= 0.0)
        {
            throw std::runtime_error(
                "Invalid density while planning fluid cells");
        }

        const auto velocity = pdfs->getVelocity(
            cell.x(), cell.y(), cell.z());

        const Eigen::Vector3d velocityLB{
            static_cast<double>(velocity[0]),
            static_cast<double>(velocity[1]),
            static_cast<double>(velocity[2])
        };

        const double proposedFill =
            oldFill + massDeltaLB / density;

        plannedFluidCells.push_back({
            recipient.block,
            cell,
            proposedFill,
            density,
            velocityLB,
            false
        });
    }

    uint_t invalidPlannedCells = 0;
    uint_t pdfInitializations = 0;

    double minimumPlannedFill = 1.0;
    double maximumPlannedFill = 0.0;

    for (const auto& state : plannedFluidCells)
    {
        if (state.initializePdfs)
            ++pdfInitializations;

        const bool valid =
            std::isfinite(state.fill) &&
            std::isfinite(state.density) &&
            state.density > 0.0 &&
            state.velocityLB.allFinite() &&
            state.fill >= 0.0 &&
            state.fill <= 1.0;

        if (!valid)
        {
            ++invalidPlannedCells;
            continue;
        }

        minimumPlannedFill =
            std::min(minimumPlannedFill, state.fill);

        maximumPlannedFill =
            std::max(maximumPlannedFill, state.fill);
    }

    massPlanSupported =
        massPlanSupported && invalidPlannedCells == 0;

    WALBERLA_LOG_INFO_ON_ROOT(
        "Proposed fluid states:"
        << " cells=" << plannedFluidCells.size()
        << " pdfInitializations=" << pdfInitializations
        << " invalid=" << invalidPlannedCells
        << " minFill=" << minimumPlannedFill
        << " maxFill=" << maximumPlannedFill
        << " supported=" << massPlanSupported
    );

    uint_t promotedLiquidCells = 0;

for (const auto& uncovered : uncoveredStates)
{
    // Positive-fill uncovered cells will become interface cells.
    // Exactly zero-fill uncovered cells will become gas.
    if (uncovered.fill != 0.0)
        continue;

    auto* block = uncovered.block;

    const auto* flags =
        block->getData<FlagField_T>(flagFieldID);

    const auto* fill =
        block->getData<ScalarField_T>(fillFieldID);

    const auto* pdfs =
        block->getData<PdfField_T>(pdfFieldID);

    for (auto d = LatticeModelStencil_T::beginNoCenter();
         d != LatticeModelStencil_T::end();
         ++d)
    {
        const Cell neighbor(
            uncovered.cell.x() + d.cx(),
            uncovered.cell.y() + d.cy(),
            uncovered.cell.z() + d.cz());

        if (!flags->xyzSize().contains(neighbor))
        {
            throw std::runtime_error(
                "Topology planning currently requires "
                "changed cells away from block boundaries");
        }

        const auto flag = flags->get(
            neighbor.x(), neighbor.y(), neighbor.z());

        if (!flagInfo.isLiquid(flag))
            continue;

        const Eigen::Vector3d positionLB =
            cellCenterEigen(*block, neighbor);

        // A cell covered by the new sphere becomes solid instead.
        if ((positionLB - shiftedCenter).squaredNorm()
            <= double(solidRadius) * double(solidRadius))
        {
            continue;
        }

        // Avoid duplicates. Already-planned fluid cells use
        // their planned phase rather than their original flag.
        const bool alreadyPlanned = std::any_of(
            plannedFluidCells.begin(),
            plannedFluidCells.end(),
            [&](const PlannedFluidCell& state)
            {
                return state.block == block &&
                       state.cell == neighbor;
            });

        if (alreadyPlanned)
            continue;

        const double oldFill = double(fill->get(
            neighbor.x(), neighbor.y(), neighbor.z()));

        const double density = double(pdfs->getDensity(
            neighbor.x(), neighbor.y(), neighbor.z()));

        const auto velocity = pdfs->getVelocity(
            neighbor.x(), neighbor.y(), neighbor.z());

        const Eigen::Vector3d velocityLB{
            static_cast<double>(velocity[0]),
            static_cast<double>(velocity[1]),
            static_cast<double>(velocity[2])
        };

        plannedFluidCells.push_back({
            block,
            neighbor,
            oldFill,
            density,
            velocityLB,
            false
        });

        ++promotedLiquidCells;
    }
}

WALBERLA_LOG_INFO_ON_ROOT(
    "Interface topology plan:"
    << " liquidToInterface=" << promotedLiquidCells
    << " totalPlannedFluidCells=" << plannedFluidCells.size()
);

    uint_t remapBlockCount = 0;

    for (auto it = blockForest->begin();
         it != blockForest->end(); ++it)
    {
        ++remapBlockCount;
    }

    if (numProcesses != uint_c(1) ||
        remapBlockCount != uint_c(1))
    {
        throw std::runtime_error(
            "The startup remap test currently requires "
            "one process and one block");
    }

    if (!massPlanSupported)
    {
        throw std::runtime_error(
            "Cannot apply an unsupported remapping plan");
    }



            struct StoredExcessTransfer
         {
             IBlock* block;
             Cell source;
             Cell destination;
             double massLB;
         };

            std::vector<StoredExcessTransfer> storedExcessTransfers;

            // Validate all states, including the later topology promotions.
            for (const auto& state : plannedFluidCells)
            {
                if (!std::isfinite(state.fill) ||
                    state.fill < 0.0 || state.fill > 1.0 ||
                    !std::isfinite(state.density) ||
                    state.density <= 0.0 ||
                    !state.velocityLB.allFinite())
                {
                    throw std::runtime_error(
                        "Invalid fluid state before remapping");
                }

                if (!state.initializePdfs &&
    remapExcessMassFieldID != ConstBlockDataID())
{
    const auto* storedExcess =
        state.block->getData<ScalarField_T>(
            remapExcessMassFieldID);

    const double extra = double(storedExcess->get(
        state.cell.x(), state.cell.y(), state.cell.z()));

    if (!std::isfinite(extra))
    {
        throw std::runtime_error(
            "Non-finite stored excess mass");
    }

    if (extra != 0.0)
    {
        const auto* flags =
            state.block->getData<FlagField_T>(flagFieldID);

        bool destinationFound = false;

        for (auto d = LatticeModelStencil_T::beginNoCenter();
             d != LatticeModelStencil_T::end();
             ++d)
        {
            const Cell neighbor(
                state.cell.x() + d.cx(),
                state.cell.y() + d.cy(),
                state.cell.z() + d.cz());

            if (!flags->xyzSize().contains(neighbor))
                continue;

            const auto flag = flags->get(
                neighbor.x(), neighbor.y(), neighbor.z());

            if (!flagInfo.isLiquid(flag))
                continue;

            const Eigen::Vector3d positionLB =
                cellCenterEigen(*state.block, neighbor);

            // The destination must remain outside the particle.
            if ((positionLB - shiftedCenter).squaredNorm()
                <= double(solidRadius) * double(solidRadius))
            {
                continue;
            }

            // It must also remain an unchanged liquid cell.
            const bool destinationChanges = std::any_of(
                plannedFluidCells.begin(),
                plannedFluidCells.end(),
                [&](const PlannedFluidCell& candidate)
                {
                    return candidate.block == state.block &&
                           candidate.cell == neighbor;
                });

            if (destinationChanges)
                continue;

            storedExcessTransfers.push_back({
                state.block,
                state.cell,
                neighbor,
                extra
            });

            destinationFound = true;
            break;
        }

        if (!destinationFound)
        {
            throw std::runtime_error(
                "No unchanged neighboring liquid cell can "
                "receive stored excess mass");
        }
    }
}
            }

    auto measureRemapLiquidState = [&]()
    {
        double massLB = 0.0;
        Eigen::Vector3d momentumLB = Eigen::Vector3d::Zero();

        for (auto& block : *blockForest)
        {
            const auto* flags =
                block.getData<FlagField_T>(flagFieldID);

            const auto* fill =
                block.getData<ScalarField_T>(fillFieldID);

            const auto* pdfs =
                block.getData<PdfField_T>(pdfFieldID);

            const ScalarField_T* storedExcess = nullptr;

            if (remapExcessMassFieldID != ConstBlockDataID())
            {
                storedExcess = block.getData<ScalarField_T>(
                    remapExcessMassFieldID);
            }

            for (const auto& cell : flags->xyzSize())
            {
                const auto flag =
                    flags->get(cell.x(), cell.y(), cell.z());

                if (!flagInfo.isLiquid(flag) &&
                    !flagInfo.isInterface(flag))
                {
                    continue;
                }

                const double fraction = double(
                    fill->get(cell.x(), cell.y(), cell.z()));

                massLB += fraction * double(pdfs->getDensity(
                    cell.x(), cell.y(), cell.z()));

                if (storedExcess != nullptr)
                {
                    massLB += double(storedExcess->get(
                        cell.x(), cell.y(), cell.z()));
                }

                const auto p = pdfs->getMomentumDensity(
                    cell.x(), cell.y(), cell.z());

                momentumLB += fraction * Eigen::Vector3d{
                    static_cast<double>(p[0]),
                    static_cast<double>(p[1]),
                    static_cast<double>(p[2])
                };
            }
        }

        return std::make_pair(massLB, momentumLB);
    };

    const auto beforeRemap = measureRemapLiquidState();

    // Newly covered cells become particle boundary cells.
for (const auto& entry : shifted.coveredCells)
{
    auto* handling = entry.block->getData<BoundaryHandling_T>(
        freeSurfaceBoundaryHandling->getHandlingID());

    auto* flags = handling->getFlagField();

    auto* fill =
        entry.block->getData<ScalarField_T>(fillFieldID);

    const auto& cell = entry.cell;

    if (remapExcessMassFieldID != ConstBlockDataID())
    {
        const BlockDataID writableExcessID{
            static_cast<uint_t>(remapExcessMassFieldID)
        };

        auto* storedExcess =
            entry.block->getData<ScalarField_T>(writableExcessID);

        storedExcess->get(
            cell.x(), cell.y(), cell.z()) = real_c(0);
    }

    // Gas is not part of BoundaryHandling's fluid-domain mask,
    // so remove its flag explicitly.
    if (flagInfo.isGas(
            flags->get(cell.x(), cell.y(), cell.z())))
    {
        handling->removeFlag(
            flagInfo.gasFlag,
            cell.x(), cell.y(), cell.z());
    }

    handling->forceBoundary(
        FreeSurfaceBoundaryHandling_T::ubbFlagID,
        cell.x(), cell.y(), cell.z(),
        UBB_T::Velocity(real_c(0), real_c(0), real_c(0)));

    // Match waLBerla's obstacle-fill initialization.
    // Solid cells are excluded from liquid-mass measurements.
    fill->get(cell.x(), cell.y(), cell.z()) = real_c(1);
}

// Apply uncovered-cell states, recipient changes,
// and liquid-to-interface promotions.
for (const auto& state : plannedFluidCells)
{
    auto* handling = state.block->getData<BoundaryHandling_T>(
        freeSurfaceBoundaryHandling->getHandlingID());

    auto* fill =
        state.block->getData<ScalarField_T>(fillFieldID);

    auto* pdfs =
        state.block->getData<PdfField_T>(pdfFieldID);

    const auto& cell = state.cell;

    handling->clear(cell.x(), cell.y(), cell.z());

    const bool becomesGas =
        state.initializePdfs && state.fill == 0.0;

    if (becomesGas)
    {
        handling->setFlag(
            flagInfo.gasFlag,
            cell.x(), cell.y(), cell.z());
    }
    else
    {
        handling->setDomain(
            flagInfo.interfaceFlag,
            cell.x(), cell.y(), cell.z());
    }

    fill->get(cell.x(), cell.y(), cell.z()) =
        real_c(state.fill);

    if (state.initializePdfs)
    {

        if (remapExcessMassFieldID != ConstBlockDataID())
        {
            const BlockDataID writableExcessID{
                static_cast<uint_t>(remapExcessMassFieldID)
            };

            auto* storedExcess =
                state.block->getData<ScalarField_T>(writableExcessID);

            storedExcess->get(
                cell.x(), cell.y(), cell.z()) = real_c(0);
        }

        pdfs->setToEquilibrium(
            cell,
            Vector3<real_t>(
                real_c(state.velocityLB.x()),
                real_c(state.velocityLB.y()),
                real_c(state.velocityLB.z())),
            real_c(state.density));
    }
}

for (auto& block : *blockForest)
{
    auto* handling = block.getData<BoundaryHandling_T>(
        freeSurfaceBoundaryHandling->getHandlingID());

    handling->refresh();
}

            if (!storedExcessTransfers.empty())
            {
                const BlockDataID writableExcessID{
                    static_cast<uint_t>(remapExcessMassFieldID)
                };

                for (const auto& transfer : storedExcessTransfers)
                {
                    auto* storedExcess =
                        transfer.block->getData<ScalarField_T>(
                            writableExcessID);

                    storedExcess->get(
                        transfer.destination.x(),
                        transfer.destination.y(),
                        transfer.destination.z()) += real_c(transfer.massLB);

                    storedExcess->get(
                        transfer.source.x(),
                        transfer.source.y(),
                        transfer.source.z()) = real_c(0);
                }
            }

            WALBERLA_LOG_INFO_ON_ROOT(
                "Stored excess transfers=" << storedExcessTransfers.size()
            );

    const auto afterRemap = measureRemapLiquidState();

    const Eigen::Vector3d momentumPredictionErrorLB =
        afterRemap.second -
        beforeRemap.second -
        fluidMomentumChangeLB;

    const auto geometryAfterRemap =
        detectParticleCellChanges(shiftedCenter);

    WALBERLA_LOG_INFO_ON_ROOT(
        "Applied particle remap:"
        << " massBeforeLB=" << beforeRemap.first
        << " massAfterLB=" << afterRemap.first
        << " massErrorLB="
        << afterRemap.first - beforeRemap.first
        << " momentumPredictionErrorLB="
        << momentumPredictionErrorLB.norm()
        << " geometryCovered=" << geometryAfterRemap.covered
        << " geometryUncovered=" << geometryAfterRemap.uncovered
        << " sphereCells=" << geometryAfterRemap.retained
    );

    // Synchronize all geometric representations with the remapped sphere.
    solidCenterEigen = shiftedCenter;

    solidCenter = Vector3<real_t>(
        real_c(solidCenterEigen.x()),
        real_c(solidCenterEigen.y()),
        real_c(solidCenterEigen.z()));

    solidSurface = makeSphereSurface(
        solidCenterEigen,
        double(solidRadius),
        4);

            return RemapResult{
                particleRemapImpulseLB,
                particleRemapAngularImpulseLB,
                angularReferenceLB
            };
        }; // End remapParticle.

   // communication after initialization
   Communication_T communication(blockForest, flagFieldID, fillFieldID, forceDensityFieldID);
   communication();

   PdfCommunication_T pdfCommunication(blockForest, pdfFieldID);
   pdfCommunication();

   // add bubble model
   std::shared_ptr< bubble_model::BubbleModelBase > bubbleModel = nullptr;
   if (enableBubbleModel)
   {
      const std::shared_ptr< bubble_model::BubbleModel< LatticeModelStencil_T > > bubbleModelDerived =
         std::make_shared< bubble_model::BubbleModel< LatticeModelStencil_T > >(blockForest, enableBubbleSplits);
      bubbleModelDerived->initFromFillLevelField(fillFieldID);
      bubbleModelDerived->setAtmosphere(
         Cell(domainSize[0] - uint_c(1), domainSize[1] - uint_c(1), domainSize[2] - uint_c(1)), real_c(1));

      bubbleModel = std::static_pointer_cast< bubble_model::BubbleModelBase >(bubbleModelDerived);
   }
   else { bubbleModel = std::make_shared< bubble_model::BubbleModelConstantPressure >(real_c(1)); }

   // initialize hydrostatic pressure
   // initHydrostaticPressure< PdfField_T >(blockForest, pdfFieldID, force, real_c(0.5) * dropDiameter);

   // set density in non-liquid or non-interface cells to 1 (after initializing with hydrostatic pressure)
   setDensityInNonFluidCellsToOne< FlagField_T, PdfField_T >(blockForest, flagInfo, flagFieldID, pdfFieldID);

   // create timeloop
   SweepTimeloop timeloop(blockForest, timesteps);

    const uint_t diagnosticFrequency = uint_c(50);

   Vector3<real_t> sphereForceLB(
       real_c(0), real_c(0), real_c(0)
   );

    Eigen::Vector3d sphereTorqueLB =
        Eigen::Vector3d::Zero();

   uint_t sphereLinkCount = uint_c(0);
    real_t contactLineZSum = real_c(0);
    uint_t contactLineCellCount = uint_c(0);

    auto updateSphereWallVelocity =
    [&](const Eigen::Vector3d& velocityLB,
        const Eigen::Vector3d& angularVelocityLB)
    {
        for (auto& block : *blockForest)
        {
            auto* handling = block.getData<BoundaryHandling_T>(
                freeSurfaceBoundaryHandling->getHandlingID());

            auto* flags = handling->getFlagField();

            const auto particleMask = flags->getFlag(
                FreeSurfaceBoundaryHandling_T::ubbFlagID);

            auto& ubb = handling->template getBoundaryCondition<UBB_T>(
                FreeSurfaceBoundaryHandling_T::ubbBoundaryID);

            for (const auto& cell : flags->xyzSizeWithGhostLayer())
            {
                if ((flags->get(cell.x(), cell.y(), cell.z())
                     & particleMask) == 0)
                    continue;

                const auto position =
                    blockForest->getBlockLocalCellCenter(block, cell);

                const Eigen::Vector3d displacement(
                    position[0] - solidCenter[0],
                    position[1] - solidCenter[1],
                    position[2] - solidCenter[2]);

                const Eigen::Vector3d wallVelocity =
                    velocityLB + angularVelocityLB.cross(displacement);

                const UBB_T::Velocity configuration(
                    real_c(wallVelocity.x()),
                    real_c(wallVelocity.y()),
                    real_c(wallVelocity.z()));

                ubb.registerCell(
                    particleMask,
                    cell.x(), cell.y(), cell.z(),
                    configuration);
            }
        }
    };

    double movingWallMassCorrectionLB = 0.0;

    double liquidWallFzLB = 0.0;
    double interfaceWallFzLB = 0.0;
    double weightedInterfaceWallFzLB = 0.0;

   auto measureSphereForce = [&]()
    {
        sphereForceLB = Vector3<real_t>(
            real_c(0), real_c(0), real_c(0)
        );
       sphereTorqueLB.setZero();
       liquidWallFzLB = 0.0;
       interfaceWallFzLB = 0.0;
       weightedInterfaceWallFzLB = 0.0;
        sphereLinkCount = uint_c(0);
       contactLineZSum = real_c(0);
       contactLineCellCount = uint_c(0);
       movingWallMassCorrectionLB = 0.0;

        for (auto& block : *blockForest)
        {
            const auto* pdfs =
                block.getData<PdfField_T>(pdfFieldID);

            const auto* flags =
                block.getData<FlagField_T>(flagFieldID);

            if (remapExcessMassFieldID == ConstBlockDataID())
            {
                throw std::runtime_error(
                    "Moving-wall mass correction requires "
                    "the fallback-liquid excess-mass field");
            }

            const BlockDataID writableExcessID{
                static_cast<uint_t>(remapExcessMassFieldID)
            };

            auto* wallExcessMass =
                block.getData<ScalarField_T>(writableExcessID);

            auto* handling = block.getData<BoundaryHandling_T>(
            freeSurfaceBoundaryHandling->getHandlingID());

            const auto& ubb =
                handling->template getBoundaryCondition<UBB_T>(
                    FreeSurfaceBoundaryHandling_T::ubbBoundaryID);

            const auto particleBoundaryMask = flags->getFlag(
                FreeSurfaceBoundaryHandling_T::ubbFlagID
            );

            // Visit owned fluid cells only, avoiding duplicate contributions
            // from ghost cells.
            for (const auto& cell : flags->xyzSize())
            {
                const auto fluidFlag = flags->get(
                    cell.x(), cell.y(), cell.z()
                );

                if (!flagInfo.isLiquid(fluidFlag) &&
                    !flagInfo.isInterface(fluidFlag))
                {
                    continue;
                }

                bool touchesSphere = false;

                for (auto d = LatticeModelStencil_T::beginNoCenter();
                     d != LatticeModelStencil_T::end();
                     ++d)
                {
                    const Cell neighbor(
                        cell.x() + d.cx(),
                        cell.y() + d.cy(),
                        cell.z() + d.cz()
                    );

                    const auto neighborFlag = flags->get(
                        neighbor.x(),
                        neighbor.y(),
                        neighbor.z()
                    );

                    if ((neighborFlag & particleBoundaryMask) == 0)
                        continue;

                    // The floor also uses NoSlip. Count only solid cells
                    // belonging to our sphere.
                    const auto neighborPosition =
                        blockForest->getBlockLocalCellCenter(
                            block, neighbor
                        );

                    if ((neighborPosition - solidCenter).length()
                        > solidRadius)
                    {
                        continue;
                    }

                    touchesSphere = true;

                    const real_t outgoingPDF = pdfs->get(
    cell.x(), cell.y(), cell.z(), d.toIdx());

                    const real_t density = pdfs->getDensity(
                        cell.x(), cell.y(), cell.z());

                    // Read the exact velocity used by the boundary condition.
                    const auto& wallVelocity = ubb.getValue(
                        neighbor.x(), neighbor.y(), neighbor.z());

                    const real_t directionDotVelocity =
                        real_c(d.cx()) * wallVelocity[0] +
                        real_c(d.cy()) * wallVelocity[1] +
                        real_c(d.cz()) * wallVelocity[2];

                    const real_t weight =
                        real_c(LatticeModel_T::w[d.toIdx()]);

                    // Same reflected PDF as waLBerla's compressible UBB.
                    const real_t reflectedPDF =
                        outgoingPDF -
                        real_c(6) * density * weight * directionDotVelocity;

                    // Full-liquid cells acquire the moving-wall PDF mass change
                    // during streaming. Compensate it through excess-mass storage.
                    if (flagInfo.isLiquid(fluidFlag))
                    {
                        const real_t massCorrection =
                            outgoingPDF - reflectedPDF;

                        wallExcessMass->get(
                            cell.x(), cell.y(), cell.z()) += massCorrection;

                        movingWallMassCorrectionLB += double(massCorrection);
                    }

                    // Subtract the ambient reference pressure.
                    const real_t rhoAmbient = real_c(1);
                    const real_t momentumTransfer =
                        outgoingPDF + reflectedPDF -
                        real_c(2) * weight * rhoAmbient;

                    const Eigen::Vector3d linkForce(
                        double(momentumTransfer) * double(d.cx()),
                        double(momentumTransfer) * double(d.cy()),
                        double(momentumTransfer) * double(d.cz()));

                    sphereForceLB[0] += real_c(linkForce.x());
                    sphereForceLB[1] += real_c(linkForce.y());
                    sphereForceLB[2] += real_c(linkForce.z());

                    if (flagInfo.isLiquid(fluidFlag)) {
                        liquidWallFzLB += linkForce.z();
                    } else {
                        interfaceWallFzLB += linkForce.z();

                        const auto* fill =
                            block.getData<ScalarField_T>(fillFieldID);

                        const double fraction = double(fill->get(
                            cell.x(), cell.y(), cell.z()));

                        weightedInterfaceWallFzLB += fraction * linkForce.z();
                    }

                    const auto fluidPosition =
                        blockForest->getBlockLocalCellCenter(block, cell);

                    const Eigen::Vector3d boundaryPosition(
                        0.5 * double(fluidPosition[0] + neighborPosition[0]),
                        0.5 * double(fluidPosition[1] + neighborPosition[1]),
                        0.5 * double(fluidPosition[2] + neighborPosition[2]));

                    sphereTorqueLB +=
                        (boundaryPosition - solidCenterEigen).cross(linkForce);

                    ++sphereLinkCount;
                }

                if (touchesSphere && flagInfo.isInterface(fluidFlag))
                {
                    const auto position =
                        blockForest->getBlockLocalCellCenter(block, cell);

                    const auto displacement = position - solidCenter;
                    const real_t distance = displacement.length();

                    if (distance > real_c(0))
                    {
                        // Project the interface-cell center onto the analytic sphere.
                        const real_t projectedZ =
                            solidCenter[2]
                            + solidRadius * displacement[2] / distance;

                        contactLineZSum += projectedZ;
                        ++contactLineCellCount;
                    }
                }
            }
        }
    };

    timeloop.addFuncBeforeTimeStep(
    [&]()
        {
            const Eigen::Vector3d velocitySI =
                dem.particleVelocitySI();

            const Eigen::Vector3d angularVelocitySI =
                dem.particleAngularVelocitySI();

            if (!velocitySI.allFinite() ||
                !angularVelocitySI.allFinite())
            {
                throw std::runtime_error(
                    "DEM supplied non-finite particle velocity");
            }

            const Eigen::Vector3d angularVelocityLB =
                angularVelocitySI * units.dtSI;

            const Eigen::Vector3d velocityLB =
                velocitySI * (units.dtSI / units.dxSI);

            updateSphereWallVelocity(
                velocityLB,
                angularVelocityLB);

            measureSphereForce();
        },
        "Apply DEM wall velocity and measure fluid load"
    );

   // Laplace pressure = 2 * surface tension * curvature; curvature computation is not necessary with 0 surface tension
   bool computeCurvature = false;
   if (!realIsEqual(surfaceTension, real_c(0), real_c(1e-14))) { computeCurvature = true; }

   // add surface geometry handler
   const SurfaceGeometryHandler< LatticeModel_T, FlagField_T, ScalarField_T, VectorField_T > geometryHandler(
      blockForest, freeSurfaceBoundaryHandling, fillFieldID, curvatureModel, computeCurvature, enableWetting,
      contactAngle);

    auto measureLiquidMomentumLB = [&]() -> Eigen::Vector3d
    {
        Eigen::Vector3d momentum = Eigen::Vector3d::Zero();

        for (auto& block : *blockForest)
        {
            const auto* flags =
                block.getData<FlagField_T>(flagFieldID);

            const auto* fill =
                block.getData<ScalarField_T>(fillFieldID);

            const auto* pdfs =
                block.getData<PdfField_T>(pdfFieldID);

            for (const auto& cell : flags->xyzSize())
            {
                const auto flag =
                    flags->get(cell.x(), cell.y(), cell.z());

                if (!flagInfo.isLiquid(flag) &&
                    !flagInfo.isInterface(flag))
                    continue;

                const double fraction = double(
                    fill->get(cell.x(), cell.y(), cell.z()));

                const auto p = pdfs->getMomentumDensity(
                    cell.x(), cell.y(), cell.z());

                momentum += fraction * Eigen::Vector3d{
                    static_cast<double>(p[0]),
                    static_cast<double>(p[1]),
                    static_cast<double>(p[2])
                };
            }
        }

        return momentum;
    };

    Eigen::Vector3d momentumAfterGeometryLB =
    Eigen::Vector3d::Zero();

   geometryHandler.addSweeps(timeloop);

    timeloop.add() << Sweep(
    [&](IBlock*)
    {
        momentumAfterGeometryLB =
            measureLiquidMomentumLB();
    },
    "Audit momentum after geometry"
);

   const ConstBlockDataID curvatureFieldID = geometryHandler.getConstCurvatureFieldID();
   const ConstBlockDataID normalFieldID    = geometryHandler.getConstNormalFieldID();

    timeloop.addFuncBeforeTimeStep(
    [&]() {
        interfacePlanes.clear();
    },
    "Reset interface-plane samples"
);

timeloop.add() << Sweep(
    [&](IBlock* block)
    {
        using Vec = Eigen::Vector3d;

        const auto* flags =
            block->getData<FlagField_T>(flagFieldID);

        const auto* fill =
            block->getData<ScalarField_T>(fillFieldID);

        const auto* normals =
            block->getData<VectorField_T>(normalFieldID);

        for (const auto& cell : flags->xyzSize())
        {
            const auto flag = flags->get(
                cell.x(), cell.y(), cell.z()
            );

            if (!flagInfo.isInterface(flag))
                continue;

            const auto cellCenter =
                blockForest->getBlockLocalCellCenter(
                    *block, cell
                );

            const Vec center(
                cellCenter[0],
                cellCenter[1],
                cellCenter[2]
            );

            // Collect planes within three lattice spacings
            // of the analytic sphere surface.
            const double surfaceDistance = std::abs(
                (center - solidCenterEigen).norm()
                - double(solidRadius)
            );

            if (surfaceDistance > 3.0)
                continue;

            const real_t fraction = fill->get(
                cell.x(), cell.y(), cell.z()
            );

            if (!std::isfinite(fraction) ||
                fraction <= real_c(1e-6) ||
                fraction >= real_c(1) - real_c(1e-6))
            {
                continue;
            }

            auto normal = normals->get(
                cell.x(), cell.y(), cell.z()
            );

            const Vec normalEigen(
                normal[0], normal[1], normal[2]
            );

            if (!normalEigen.allFinite() ||
                normalEigen.squaredNorm() < 1e-14)
            {
                continue;
            }

            normal = normal.getNormalizedOrZero();

            const auto localPoint =
                walberla::free_surface::getInterfacePoint(
                    normal, fraction
                );

            // getInterfacePoint uses [0,1]^3 coordinates.
            // This demo has unit lattice spacing.
            const Vec globalPoint = center + Vec(
                localPoint[0] - real_c(0.5),
                localPoint[1] - real_c(0.5),
                localPoint[2] - real_c(0.5)
            );

            if (!globalPoint.allFinite())
                continue;

            interfacePlanes.push_back({
                globalPoint,
                normalEigen.normalized()
            });
        }
    },
    "Collect interface planes near sphere"
);

#if 0
    timeloop.add() << Sweep(
            [&](IBlock* block)
            {
                using EVec = Eigen::Vector3d;

                const auto* flags =
                    block->getData<FlagField_T>(flagFieldID);

                const auto* fill =
                    block->getData<ScalarField_T>(fillFieldID);

                const auto* normals =
                    block->getData<VectorField_T>(normalFieldID);

                const auto noSlipMask = flags->getFlag(
                    FreeSurfaceBoundaryHandling_T::noSlipFlagID
                );

                const EVec particleCenter(
                    solidCenter[0],
                    solidCenter[1],
                    solidCenter[2]
                );

                for (const auto& cell : flags->xyzSize())
                {
                    const auto flag = flags->get(
                        cell.x(), cell.y(), cell.z()
                    );

                    if (!flagInfo.isInterface(flag))
                        continue;

                    const real_t fraction = fill->get(
                        cell.x(), cell.y(), cell.z()
                    );

                    // Do not reconstruct a plane from empty/full cells.
                    if (fraction <= real_c(1e-6) ||
                        fraction >= real_c(1) - real_c(1e-6))
                    {
                        continue;
                    }

                    auto normal = normals->get(
                        cell.x(), cell.y(), cell.z()
                    );

                    if (!std::isfinite(normal[0]) ||
                        !std::isfinite(normal[1]) ||
                        !std::isfinite(normal[2]) ||
                        normal.sqrLength() < real_c(1e-14))
                    {
                        ++invalidNormals;
                        continue;
                    }

                    normal = normal.getNormalizedOrZero();

                    const auto localPoint =
                        walberla::free_surface::getInterfacePoint(
                            normal, fraction
                        );

                    const EVec ni(
                        normal[0], normal[1], normal[2]
                    );

                    const EVec planePoint(
                        localPoint[0], localPoint[1], localPoint[2]
                    );

                    const auto cellCenter =
                        blockForest->getBlockLocalCellCenter(
                            *block, cell
                        );

                    // Current demo uses unit lattice spacing.
                    const EVec cellOrigin(
                        cellCenter[0] - real_c(0.5),
                        cellCenter[1] - real_c(0.5),
                        cellCenter[2] - real_c(0.5)
                    );

                    // Only the six face-sharing neighbors define cube faces.
                    for (int axis = 0; axis < 3; ++axis)
                    {
                        for (int side : {-1, 1})
                        {
                            Cell neighbor = cell;
                            neighbor[axis] += side;

                            const auto neighborFlag = flags->get(
                                neighbor.x(),
                                neighbor.y(),
                                neighbor.z()
                            );

                            if ((neighborFlag & noSlipMask) == 0)
                                continue;

                            const auto neighborPosition =
                                blockForest->getBlockLocalCellCenter(
                                    *block, neighbor
                                );

                            // Exclude floor and other NoSlip boundaries.
                            if ((neighborPosition - solidCenter).length()
                                > solidRadius)
                            {
                                continue;
                            }

                            EVec a, b;

                            const int result = interfacePlaneFaceSegment(
                                ni, planePoint, axis, side, a, b
                            );

                            if (result < 0) {
                                ++degenerateFaces;
                                continue;
                            }

                            if (result == 0)
                                continue;

                            // Normal from the solid voxel toward this fluid cell.
                            EVec ns = EVec::Zero();
                            ns[axis] = -double(side);

                            const EVec globalA = cellOrigin + a;
                            const EVec globalB = cellOrigin + b;

                            const double length = (globalB - globalA).norm();

                            const EVec midpoint = 0.5 * (globalA + globalB);

                            const EVec sphereNormal =
                                (midpoint - particleCenter).normalized();

                            const double normalAngle = std::acos(
                                std::clamp(ns.dot(sphereNormal), -1.0, 1.0)
                            );

                            reconstructedLength += length;
                            normalMismatchIntegral += length * normalAngle;

                            const auto load = contactLineSegmentLoad(
                                cellOrigin + a,
                                cellOrigin + b,
                                ns,
                                ni,
                                particleCenter,
                                double(surfaceTension)
                            );

                            reconstructedLineLoad.force += load.force;
                            reconstructedLineLoad.torque += load.torque;

                            ++reconstructedSegments;
                        }
                    }
                }
            },
            "Reconstruct contact-line segments"
        );
#endif

   // add boundary handling for standard boundaries and free surface boundaries
   const AuditedSurfaceDynamicsHandler< LatticeModel_T, FlagField_T, ScalarField_T, VectorField_T > dynamicsHandler(
      blockForest, pdfFieldID, flagFieldID, fillFieldID, forceDensityFieldID, normalFieldID, curvatureFieldID,
      freeSurfaceBoundaryHandling, bubbleModel, pdfReconstructionModel, pdfRefillingModel, excessMassDistributionModel,
      relaxationRate, force, surfaceTension, useSimpleMassExchange, cellConversionThreshold,
      cellConversionForceThreshold);

    remapExcessMassFieldID =
    dynamicsHandler.getConstExcessMassFieldID();

    // Synchronize fields after the imposed startup displacement.
    communication();
    pdfCommunication();

    uint_t auditedStep = uint_c(0);

    Eigen::Vector3d previousStageMomentumLB =
        Eigen::Vector3d::Zero();

    struct MomentumCellSample
    {
        IBlock* block;
        Cell cell;
        double fill;
    };

    struct MomentumLinkSample
    {
        IBlock* block;
        Cell cell;
        uint_t incomingIndex;
        double outgoing;
        double fill;
        double cz;
        double weight;
        int category; // 0: wall, 1: gas, 2: fluid
    };

    std::vector<MomentumLinkSample> momentumLinkSamples;

    std::vector<MomentumCellSample> momentumCellSamples;
    double momentumBeforeStreamingZLB = 0.0;

    dynamicsHandler.addSweeps(
        timeloop,
        [&](const char* stage)
        {
            if (auditedStep % diagnosticFrequency != uint_c(0))
                return;

            const Eigen::Vector3d currentMomentumLB =
                measureLiquidMomentumLB();

            if (std::string(stage) == "Audit: after boundary")
            {
                momentumCellSamples.clear();
                momentumLinkSamples.clear();
                momentumBeforeStreamingZLB = currentMomentumLB.z();

                for (auto& block : *blockForest)
                {
                    const auto* flags =
                        block.getData<FlagField_T>(flagFieldID);

                    const auto* fill =
                        block.getData<ScalarField_T>(fillFieldID);

                    for (const auto& cell : flags->xyzSize())
                    {
                        const auto flag = flags->get(
                            cell.x(), cell.y(), cell.z());

                        if (!flagInfo.isLiquid(flag) &&
                            !flagInfo.isInterface(flag))
                        {
                            continue;
                        }

                        momentumCellSamples.push_back({
                            &block,
                            cell,
                            double(fill->get(
                                cell.x(), cell.y(), cell.z()))
                        });

                        const auto* pdfs =
    block.getData<PdfField_T>(pdfFieldID);

const double fraction = double(fill->get(
    cell.x(), cell.y(), cell.z()));

for (auto d = LatticeModelStencil_T::beginNoCenter();
     d != LatticeModelStencil_T::end(); ++d)
{
    const auto neighborFlag = flags->get(
        cell.x() + d.cx(),
        cell.y() + d.cy(),
        cell.z() + d.cz());

    int category = 2;

    if (flagInfo.isObstacle(neighborFlag))
        category = 0;
    else if (flagInfo.isGas(neighborFlag))
        category = 1;

    momentumLinkSamples.push_back({
        &block,
        cell,
        uint_c(d.toInvIdx()),
        double(pdfs->get(
            cell.x(), cell.y(), cell.z(), d.toIdx())),
        fraction,
        double(d.cz()),
        double(LatticeModel_T::w[d.toIdx()]),
        category
    });
}
                    }
                }
            }
            else if (
                std::string(stage) ==
                "Audit: after streaming and reconstruction")
            {
                double oldFillNewMomentumZLB = 0.0;

                for (const auto& sample : momentumCellSamples)
                {
                    const auto* pdfs =
                        sample.block->getData<PdfField_T>(pdfFieldID);

                    const auto momentum = pdfs->getMomentumDensity(
                        sample.cell.x(),
                        sample.cell.y(),
                        sample.cell.z());

                    oldFillNewMomentumZLB +=
                        sample.fill * double(momentum[2]);
                }

                const double pdfContributionZLB =
                    oldFillNewMomentumZLB - momentumBeforeStreamingZLB;

                const double fillContributionZLB =
                    currentMomentumLB.z() - oldFillNewMomentumZLB;

                WALBERLA_LOG_INFO_ON_ROOT(
                    "Streaming momentum split: step=" << auditedStep
                    << " pdfContributionZLB=" << pdfContributionZLB
                    << " fillContributionZLB=" << fillContributionZLB);

                double linkDeltaPzLB[3] = {0.0, 0.0, 0.0};

for (const auto& sample : momentumLinkSamples)
{
    const auto* pdfs =
        sample.block->getData<PdfField_T>(pdfFieldID);

    const double incoming = double(pdfs->get(
        sample.cell.x(),
        sample.cell.y(),
        sample.cell.z(),
        sample.incomingIndex));

    // Ambient density is 1 in the current configuration.
    // This reference cancels when summed over all directions.
    const double reference = 2.0 * sample.weight;

    linkDeltaPzLB[sample.category] -=
        sample.fill *
        (incoming + sample.outgoing - reference) *
        sample.cz;
}

WALBERLA_LOG_INFO_ON_ROOT(
    "Streaming link budget: step=" << auditedStep
    << " wallDeltaPzLB=" << linkDeltaPzLB[0]
    << " gasDeltaPzLB=" << linkDeltaPzLB[1]
    << " internalDeltaPzLB=" << linkDeltaPzLB[2]
    << " closureErrorLB="
    << (
        linkDeltaPzLB[0]
        + linkDeltaPzLB[1]
        + linkDeltaPzLB[2]
        - pdfContributionZLB
    ));
            }

            const Eigen::Vector3d changeLB =
                currentMomentumLB - previousStageMomentumLB;

            WALBERLA_LOG_INFO_ON_ROOT(
                "Momentum stage: step=" << auditedStep
                << " stage=" << stage
                << " deltaPzLB=" << changeLB.z()
                << " PzLB=" << currentMomentumLB.z());

            previousStageMomentumLB = currentMomentumLB;
        });

   // add load balancing
   LoadBalancer< FlagField_T, CommunicationStencil_T, LatticeModelStencil_T > loadBalancer(
      blockForest, communication, pdfCommunication, bubbleModel, uint_c(50), uint_c(10), uint_c(5),
      loadBalancingFrequency, printLoadBalancingStatistics);
   timeloop.addFuncAfterTimeStep(loadBalancer, "Sweep: load balancing");

   // evaluate height of the drop
   const std::shared_ptr< real_t > dropHeight = std::make_shared< real_t >(real_c(0));
   const DropHeightEvaluator< FreeSurfaceBoundaryHandling_T > dropHeightEvaluator(
      blockForest, freeSurfaceBoundaryHandling, fillFieldID, dropHeight, evaluationFrequency);
   timeloop.addFuncAfterTimeStep(dropHeightEvaluator, "Evaluator: drop height");

   // add evaluator for total and excessive mass (mass that is currently undistributed)
    const auto totalMass =
     std::make_shared<real_t>(real_c(0));

    const auto excessMass =
        std::make_shared<real_t>(real_c(0));

    TotalMassComputer<
        FreeSurfaceBoundaryHandling_T,
        PdfField_T,
        FlagField_T,
        ScalarField_T
    > totalMassComputer(
        blockForest,
        freeSurfaceBoundaryHandling,
        pdfFieldID,
        fillFieldID,
        dynamicsHandler.getConstExcessMassFieldID(),
        uint_c(1),  // Update every LB step.
        totalMass,
        excessMass);

    // Measure the initial state before advancing LB.
    totalMassComputer();

    const double initialLiquidMassLB = double(*totalMass);

    if (!std::isfinite(initialLiquidMassLB) ||
        initialLiquidMassLB <= 0.0)
    {
        throw std::runtime_error(
            "Initial liquid mass must be finite and positive");
    }

    WALBERLA_LOG_INFO_ON_ROOT(
        "Initial liquid mass LB=" << initialLiquidMassLB
    );

    timeloop.addFuncAfterTimeStep(
        totalMassComputer,
        "Evaluator: total mass");

   // add VTK output
   addVTKOutput< LatticeModel_T, FreeSurfaceBoundaryHandling_T, PdfField_T, FlagField_T, ScalarField_T, VectorField_T >(
      blockForest, timeloop, walberlaEnv.config(), flagInfo, pdfFieldID, flagFieldID, fillFieldID, forceDensityFieldID,
      geometryHandler.getCurvatureFieldID(), geometryHandler.getNormalFieldID(),
      geometryHandler.getObstNormalFieldID());

   // add triangle mesh output of free surface
   CapillarySurfaceMeshWriter< ScalarField_T, FlagField_T > surfaceMeshWriter(
      blockForest, fillFieldID, flagFieldID, flagIDs::liquidInterfaceGasFlagIDs, real_c(0), walberlaEnv.config());
   surfaceMeshWriter(); // write initial mesh
   timeloop.addFuncAfterTimeStep(surfaceMeshWriter, "Writer: surface mesh");

    const uint_t particleWriteFrequency =
    walberlaEnv.config()
        ->getOneBlock("MeshOutputParameters")
        .getParameter<uint_t>("writeFrequency");

uint_t particleOutputCounter = 0;

auto writeParticleMesh = [&]()
{
    if (particleWriteFrequency == uint_c(0))
        return;

    const uint_t frame = particleOutputCounter++;

    if (frame % particleWriteFrequency != uint_c(0))
        return;

    const std::filesystem::path folder("particle-out");
    std::filesystem::create_directories(folder);

    std::ostringstream filename;
    filename << "particle_"
             << std::setfill('0') << std::setw(12)
             << frame << ".vtp";

    std::ofstream out(folder / filename.str());

    if (!out)
    {
        throw std::runtime_error(
            "Cannot open particle visualization output");
    }

    out << std::setprecision(17);

    const std::size_t triangleCount = solidSurface.size();
    const std::size_t pointCount = 3 * triangleCount;

    out << "<?xml version=\"1.0\"?>\n"
        << "<VTKFile type=\"PolyData\" version=\"0.1\" "
           "byte_order=\"LittleEndian\">\n"
        << "<PolyData>\n"
        << "<Piece NumberOfPoints=\"" << pointCount
        << "\" NumberOfVerts=\"0\" NumberOfLines=\"0\" "
           "NumberOfStrips=\"0\" NumberOfPolys=\""
        << triangleCount << "\">\n";

    out << "<Points>\n"
        << "<DataArray type=\"Float64\" "
           "NumberOfComponents=\"3\" format=\"ascii\">\n";

    for (const auto& triangle : solidSurface)
    {
        for (const auto& point : triangle)
        {
            out << point.x() << ' '
                << point.y() << ' '
                << point.z() << '\n';
        }
    }

    out << "</DataArray>\n</Points>\n"
        << "<Polys>\n"
        << "<DataArray type=\"Int64\" "
           "Name=\"connectivity\" format=\"ascii\">\n";

    for (std::size_t i = 0; i < triangleCount; ++i)
    {
        out << 3 * i << ' '
            << 3 * i + 1 << ' '
            << 3 * i + 2 << '\n';
    }

    out << "</DataArray>\n"
        << "<DataArray type=\"Int64\" "
           "Name=\"offsets\" format=\"ascii\">\n";

    for (std::size_t i = 0; i < triangleCount; ++i)
        out << 3 * (i + 1) << '\n';

    out << "</DataArray>\n</Polys>\n"
        << "</Piece>\n</PolyData>\n</VTKFile>\n";

    out.close();

    if (!out)
    {
        throw std::runtime_error(
            "Failed writing particle visualization output");
    }
};

// Match the liquid mesh writer's initial frame and cadence.
writeParticleMesh();

timeloop.addFuncAfterTimeStep(
    writeParticleMesh,
    "Writer: particle mesh");

   // add logging for computational performance
   const lbm::PerformanceLogger< FlagField_T > perfLogger(blockForest, flagFieldID, flagIDs::liquidInterfaceFlagIDs,
                                                          performanceLogFrequency);
   // timeloop.addFuncAfterTimeStep(perfLogger, "Evaluator: performance logging");

   WcTimingPool timingPool;

    Eigen::Vector3d forceSumSI =
        Eigen::Vector3d::Zero();

    Eigen::Vector3d torqueSumSI =
        Eigen::Vector3d::Zero();

    if (!dem.advanceDemStep ||
        !dem.applyRemapImpulseSI ||
        !dem.submitFluidLoadSI)
    {
        throw std::runtime_error(
            "Missing DEM coupling callbacks");
    }

   real_t formerDropHeight = real_c(0);
   for (uint_t t = uint_c(0); t != timesteps; ++t)
   {
       const Eigen::Vector3d momentumBeforeLB =
    measureLiquidMomentumLB();

       auditedStep = t;
       previousStageMomentumLB = momentumBeforeLB;

       timeloop.singleStep(timingPool, true);

       const Eigen::Vector3d fluidMomentumChangeLB =
           measureLiquidMomentumLB() - momentumBeforeLB;

       // Diagnostic reaction equivalent to the measured fluid change.
       // One LB timestep equals one in lattice units.
       const Eigen::Vector3d reactionEquivalentLB =
           -fluidMomentumChangeLB;

       const auto wallAngleLine = reconstructSurfaceLine(
                    solidSurface,
                    interfacePlanes,
                    solidCenterEigen,
                    double(surfaceTension),
                    3.0,
                    double(contactAngle)
                );

       const Eigen::Vector3d pressureForceLB(
            sphereForceLB[0],
            sphereForceLB[1],
            sphereForceLB[2]
        );

       const Eigen::Vector3d totalFluidForceLB =
           pressureForceLB + wallAngleLine.load.force;

       const Eigen::Vector3d totalFluidTorqueLB =
           sphereTorqueLB + wallAngleLine.load.torque;

       if (t % diagnosticFrequency == uint_c(0))
       {
           WALBERLA_LOG_INFO_ON_ROOT(
               "Fixed-particle force audit:"
               << " step=" << t
               << " fluidDeltaPzLB=" << fluidMomentumChangeLB.z()
               << " reactionEquivalentFzLB=" << reactionEquivalentLB.z()
               << " pressureFzLB=" << pressureForceLB.z()
               << " reconstructedLineFzLB=" << wallAngleLine.load.force.z()
               << " inferredLineFzLB="
               << reactionEquivalentLB.z() - pressureForceLB.z()
               << " balanceResidualZLB="
               << fluidMomentumChangeLB.z() + totalFluidForceLB.z()
               << " liquidWallFzLB=" << liquidWallFzLB
               << " interfaceWallFzLB=" << interfaceWallFzLB
               << " weightedInterfaceWallFzLB="
                << weightedInterfaceWallFzLB
                << " fillWeightedResidualZLB="
                << (
                    fluidMomentumChangeLB.z()
                    + liquidWallFzLB
                    + weightedInterfaceWallFzLB
                    + wallAngleLine.load.force.z()
                )
           );
       }

        const Eigen::Vector3d totalFluidForceSI =
            units.forceScaleSI * totalFluidForceLB;

       const Eigen::Vector3d totalFluidTorqueSI =
           units.torqueScaleSI * totalFluidTorqueLB;

       forceSumSI += totalFluidForceSI;
       torqueSumSI += totalFluidTorqueSI;

       if (t % diagnosticFrequency == uint_c(0))
       {
           if (contactLineCellCount > uint_c(0))
           {
               const real_t pi = std::acos(real_c(-1));

               const real_t contactLineZ =
                   contactLineZSum / real_c(contactLineCellCount);

               const real_t cosBeta = std::clamp(
                   (contactLineZ - solidCenter[2]) / solidRadius,
                   real_c(-1),
                   real_c(1)
               );

               const real_t beta = std::acos(cosBeta);
               const real_t theta =
                   contactAngle * pi / real_c(180);

               const real_t contactRadius =
                   solidRadius * std::sin(beta);

               const Eigen::Vector3d probe =
                solidCenterEigen + Eigen::Vector3d(
                    double(contactRadius),
                    0.0,
                    double(contactLineZ - solidCenter[2])
                );

               const auto fieldSample = sampleInterfaceField(
                   probe,
                   interfacePlanes,
                   3.0 // Support radius in lattice units.
               );

               const auto fittedLine = reconstructSurfaceLine(
                    solidSurface,
                    interfacePlanes,
                    solidCenterEigen,
                    double(surfaceTension),
                    3.0
                );

               const std::filesystem::path outputDirectory =
                    "contact-line-out";

               std::filesystem::create_directories(outputDirectory);

               std::ostringstream filename;
               filename << "contact_line_"
                        << std::setfill('0') << std::setw(12)
                        << t << ".vtp";

               writeContactLineVtp(
                   fittedLine,
                   outputDirectory / filename.str()
               );

               const real_t lineForceZ =
                   real_c(2) * pi
                   * contactRadius
                   * surfaceTension
                   * std::sin(beta + theta);

               const double liquidMassChangeLB =
                double(*totalMass) - initialLiquidMassLB;

               const double relativeLiquidMassChange =
                   liquidMassChangeLB / initialLiquidMassLB;

               WALBERLA_LOG_INFO_ON_ROOT(
                    "step=" << t
                    << " planes=" << interfacePlanes.size()
                    << " pressureFz=" << sphereForceLB[2]
                    << " ringLineFz=" << lineForceZ
                    << " ringTotalFz=" << sphereForceLB[2] + lineForceZ
                    << " fitSamples="
                    << (fieldSample ? fieldSample->sampleCount : 0)
                    << " fitValue="
                    << (fieldSample
                        ? std::to_string(fieldSample->value)
                        : std::string("unknown"))
                    << " fittedSegments=" << fittedLine.segments.size()
                    << " fittedLength=" << fittedLine.length
                    << " fittedLineFz=" << fittedLine.load.force.z()
                    << " fittedTotalFz="
                    << double(sphereForceLB[2]) + fittedLine.load.force.z()
                    << " badEnds=" << fittedLine.badEndpoints
                    << " unsupported=" << fittedLine.unsupportedTriangles
                    << " degenerate=" << fittedLine.degenerateTriangles
                    << " invalid=" << fittedLine.invalidSegments
                    << " wallLineFz=" << wallAngleLine.load.force.z()
                    << " wallTotalFz="
                    << double(sphereForceLB[2]) + wallAngleLine.load.force.z()
                    << " wallInvalid=" << wallAngleLine.invalidSegments
                    << " totalForceLB=" << totalFluidForceLB.transpose()
                    << " totalTorqueNormLB=" << totalFluidTorqueLB.norm()
                    << " forceSI_N=" << totalFluidForceSI.transpose()
                    << " torqueSI_Nm=" << totalFluidTorqueSI.transpose()
                    << " hydroTorqueZLB=" << sphereTorqueLB.z()
                    << " liquidMassLB=" << *totalMass
                    << " relativeMassChange=" << relativeLiquidMassChange
                    << " excessMassLB=" << *excessMass
                    << " wallMassCorrectionLB=" << movingWallMassCorrectionLB
                    << " dynamicsDeltaPzLB="
<< (
    measureLiquidMomentumLB().z()
    - momentumAfterGeometryLB.z()
)
<< " beforeDynamicsDeltaPzLB="
<< (
    fluidMomentumChangeLB.z()
    - (
        measureLiquidMomentumLB().z()
        - momentumAfterGeometryLB.z()
    )
)
                );
           }
           else
           {
               WALBERLA_LOG_INFO_ON_ROOT(
                   "step=" << t
                   << " pressureFz=" << sphereForceLB[2]
                   << " ring estimate unavailable"
               );
           }
       }

       if ((t + uint_c(1)) % units.lbStepsPerDemStep == uint_c(0))
{
    const double count = double(units.lbStepsPerDemStep);

    // Apply the mean fluid load over this DEM interval.
    dem.submitFluidLoadSI(
        forceSumSI / count,
        torqueSumSI / count);

    forceSumSI.setZero();
    torqueSumSI.setZero();

    // Advance the actual DEM integrator.
           // Diagnostic: hold DEM fixed while checking the fluid-force balance.
           // dem.advanceDemStep();

    const Eigen::Vector3d newCenterLB =
        demPositionToLb(dem.particlePositionSI());

    const Eigen::Vector3d displacementLB =
        newCenterLB - solidCenterEigen;

    if (!newCenterLB.allFinite() ||
        displacementLB.norm() > 0.1)
    {
        throw std::runtime_error(
            "DEM displacement exceeds 0.1 lattice cells "
            "per coupling interval; reduce the DEM timestep");
    }

    RemapResult result{
        Eigen::Vector3d::Zero(),
        Eigen::Vector3d::Zero(),
        solidCenterEigen
    };

    const auto changes =
        detectParticleCellChanges(newCenterLB);

    if (changes.covered != 0 || changes.uncovered != 0)
    {
        // Particle motion changes the occupied lattice cells.
        result = remapParticle(newCenterLB);

        communication();
        pdfCommunication();
    }
    else
    {
        // Subcell motion: update the analytic geometry without
        // modifying fluid cells or rebuilding their PDFs.
        solidCenterEigen = newCenterLB;

        solidCenter = Vector3<real_t>(
            real_c(newCenterLB.x()),
            real_c(newCenterLB.y()),
            real_c(newCenterLB.z()));

        solidSurface = makeSphereSurface(
            solidCenterEigen,
            double(solidRadius),
            4);
    }

    // The remapper reports angular impulse about the OLD center.
    // Convert it to angular impulse about the CURRENT center.
    const Eigen::Vector3d angularImpulseAtNewCenterLB =
        result.angularImpulseAboutOldCenterLB -
        (newCenterLB - result.oldCenterLB)
            .cross(result.linearImpulseLB);

    const Eigen::Vector3d linearImpulseSI =
        (units.forceScaleSI * units.dtSI) *
        result.linearImpulseLB;

    const Eigen::Vector3d angularImpulseSI =
        (units.torqueScaleSI * units.dtSI) *
        angularImpulseAtNewCenterLB;

    dem.applyRemapImpulseSI(
        linearImpulseSI,
        angularImpulseSI);

    // Refresh wall velocities after DEM stepping and impulses.
    const Eigen::Vector3d velocityLB =
        dem.particleVelocitySI() *
        (units.dtSI / units.dxSI);

    const Eigen::Vector3d angularVelocityLB =
        dem.particleAngularVelocitySI() * units.dtSI;

    if (!velocityLB.allFinite() ||
        !angularVelocityLB.allFinite())
    {
        throw std::runtime_error(
            "Non-finite DEM velocity after coupling");
    }

    updateSphereWallVelocity(
        velocityLB,
        angularVelocityLB);

    // Include any remapping in the current mass measurement.
    totalMassComputer();

           if (t % diagnosticFrequency == uint_c(0))
{
    double fluidMassLB = 0.0;
    double fluidMassWeightedZLB = 0.0;

    Eigen::Vector3d fluidMomentumLB =
        Eigen::Vector3d::Zero();

    for (auto& block : *blockForest)
    {
        const auto* flags =
            block.getData<FlagField_T>(flagFieldID);

        const auto* fill =
            block.getData<ScalarField_T>(fillFieldID);

        const auto* pdfs =
            block.getData<PdfField_T>(pdfFieldID);

        for (const auto& cell : flags->xyzSize())
        {
            const auto flag =
                flags->get(cell.x(), cell.y(), cell.z());

            if (!flagInfo.isLiquid(flag) &&
                !flagInfo.isInterface(flag))
                continue;

            const double fraction = double(
                fill->get(cell.x(), cell.y(), cell.z()));

            const double massLB = fraction * double(
                pdfs->getDensity(cell.x(), cell.y(), cell.z()));

            const auto position =
                blockForest->getBlockLocalCellCenter(block, cell);

            const auto momentum = pdfs->getMomentumDensity(
                cell.x(), cell.y(), cell.z());

            fluidMassLB += massLB;
            fluidMassWeightedZLB += massLB * double(position[2]);

            fluidMomentumLB += fraction * Eigen::Vector3d{
                static_cast<double>(momentum[0]),
                static_cast<double>(momentum[1]),
                static_cast<double>(momentum[2])
            };
        }
    }

    const Eigen::Vector3d fluidMomentumSI =
        (units.forceScaleSI * units.dtSI) * fluidMomentumLB;

    const Eigen::Vector3d particleMomentumSI =
        dem.particleMassSI * dem.particleVelocitySI();

    const Eigen::Vector3d combinedMomentumSI =
        fluidMomentumSI + particleMomentumSI;

    WALBERLA_LOG_INFO_ON_ROOT(
        "Momentum balance:"
        << " particleZLB=" << solidCenterEigen.z()
        << " fluidCOMzLB="
        << fluidMassWeightedZLB / fluidMassLB
        << " particlePzSI=" << particleMomentumSI.z()
        << " fluidPzSI=" << fluidMomentumSI.z()
        << " combinedPzSI=" << combinedMomentumSI.z()
    );
}
}

      // check convergence
      if (t % evaluationFrequency == uint_c(0))
      {
         // WALBERLA_LOG_DEVEL_ON_ROOT("time step = " << t << "\n\t\tdrop height = " << *dropHeight)
         if (std::abs(formerDropHeight - *dropHeight) / *dropHeight < convergenceThreshold)
         {
            WALBERLA_LOG_DEVEL_ON_ROOT("Final converged drop height=" << *dropHeight);
            return EXIT_SUCCESS;
         }
         formerDropHeight = *dropHeight;
      }

      // if (t % performanceLogFrequency == uint_c(0) && t > uint_c(0)) { timingPool.logResultOnRoot(); }
   }

   WALBERLA_LOG_DEVEL_ON_ROOT("Final non-converged drop height=" << *dropHeight);

   return EXIT_SUCCESS;
}
} // namespace DropWetting
} // namespace free_surface
} // namespace walberla

int runCapillaryLb(
    int argc,
    char** argv,
    const DemLbParameters& dem
) {
    return walberla::free_surface::DropWetting::runWithDem(
        argc, argv, dem
    );
}
