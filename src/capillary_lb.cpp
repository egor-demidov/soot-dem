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
//! \file capillary_lb.cpp
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//
// Single-particle DEM/free-surface-LB prototype. The continuum contact-line
// estimate is retained; an explicit PDF momentum projection repairs its
// mismatch with the native free-surface update. Conservation is enforced
// numerically and must not be interpreted as validation of capillary-force
// accuracy.
//======================================================================================================================

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <memory>
#include <queue>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "blockforest/Initialization.h"

#include "core/Environment.h"

#include "lbm/field/AddToStorage.h"
#include "lbm/free_surface/LoadBalancing.h" // createNonUniformBlockForest
#include "lbm/free_surface/SurfaceMeshWriter.h"
#include "lbm/free_surface/TotalMassComputer.h"
#include "lbm/free_surface/VtkWriter.h"
#include "lbm/free_surface/bubble_model/Geometry.h"
#include "lbm/free_surface/dynamics/SurfaceDynamicsHandler.h"
#include "lbm/free_surface/surface_geometry/SurfaceGeometryHandler.h"
#include "lbm/free_surface/surface_geometry/Utility.h"
#include "lbm/lattice_model/D3Q19.h"

#include "contact_line_force.h"
#include "contact_line_output.h"
#include "dem_lb_parameters.h"
#include "interface_plane_field.h"
#include "lb_units.h"
#include "capillary_setup.h"
#include "sphere_surface.h"
#include "surface_contact_line.h"

#include "momentum_projection.h"

namespace walberla::free_surface {
template <typename ScalarField_T, typename FlagField_T> class CapillarySurfaceMeshWriter {
  public:
    CapillarySurfaceMeshWriter(const std::weak_ptr<StructuredBlockForest> &blockForest,
                               const ConstBlockDataID &fillFieldID, const ConstBlockDataID &flagFieldID,
                               const Set<FlagUID> &liquidInterfaceGasFlagIDSet, real_t obstacleFillLevel,
                               uint_t writeFrequency, const std::string &baseFolder)
        : blockForest_(blockForest), fillFieldID_(fillFieldID), flagFieldID_(flagFieldID),
          liquidInterfaceGasFlagIDSet_(liquidInterfaceGasFlagIDSet), obstacleFillLevel_(obstacleFillLevel),
          writeFrequency_(writeFrequency), baseFolder_(baseFolder), executionCounter_(uint_c(0)),
          fillFieldCloneID_((abortIfNullptr::abortIfNullptr(blockForest),
                             field::addCloneToStorage<ScalarField_T>(blockForest_.lock(), fillFieldID_,
                                                                     "Fill level field clone"))) {}

    // config block must be named "MeshOutputParameters"
    CapillarySurfaceMeshWriter(const std::weak_ptr<StructuredBlockForest> &blockForest,
                               const ConstBlockDataID &fillFieldID, const ConstBlockDataID &flagFieldID,
                               const Set<FlagUID> &liquidInterfaceGasFlagIDSet, real_t obstacleFillLevel,
                               const std::weak_ptr<Config> &config)
        : CapillarySurfaceMeshWriter(
              blockForest, fillFieldID, flagFieldID, liquidInterfaceGasFlagIDSet, obstacleFillLevel,
              (abortIfNullptr::abortIfNullptr(config),
               config.lock()->getOneBlock("MeshOutputParameters").getParameter<uint_t>("writeFrequency")),
              (abortIfNullptr::abortIfNullptr(config),
               config.lock()->getOneBlock("MeshOutputParameters").getParameter<std::string>("baseFolder"))) {}

    void operator()() {
        if (writeFrequency_ == uint_c(0)) {
            return;
        }

        if (executionCounter_ == uint_c(0)) {
            createBaseFolder();
            writeMesh();
        } else {
            writeMesh();
        }

        ++executionCounter_;
    }

  private:
    void createBaseFolder() const {
        WALBERLA_ROOT_SECTION() {
            const std::filesystem::path basePath(baseFolder_);
            std::filesystem::create_directories(basePath);
            for (const auto &entry : std::filesystem::directory_iterator(basePath)) {
                const auto name = entry.path().filename().string();
                if (entry.is_regular_file() && name.starts_with("simulation_step_") &&
                    entry.path().extension() == ".vtp")
                    std::filesystem::remove(entry.path());
            }
        }
        WALBERLA_MPI_BARRIER();
    }

    void writeMesh() {
        // only write mesh in given frequency
        if (executionCounter_ % writeFrequency_ != uint_c(0)) {
            return;
        }

        // rank=0 is just an arbitrary choice here
        const int targetRank = 0;

        auto blockForest = blockForest_.lock();
        WALBERLA_CHECK_NOT_NULLPTR(blockForest);

        // update clone of fill level field and set fill level of all non-liquid,
        // -interface, or -gas cells to zero
        updateFillFieldClone(blockForest);

        const auto surfaceMesh = postprocessing::realFieldToSurfaceMesh<ScalarField_T>(
            blockForest, fillFieldCloneID_, real_c(0.5), uint_c(0), true, targetRank, MPI_COMM_WORLD);

        surfaceMesh->translate(Vector3<real_t>(real_c(0.5) * blockForest->dx(),
                                               real_c(0.5) * blockForest->dy(),
                                               real_c(0.5) * blockForest->dz()));

        WALBERLA_EXCLUSIVE_WORLD_SECTION(targetRank) {
            std::ostringstream filename;
            filename << "simulation_step_" << std::setfill('0') << std::setw(12) << executionCounter_
                     << ".vtp";

            const auto path = std::filesystem::path(baseFolder_) / filename.str();

            std::ofstream output(path);
            if (!output.is_open()) {
                WALBERLA_ABORT("Cannot open surface mesh output file.");
            }

            geometry::writeMeshVtp(output, *surfaceMesh);
        }
    }

    // update clone of fill level field and set fill level of all non-liquid,
    // -interface, or -gas cells to zero; explicitly use shared_ptr instead of
    // weak_ptr to avoid checking the latter's validity (is done in writeMesh()
    // already)
    void updateFillFieldClone(const std::shared_ptr<StructuredBlockForest> &blockForest) {
        for (auto block = blockForest->begin(); block != blockForest->end(); ++block) {
            auto *output = block->template getData<ScalarField_T>(fillFieldCloneID_);

            const auto *input = block->template getData<ScalarField_T>(fillFieldID_);

            const auto *flags = block->template getData<FlagField_T>(flagFieldID_);

            const auto fluidMask = flags->getMask(liquidInterfaceGasFlagIDSet_);

            const cell_idx_t nx = cell_idx_c(input->xSize());
            const cell_idx_t ny = cell_idx_c(input->ySize());
            const cell_idx_t nz = cell_idx_c(input->zSize());

            const cell_idx_t inputGhost = cell_idx_c(input->nrOfGhostLayers());
            const cell_idx_t flagGhost = cell_idx_c(flags->nrOfGhostLayers());

            const cell_idx_t ng = inputGhost < flagGhost ? inputGhost : flagGhost;

            for (const auto &cell : output->xyzSizeWithGhostLayer()) {
                real_t sum = real_c(0);
                uint_t validCount = uint_c(0);

                // Value at a grid corner: average its eight adjacent cells.
                for (cell_idx_t dz = -1; dz <= 0; ++dz)
                    for (cell_idx_t dy = -1; dy <= 0; ++dy)
                        for (cell_idx_t dx = -1; dx <= 0; ++dx) {
                            const cell_idx_t x = cell.x() + dx;
                            const cell_idx_t y = cell.y() + dy;
                            const cell_idx_t z = cell.z() + dz;

                            if (x < -ng || x >= nx + ng || y < -ng || y >= ny + ng || z < -ng ||
                                z >= nz + ng) {
                                continue;
                            }

                            // Exclude solid bookkeeping fill values from the mesh.
                            if ((flags->get(x, y, z) & fluidMask) != 0) {
                                sum += input->get(x, y, z);
                                ++validCount;
                            }
                        }

                output->get(cell.x(), cell.y(), cell.z()) =
                    validCount > uint_c(0) ? sum / real_c(validCount) : real_c(0);
            }
        }
    }

    std::weak_ptr<StructuredBlockForest> blockForest_;
    ConstBlockDataID fillFieldID_;
    ConstBlockDataID flagFieldID_;
    Set<FlagUID> liquidInterfaceGasFlagIDSet_;
    real_t obstacleFillLevel_;
    uint_t writeFrequency_;
    std::string baseFolder_;
    uint_t executionCounter_;

    BlockDataID fillFieldCloneID_;
}; // class SurfaceMeshWriter
} // namespace walberla::free_surface

namespace walberla {
namespace free_surface {
namespace DropWetting {
using ScalarField_T = GhostLayerField<real_t, 1>;
using VectorField_T = GhostLayerField<Vector3<real_t>, 1>;

using CollisionModel_T = lbm::collision_model::SRT;
using ForceModel_T = lbm::force_model::GuoField<VectorField_T>;
using LatticeModel_T = lbm::D3Q19<CollisionModel_T, true, ForceModel_T, 2>;
using LatticeModelStencil_T = LatticeModel_T::Stencil;
using PdfField_T = lbm::PdfField<LatticeModel_T>;
using PdfCommunication_T = blockforest::SimpleCommunication<LatticeModelStencil_T>;

// the geometry computations in SurfaceGeometryHandler require meaningful values
// in the ghost layers in corner directions (flag field and fill level field);
// this holds, even if the lattice model uses a D3Q19 stencil
using CommunicationStencil_T =
    std::conditional_t<LatticeModel_T::Stencil::D == uint_t{2}, stencil::D2Q9, stencil::D3Q27>;
using Communication_T = blockforest::SimpleCommunication<CommunicationStencil_T>;

using flag_t = uint32_t;
using FlagField_T = FlagField<flag_t>;
using FreeSurfaceBoundaryHandling_T = FreeSurfaceBoundaryHandling<LatticeModel_T, FlagField_T, ScalarField_T>;

int runWithDem(int argc, char **argv, const DemLbParameters &dem) {
    Environment walberlaEnv(argc, argv);

    if (argc < 2) {
        WALBERLA_ABORT("Please specify a parameter file as input argument.")
    }

    const LbUnits units = makeLbUnits(dem);

    if (!dem.particlePositionSI) {
        throw std::runtime_error("DEM particle-position callback is missing");
    }

    if (dem.particleCount == 0 || !dem.setup) throw std::runtime_error("Missing input geometry");
    // Use particle 0's initial position as the common reference.
    const Eigen::Vector3d demReferencePositionSI = dem.particlePositionSI(0);

    if (!dem.particleVelocitySI || !dem.particleAngularVelocitySI) {
        throw std::runtime_error("DEM velocity callbacks are missing");
    }

    const Eigen::Vector3d initialVelocityLB = dem.particleVelocitySI(0) * (units.dtSI / units.dxSI);

    const Eigen::Vector3d initialAngularVelocityLB = dem.particleAngularVelocitySI(0) * units.dtSI;

    WALBERLA_LOG_INFO_ON_ROOT("DEM initial velocity [m/s] = "
                              << dem.particleVelocitySI(0).transpose()
                              << ", velocity LB = " << initialVelocityLB.transpose()
                              << ", angular velocity LB = " << initialAngularVelocityLB.transpose());

    WALBERLA_LOG_INFO_ON_ROOT("Derived LB scales:" << " dxSI=" << units.dxSI << " dtSI=" << units.dtSI
                                                   << " LBstepsPerDEMstep=" << units.lbStepsPerDemStep
                                                   << " viscosityLB=" << units.viscosityLB
                                                   << " omega=" << units.relaxationRate
                                                   << " surfaceTensionLB=" << units.surfaceTensionLB
                                                   << " forceScaleSI=" << units.forceScaleSI
                                                   << " torqueScaleSI=" << units.torqueScaleSI);

    // get block forest parameters from parameter file
    auto blockForestParameters = walberlaEnv.config()->getOneBlock("BlockForestParameters");
    const Vector3<uint_t> cellsPerBlock =
        blockForestParameters.getParameter<Vector3<uint_t>>("cellsPerBlock");
    const Vector3<bool> periodicity = blockForestParameters.getParameter<Vector3<bool>>("periodicity");
    const uint_t loadBalancingFrequency =
        blockForestParameters.getParameter<uint_t>("loadBalancingFrequency");

    // read domain parameters from parameter file
    const auto domainParameters = walberlaEnv.config()->getOneBlock("DomainParameters");
    const real_t dropDiameter = domainParameters.getParameter<real_t>("dropDiameter");
    const Vector3<real_t> domainSizeFactor =
        domainParameters.getParameter<Vector3<real_t>>("domainSizeFactor");

    // define domain size
    Vector3<uint_t> domainSize = domainSizeFactor * dropDiameter;
    domainSize[0] = uint_c(domainSizeFactor[0] * dropDiameter);
    domainSize[1] = uint_c(domainSizeFactor[1] * dropDiameter);
    domainSize[2] = uint_c(domainSizeFactor[2] * dropDiameter);

    // compute number of blocks as defined by domainSize and cellsPerBlock
    Vector3<uint_t> numBlocks;
    numBlocks[0] = uint_c(std::ceil(real_c(domainSize[0]) / real_c(cellsPerBlock[0])));
    numBlocks[1] = uint_c(std::ceil(real_c(domainSize[1]) / real_c(cellsPerBlock[1])));
    numBlocks[2] = uint_c(std::ceil(real_c(domainSize[2]) / real_c(cellsPerBlock[2])));

    // get number of (MPI) processes
    uint_t numProcesses = uint_c(MPIManager::instance()->numProcesses());
    WALBERLA_CHECK_LESS_EQUAL(numProcesses, numBlocks[0] * numBlocks[1] * numBlocks[2],
                              "The number of MPI processes is greater than the number of blocks as "
                              "defined by "
                              "\"domainSize/cellsPerBlock\". This would result in unused MPI "
                              "processes. Either decrease "
                              "the number of MPI processes or increase \"cellsPerBlock\".")

    Vector3<uint_t> realDomainSize;
    realDomainSize[0] = cellsPerBlock[0] * numBlocks[0];
    realDomainSize[1] = cellsPerBlock[1] * numBlocks[1];
    realDomainSize[2] = cellsPerBlock[2] * numBlocks[2];

    if (domainSize[0] != realDomainSize[0] && periodicity[0]) {
        WALBERLA_ABORT("The specified domain size in x-direction can not be "
                       "obtained with the number of "
                       "blocks you specified.")
    }
    if (domainSize[1] != realDomainSize[1] && periodicity[1]) {
        WALBERLA_ABORT("The specified domain size in y-direction can not be "
                       "obtained with the number of "
                       "blocks you specified.")
    }
    // the z-direction must be no slip in this setup and is simply extended (see
    // below) to obtain the specified size
    int boundaryThicknessZ = int_c(realDomainSize[2]) - int_c(domainSize[2]);
    if (boundaryThicknessZ < 0) {
        WALBERLA_ABORT("Something went wrong: the resulting domain size in "
                       "z-direction is less than specified.")
    }

    // read physics parameters from parameter file
    const auto physicsParameters = walberlaEnv.config()->getOneBlock("PhysicsParameters");
    const real_t relaxationRate = real_c(units.relaxationRate);
    const real_t surfaceTension = real_c(units.surfaceTensionLB);
    const Vector3<real_t> force = physicsParameters.getParameter<Vector3<real_t>>("force");

    const bool enableWetting = physicsParameters.getParameter<bool>("enableWetting");
    const real_t contactAngle = enableWetting ? real_c(dem.liquidContactAngleDeg) : real_c(90);

    WALBERLA_LOG_INFO_ON_ROOT("Effective contact angle [deg] = " << contactAngle
                                                                 << ", enableWetting = " << enableWetting);

    const bool advanceDEM = physicsParameters.getParameter<bool>("advanceDEM", true);
    const uint_t diagnosticFrequency =
        physicsParameters.getParameter<uint_t>("diagnosticFrequency", uint_c(50));
    const double maxCorrectionSpeedLB = physicsParameters.getParameter<double>("maxCorrectionSpeedLB", 0.05);
    if (!std::isfinite(maxCorrectionSpeedLB) || maxCorrectionSpeedLB <= 0.0)
        throw std::runtime_error("maxCorrectionSpeedLB must be finite and positive");
    if (numProcesses != uint_c(1) || numBlocks[0] != uint_c(1) || numBlocks[1] != uint_c(1) ||
        numBlocks[2] != uint_c(1) || loadBalancingFrequency != uint_c(0))
        throw std::runtime_error("Conservation projection requires one process, "
                                 "one block, and no load balancing");
    if (force.sqrLength() != real_c(0))
        throw std::runtime_error("Conservation projection currently requires zero "
                                 "external acceleration");

    const CollisionModel_T collisionModel = CollisionModel_T(relaxationRate);

    const uint_t timesteps = physicsParameters.getParameter<uint_t>("timesteps");

    if (advanceDEM && timesteps % units.lbStepsPerDemStep != uint_c(0))
        throw std::runtime_error("timesteps must contain a whole number of DEM coupling intervals");
    if (!demReferencePositionSI.allFinite() || !initialVelocityLB.allFinite() ||
        !initialAngularVelocityLB.allFinite())
        throw std::runtime_error("Non-finite initial DEM state");

    // read model parameters from parameter file
    const auto modelParameters = walberlaEnv.config()->getOneBlock("ModelParameters");
    const std::string pdfReconstructionModel =
        modelParameters.getParameter<std::string>("pdfReconstructionModel");
    const std::string pdfRefillingModel = modelParameters.getParameter<std::string>("pdfRefillingModel");
    const std::string excessMassDistributionModel =
        modelParameters.getParameter<std::string>("excessMassDistributionModel");
    const std::string curvatureModel = modelParameters.getParameter<std::string>("curvatureModel");
    const bool useSimpleMassExchange = modelParameters.getParameter<bool>("useSimpleMassExchange");
    const bool enableBubbleModel = modelParameters.getParameter<bool>("enableBubbleModel");
    const bool enableBubbleSplits = modelParameters.getParameter<bool>("enableBubbleSplits");
    const real_t cellConversionThreshold = modelParameters.getParameter<real_t>("cellConversionThreshold");
    const real_t cellConversionForceThreshold =
        modelParameters.getParameter<real_t>("cellConversionForceThreshold");

    // create non-uniform block forest (non-uniformity required for load
    // balancing)
    const std::shared_ptr<StructuredBlockForest> blockForest =
        createNonUniformBlockForest(domainSize, cellsPerBlock, numBlocks, periodicity);

    // add force field
    const BlockDataID forceDensityFieldID =
        field::addToStorage<VectorField_T>(blockForest, "Force field", force, field::fzyx, uint_c(1));

    // create lattice model
    const LatticeModel_T latticeModel = LatticeModel_T(collisionModel, ForceModel_T(forceDensityFieldID));

    // add pdf field
    const BlockDataID pdfFieldID =
        lbm::addPdfFieldToStorage(blockForest, "PDF field", latticeModel, field::fzyx);

    // add fill level field (initialized with 0, i.e., gas everywhere)
    const BlockDataID fillFieldID = field::addToStorage<ScalarField_T>(blockForest, "Fill level field",
                                                                       real_c(0.0), field::fzyx, uint_c(2));

    // One moving particle, initially anchored at the domain center.
    const real_t solidRadius = real_c(0.5) * real_c(dem.cellsPerParticleDiameter);

    // Map the initial DEM position onto the existing LB placement.
    auto demPositionToLb = [&](const Eigen::Vector3d &positionSI) -> Eigen::Vector3d {
        return (positionSI - dem.setup->lower) / units.dxSI;
    };

    struct ParticleLbState {
        Eigen::Vector3d center;
        std::vector<SphereTriangle> surface;
        int surfaceLevel = 4;
    };

    std::vector<ParticleLbState> particles;
    particles.reserve(dem.particleCount);

    for (std::size_t i = 0; i < dem.particleCount; ++i) {
        const Eigen::Vector3d centerLB = demPositionToLb(dem.particlePositionSI(i));

        for (int axis = 0; axis < 3; ++axis) {
            const double margin = double(solidRadius) + 2.0;

            if (centerLB[axis] <= margin || centerLB[axis] >= double(domainSize[axis]) - margin) {
                throw std::runtime_error("Initial particle is too close to a domain boundary");
            }
        }

        particles.push_back({centerLB, makeSphereSurface(centerLB, double(solidRadius), 4)});

        WALBERLA_LOG_INFO_ON_ROOT("Initial particle " << i << " center LB = " << centerLB.transpose());
    }

    // Existing single-particle code still uses these during the refactor.
    Eigen::Vector3d solidCenterEigen = particles.at(0).center;

    Vector3<real_t> solidCenter(real_c(solidCenterEigen.x()), real_c(solidCenterEigen.y()),
                                real_c(solidCenterEigen.z()));

    // add boundary handling
    const std::shared_ptr<FreeSurfaceBoundaryHandling_T> freeSurfaceBoundaryHandling =
        std::make_shared<FreeSurfaceBoundaryHandling_T>(blockForest, pdfFieldID, fillFieldID);
    const BlockDataID flagFieldID = freeSurfaceBoundaryHandling->getFlagFieldID();
    const typename FreeSurfaceBoundaryHandling_T::FlagInfo_T &flagInfo =
        freeSurfaceBoundaryHandling->getFlagInfo();

    // Subcell sampling of the requested liquid region minus particle interiors.
    for (auto &block : *blockForest) {
        auto *fill = block.getData<ScalarField_T>(fillFieldID);
        double liquidVolume = 0;
        for (const auto &cell : fill->xyzSize()) {
            const auto p = blockForest->getBlockLocalCellCenter(block,cell);
            const Eigen::Vector3d position{double(p[0]),double(p[1]),double(p[2])};
            const auto centerSI = dem.setup->lower + units.dxSI * position;
            const double value = dem.setup->solid(centerSI) ? 0.0 : dem.setup->fill(centerSI);
            fill->get(cell.x(),cell.y(),cell.z()) = real_c(value);
            liquidVolume += value;
        }
        if (!(liquidVolume>0)) throw std::runtime_error("Initialized liquid is empty at this lattice resolution");
        WALBERLA_LOG_INFO_ON_ROOT("Initialized liquid volume SI=" << liquidVolume*std::pow(units.dxSI,3));
    }

    // initialize boundary conditions from config file
    const auto boundaryParameters = walberlaEnv.config()->getOneBlock("BoundaryParameters");
    freeSurfaceBoundaryHandling->initFromConfig(boundaryParameters);

    auto solidSurface = makeSphereSurface(solidCenterEigen, double(solidRadius), 4);

    struct InterfacePlaneSample {
        Eigen::Vector3d point;
        Eigen::Vector3d normal; // Points from liquid toward gas.
    };

    std::vector<InterfacePlaneSample> interfacePlanes;

    using BoundaryHandling_T = FreeSurfaceBoundaryHandling_T::BoundaryHandling_T;
    using UBB_T = FreeSurfaceBoundaryHandling_T::UBB_T;

    auto particleOwnerAt = [&](const Eigen::Vector3d &position) -> std::size_t {
        std::size_t owner = particles.size();
        double nearestDistanceSquared = double(solidRadius) * double(solidRadius);

        for (std::size_t i = 0; i < particles.size(); ++i) {
            const double distanceSquared = (position - particles[i].center).squaredNorm();

            if (distanceSquared <= nearestDistanceSquared) {
                owner = i;
                nearestDistanceSquared = distanceSquared;
            }
        }

        return owner;
    };

    auto isSolidAfterMove = [&](const Eigen::Vector3d &position, std::size_t movingParticle,
                                const Eigen::Vector3d &proposedCenter) -> bool {
        const double radiusSquared = double(solidRadius) * double(solidRadius);

        for (std::size_t i = 0; i < particles.size(); ++i) {
            const Eigen::Vector3d &center = i == movingParticle ? proposedCenter : particles[i].center;

            if ((position - center).squaredNorm() <= radiusSquared) {
                return true;
            }
        }

        return false;
    };

    for (auto &block : *blockForest) {
        auto *handling = block.getData<BoundaryHandling_T>(freeSurfaceBoundaryHandling->getHandlingID());

        auto *flags = handling->getFlagField();

        for (const auto &cell : flags->xyzSizeWithGhostLayer()) {
            const auto p = blockForest->getBlockLocalCellCenter(block, cell);

            const Eigen::Vector3d position{double(p[0]), double(p[1]), double(p[2])};

            const std::size_t owner = particleOwnerAt(position);

            if (owner == particles.size()) {
                continue;
            }

            Eigen::Vector3d velocityLB = Eigen::Vector3d::Zero();
            Eigen::Vector3d angularVelocityLB = Eigen::Vector3d::Zero();

            if (advanceDEM) {
                velocityLB = dem.particleVelocitySI(owner) * (units.dtSI / units.dxSI);

                angularVelocityLB = dem.particleAngularVelocitySI(owner) * units.dtSI;
            }

            const Eigen::Vector3d wallVelocityLB =
                velocityLB + angularVelocityLB.cross(position - particles[owner].center);

            handling->forceBoundary(FreeSurfaceBoundaryHandling_T::ubbFlagID, cell.x(), cell.y(), cell.z(),
                                    UBB_T::Velocity(real_c(wallVelocityLB.x()), real_c(wallVelocityLB.y()),
                                                    real_c(wallVelocityLB.z())));
        }
    }

    // IMPORTANT REMARK: this must be only called after every solid flag has been
    // set; otherwise, the boundary handling might not detect solid flags
    // correctly
    freeSurfaceBoundaryHandling->initFlagsFromFillLevel();

    ConstBlockDataID remapExcessMassFieldID = ConstBlockDataID();

    auto cellCenterEigen = [&](IBlock &block, const Cell &cell) -> Eigen::Vector3d {
        const auto p = blockForest->getBlockLocalCellCenter(block, cell);

        return Eigen::Vector3d{static_cast<double>(p[0]), static_cast<double>(p[1]),
                               static_cast<double>(p[2])};
    };

    struct ParticleChangedCell {
        IBlock *block;
        Cell cell;
    };

    struct CoveredLiquidCell {
        IBlock *block;
        Cell cell;

        double massLB;
        Eigen::Vector3d momentumLB;
    };

    struct ParticleCellChanges {
        uint_t covered = 0;
        uint_t uncovered = 0;
        uint_t retained = 0;

        std::vector<ParticleChangedCell> coveredCells;
        std::vector<ParticleChangedCell> uncoveredCells;

        double coveredLiquidVolumeLB = 0.0;
        double coveredLiquidMassLB = 0.0;

        std::vector<CoveredLiquidCell> coveredLiquidCells;

        Eigen::Vector3d coveredLiquidMomentumLB = Eigen::Vector3d::Zero();
    };

    auto detectParticleCellChanges = [&](std::size_t movingParticle,
                                         const Eigen::Vector3d &proposedCenterLB) {
        if (movingParticle >= particles.size()) {
            throw std::out_of_range("Invalid moving particle index");
        }

        ParticleCellChanges changes;

        for (auto &block : *blockForest) {
            const auto *flags = block.getData<FlagField_T>(flagFieldID);

            const auto *fill = block.getData<ScalarField_T>(fillFieldID);

            const auto *pdfs = block.getData<PdfField_T>(pdfFieldID);

            const ScalarField_T *storedExcess = nullptr;

            if (remapExcessMassFieldID != ConstBlockDataID()) {
                storedExcess = block.getData<ScalarField_T>(remapExcessMassFieldID);
            }

            const auto particleMask = flags->getFlag(FreeSurfaceBoundaryHandling_T::ubbFlagID);

            // Count owned cells only; ghost cells would duplicate counts.
            for (const auto &cell : flags->xyzSize()) {
                const bool wasParticle = (flags->get(cell.x(), cell.y(), cell.z()) & particleMask) != 0;

                const auto position = blockForest->getBlockLocalCellCenter(block, cell);

                const Eigen::Vector3d positionLB{double(position[0]), double(position[1]),
                                                 double(position[2])};

                const bool becomesParticle = isSolidAfterMove(positionLB, movingParticle, proposedCenterLB);

                if (!wasParticle && becomesParticle) {
                    ++changes.covered;
                    changes.coveredCells.push_back({&block, cell});

                    const auto cellFlag = flags->get(cell.x(), cell.y(), cell.z());

                    if (flagInfo.isLiquid(cellFlag) || flagInfo.isInterface(cellFlag)) {
                        const double fraction = double(fill->get(cell.x(), cell.y(), cell.z()));

                        const double density = double(pdfs->getDensity(cell.x(), cell.y(), cell.z()));

                        const double extraMassLB =
                            storedExcess != nullptr ? double(storedExcess->get(cell.x(), cell.y(), cell.z()))
                                                    : 0.0;

                        const double liquidMassLB = fraction * density + extraMassLB;

                        // Momentum density of the underlying LB liquid.
                        const auto momentumDensity = pdfs->getMomentumDensity(cell.x(), cell.y(), cell.z());

                        // Only the liquid-filled part contributes liquid momentum.
                        const Eigen::Vector3d liquidMomentumLB(fraction * double(momentumDensity[0]),
                                                               fraction * double(momentumDensity[1]),
                                                               fraction * double(momentumDensity[2]));

                        if (!std::isfinite(liquidMassLB) || !liquidMomentumLB.allFinite()) {
                            throw std::runtime_error("Non-finite liquid state in a covered cell");
                        }

                        changes.coveredLiquidCells.push_back({&block, cell, liquidMassLB, liquidMomentumLB});

                        changes.coveredLiquidVolumeLB += fraction;
                        changes.coveredLiquidMassLB += liquidMassLB;
                        changes.coveredLiquidMomentumLB += liquidMomentumLB;
                    }
                } else if (wasParticle && !becomesParticle) {
                    ++changes.uncovered;
                    changes.uncoveredCells.push_back({&block, cell});
                } else if (wasParticle && becomesParticle) {
                    ++changes.retained;
                }
            }
        }

        return changes;
    };

    struct RemapResult {
        Eigen::Vector3d linearImpulseLB;
        Eigen::Vector3d angularImpulseAboutOldCenterLB;
        Eigen::Vector3d oldCenterLB;
    };

    auto remapParticle = [&](std::size_t movingParticle,
                             const Eigen::Vector3d &shiftedCenter) -> RemapResult {
        if (movingParticle >= particles.size()) {
            throw std::out_of_range("Invalid moving particle index");
        }

        if (!shiftedCenter.allFinite()) {
            throw std::runtime_error("Non-finite proposed particle center");
        }

        const auto shifted = detectParticleCellChanges(movingParticle, shiftedCenter);

        struct InterfaceRecipient {
            IBlock *block;
            Cell cell;
            Eigen::Vector3d positionLB;
        };

        std::vector<InterfaceRecipient> interfaceRecipients;

        for (auto &block : *blockForest) {
            const auto *flags = block.getData<FlagField_T>(flagFieldID);

            for (const auto &cell : flags->xyzSize()) {
                const auto cellFlag = flags->get(cell.x(), cell.y(), cell.z());

                if (!flagInfo.isInterface(cellFlag))
                    continue;

                const auto p = blockForest->getBlockLocalCellCenter(block, cell);

                const Eigen::Vector3d positionLB{static_cast<double>(p[0]), static_cast<double>(p[1]),
                                                 static_cast<double>(p[2])};

                // A recipient must remain outside the proposed particle.
                if (isSolidAfterMove(positionLB, movingParticle, shiftedCenter)) {
                    continue;
                }

                interfaceRecipients.push_back({&block, cell, positionLB});
            }
        }

        auto findInterfaceThroughLiquid = [&](const auto &source,
                                              double minimumDonorMassLB = 0.0) -> std::vector<std::size_t> {
            std::vector<std::size_t> found;

            double reachedInterfaceMassLB = 0.0;

            const auto *flags = source.block->template getData<FlagField_T>(flagFieldID);

            const auto *fill = source.block->template getData<ScalarField_T>(fillFieldID);

            const auto *pdfs = source.block->template getData<PdfField_T>(pdfFieldID);

            const auto ownedCells = flags->xyzSize();

            using CellKey = std::array<cell_idx_t, 3>;

            std::set<CellKey> visited;
            std::queue<Cell> frontier;

            visited.insert({source.cell.x(), source.cell.y(), source.cell.z()});

            // Start at the covered cell, then enter surviving fluid.
            frontier.push(source.cell);

            while (!frontier.empty()) {
                const std::size_t layerSize = frontier.size();

                for (std::size_t n = 0; n < layerSize; ++n) {
                    const Cell current = frontier.front();
                    frontier.pop();

                    for (std::size_t j = 0; j < interfaceRecipients.size(); ++j) {
                        const auto &recipient = interfaceRecipients[j];

                        if (recipient.block == source.block && recipient.cell == current) {
                            found.push_back(j);

                            const double donorDensity =
                                double(pdfs->getDensity(current.x(), current.y(), current.z()));
                            const double donorMass =
                                double(fill->get(current.x(), current.y(), current.z())) * donorDensity;
                            if (!std::isfinite(donorMass) || !std::isfinite(donorDensity) ||
                                donorDensity <= 0.0)
                                throw std::runtime_error("Invalid donor state during interface search");
                            reachedInterfaceMassLB += std::max(0.0, donorMass);

                            break;
                        }
                    }

                    for (auto d = LatticeModelStencil_T::beginNoCenter(); d != LatticeModelStencil_T::end();
                         ++d) {
                        const Cell next(current.x() + d.cx(), current.y() + d.cy(), current.z() + d.cz());

                        if (!ownedCells.contains(next))
                            continue;

                        const CellKey key{next.x(), next.y(), next.z()};

                        if (!visited.insert(key).second)
                            continue;

                        const auto flag = flags->get(next.x(), next.y(), next.z());

                        if (!flagInfo.isLiquid(flag) && !flagInfo.isInterface(flag)) {
                            continue;
                        }

                        const auto p = blockForest->getBlockLocalCellCenter(*source.block, next);

                        const Eigen::Vector3d positionLB{static_cast<double>(p[0]), static_cast<double>(p[1]),
                                                         static_cast<double>(p[2])};

                        // Never route through the proposed solid.
                        if (isSolidAfterMove(positionLB, movingParticle, shiftedCenter)) {
                            continue;
                        }

                        frontier.push(next);
                    }
                }

                // Finish the whole layer to collect equally near recipients.
                if (!found.empty() && reachedInterfaceMassLB >= minimumDonorMassLB) {
                    return found;
                }
            }

            return found;
        };

        struct CoveredCellRoute {
            std::size_t sourceIndex;
            std::vector<std::size_t> recipientIndices;
        };

        std::vector<CoveredCellRoute> coveredCellRoutes;

        double locallyRoutableMassLB = 0.0;
        double unroutedMassLB = 0.0;

        for (std::size_t i = 0; i < shifted.coveredLiquidCells.size(); ++i) {
            const auto &source = shifted.coveredLiquidCells[i];

            const auto p = blockForest->getBlockLocalCellCenter(*source.block, source.cell);

            const Eigen::Vector3d sourcePositionLB{static_cast<double>(p[0]), static_cast<double>(p[1]),
                                                   static_cast<double>(p[2])};

            CoveredCellRoute route;
            route.sourceIndex = i;

            for (std::size_t j = 0; j < interfaceRecipients.size(); ++j) {
                const double distanceSquared =
                    (interfaceRecipients[j].positionLB - sourcePositionLB).squaredNorm();

                // D3Q19 neighbors: face or edge neighbors.
                if (distanceSquared > 0.5 && distanceSquared <= 2.0 + 1e-12) {
                    route.recipientIndices.push_back(j);
                }
            }

            if (route.recipientIndices.empty()) {
                route.recipientIndices = findInterfaceThroughLiquid(source);
            }

            if (route.recipientIndices.empty())
                unroutedMassLB += source.massLB;
            else
                locallyRoutableMassLB += source.massLB;

            coveredCellRoutes.push_back(std::move(route));
        }

        struct LiquidTransfer {
            double massLB = 0.0;
            Eigen::Vector3d momentumLB = Eigen::Vector3d::Zero();
        };

        std::vector<LiquidTransfer> recipientTransfers(interfaceRecipients.size());

        LiquidTransfer pendingTransfer;

        for (const auto &route : coveredCellRoutes) {
            const auto &source = shifted.coveredLiquidCells.at(route.sourceIndex);

            if (route.recipientIndices.empty()) {
                // Keep an explicit record of everything not yet routed.
                pendingTransfer.massLB += source.massLB;
                pendingTransfer.momentumLB += source.momentumLB;
                continue;
            }

            const double weight = 1.0 / double(route.recipientIndices.size());

            for (const std::size_t index : route.recipientIndices) {
                auto &transfer = recipientTransfers.at(index);

                transfer.massLB += weight * source.massLB;
                transfer.momentumLB += weight * source.momentumLB;
            }
        }

        LiquidTransfer allocatedTransfer;

        for (const auto &transfer : recipientTransfers) {
            allocatedTransfer.massLB += transfer.massLB;
            allocatedTransfer.momentumLB += transfer.momentumLB;
        }

        const double massBalanceErrorLB =
            allocatedTransfer.massLB + pendingTransfer.massLB - shifted.coveredLiquidMassLB;

        const Eigen::Vector3d momentumBalanceErrorLB =
            allocatedTransfer.momentumLB + pendingTransfer.momentumLB - shifted.coveredLiquidMomentumLB;

        if (std::abs(massBalanceErrorLB) > 1e-10 * std::max(1.0, shifted.coveredLiquidMassLB) ||
            momentumBalanceErrorLB.norm() > 1e-10 * std::max(1.0, shifted.coveredLiquidMomentumLB.norm()) ||
            std::abs(locallyRoutableMassLB + unroutedMassLB - shifted.coveredLiquidMassLB) >
                1e-10 * std::max(1.0, shifted.coveredLiquidMassLB))
            throw std::runtime_error("Covered-cell routing failed its conservation check");

        struct UncoveredCellState {
            IBlock *block;
            Cell cell;

            double fill;
            double density;
            Eigen::Vector3d velocityLB;
        };

        std::vector<UncoveredCellState> uncoveredStates;

        double requestedRefillMassLB = 0.0;
        uint_t unsupportedUncoveredCells = 0;

        for (const auto &entry : shifted.uncoveredCells) {
            const auto *flags = entry.block->getData<FlagField_T>(flagFieldID);

            const auto *fill = entry.block->getData<ScalarField_T>(fillFieldID);

            const auto *pdfs = entry.block->getData<PdfField_T>(pdfFieldID);

            double fractionSum = 0.0;
            double liquidMassSum = 0.0;

            Eigen::Vector3d liquidMomentumSum = Eigen::Vector3d::Zero();

            uint_t neighborCount = 0;

            for (auto d = LatticeModelStencil_T::beginNoCenter(); d != LatticeModelStencil_T::end(); ++d) {
                const Cell neighbor(entry.cell.x() + d.cx(), entry.cell.y() + d.cy(),
                                    entry.cell.z() + d.cz());

                // Single-block development implementation.
                if (!flags->xyzSize().contains(neighbor))
                    continue;

                const auto flag = flags->get(neighbor.x(), neighbor.y(), neighbor.z());

                const bool isFluid = flagInfo.isLiquid(flag) || flagInfo.isInterface(flag);

                if (!isFluid && !flagInfo.isGas(flag))
                    continue;

                const auto p = blockForest->getBlockLocalCellCenter(*entry.block, neighbor);

                const Eigen::Vector3d positionLB{static_cast<double>(p[0]), static_cast<double>(p[1]),
                                                 static_cast<double>(p[2])};

                // Only use neighbors surviving the proposed move.
                if (isSolidAfterMove(positionLB, movingParticle, shiftedCenter)) {
                    continue;
                }

                ++neighborCount;

                // Gas contributes zero liquid fraction.
                if (!isFluid)
                    continue;

                const double fraction = double(fill->get(neighbor.x(), neighbor.y(), neighbor.z()));

                const double density = double(pdfs->getDensity(neighbor.x(), neighbor.y(), neighbor.z()));

                const auto momentum = pdfs->getMomentumDensity(neighbor.x(), neighbor.y(), neighbor.z());

                fractionSum += fraction;
                liquidMassSum += fraction * density;

                liquidMomentumSum += fraction * Eigen::Vector3d{static_cast<double>(momentum[0]),
                                                                static_cast<double>(momentum[1]),
                                                                static_cast<double>(momentum[2])};
            }

            if (neighborCount == 0) {
                ++unsupportedUncoveredCells;
                continue;
            }

            const double proposedFill = std::clamp(fractionSum / double(neighborCount), 0.0, 1.0);

            double proposedDensity = 1.0;
            Eigen::Vector3d proposedVelocity = Eigen::Vector3d::Zero();

            if (fractionSum > 1e-12 && liquidMassSum > 1e-12) {
                proposedDensity = liquidMassSum / fractionSum;
                proposedVelocity = liquidMomentumSum / liquidMassSum;
            }

            uncoveredStates.push_back(
                {entry.block, entry.cell, proposedFill, proposedDensity, proposedVelocity});

            requestedRefillMassLB += proposedFill * proposedDensity;
        }

        // Use interface cells reached by the displacement-routing search.
        std::vector<bool> eligibleRecipient(interfaceRecipients.size(), false);

        for (const auto &route : coveredCellRoutes) {
            for (const std::size_t index : route.recipientIndices)
                eligibleRecipient.at(index) = true;
        }

        // Newly uncovered wet cells also need access to liquid donors.
        for (const auto &state : uncoveredStates) {
            if (state.fill <= 0.0)
                continue;

            const double requiredDonorMassLB =
                std::max(0.0, requestedRefillMassLB - shifted.coveredLiquidMassLB);

            const auto donorIndices = findInterfaceThroughLiquid(state, requiredDonorMassLB);

            for (const std::size_t index : donorIndices)
                eligibleRecipient.at(index) = true;
        }

        std::vector<double> signedInterfaceMassLB(interfaceRecipients.size(), 0.0);
        for (std::size_t j = 0; j < interfaceRecipients.size(); ++j) {
            const auto &recipient = interfaceRecipients[j];
            const auto *fill = recipient.block->getData<ScalarField_T>(fillFieldID);
            const auto *pdfs = recipient.block->getData<PdfField_T>(pdfFieldID);
            const auto &cell = recipient.cell;
            const double density = double(pdfs->getDensity(cell.x(), cell.y(), cell.z()));
            const double mass = double(fill->get(cell.x(), cell.y(), cell.z())) * density;
            if (!std::isfinite(density) || density <= 0.0 || !std::isfinite(mass))
                throw std::runtime_error("Non-finite interface mass or invalid density in remap plan");
            signedInterfaceMassLB[j] = mass;
        }
        const double remainingMassLB = shifted.coveredLiquidMassLB - requestedRefillMassLB;
        // A deficit can also come from covering a signed/negative interface
        // cell, with no newly uncovered wet cell to seed the search above.
        // Expand the SAME connected-fluid search from every changed wet cell
        // until the union of eligible positive donors can cover the NET plan.
        auto selectedDonorMass = [&]() {
            double mass = 0.0;
            for (std::size_t j = 0; j < signedInterfaceMassLB.size(); ++j)
                if (eligibleRecipient[j]) mass += std::max(0.0, signedInterfaceMassLB[j]);
            return mass;
        };
        const double donorRequirementLB = std::max(0.0, -remainingMassLB);
        const double originalDonorMassLB = selectedDonorMass();
        if (originalDonorMassLB < donorRequirementLB) {
            auto extendDonors = [&](const auto& source) {
                if (selectedDonorMass() >= donorRequirementLB) return;
                const auto indices = findInterfaceThroughLiquid(source, donorRequirementLB);
                for (const auto j : indices) eligibleRecipient.at(j) = true;
            };
            for (const auto& source : shifted.coveredLiquidCells) extendDonors(source);
            for (const auto& source : uncoveredStates)
                if (source.fill > 0.0) extendDonors(source);
            WALBERLA_LOG_DEVEL_ON_ROOT("Expanded remap donors: particle=" << movingParticle
                << " required=" << donorRequirementLB << " initial=" << originalDonorMassLB
                << " expanded=" << selectedDonorMass()
                << " coveredMass=" << shifted.coveredLiquidMassLB
                << " refillMass=" << requestedRefillMassLB);
        }
        const auto interfaceMassPlan = capillary_coupling::makeInterfaceMassPlan(
            signedInterfaceMassLB, eligibleRecipient, remainingMassLB);
        const auto &recipientMassDeltaLB = interfaceMassPlan.delta;
        bool massPlanSupported = unsupportedUncoveredCells == 0 && pendingTransfer.massLB == 0.0;

        double summedRecipientDeltaLB = 0.0;

        for (const double delta : recipientMassDeltaLB)
            summedRecipientDeltaLB += delta;

        if (std::abs(summedRecipientDeltaLB + requestedRefillMassLB - shifted.coveredLiquidMassLB) >
            1e-10 * std::max(1.0, shifted.coveredLiquidMassLB))
            throw std::runtime_error("Combined remap mass plan does not balance");

        Eigen::Vector3d fluidMomentumChangeLB = Eigen::Vector3d::Zero();

        Eigen::Vector3d fluidAngularMomentumChangeLB = Eigen::Vector3d::Zero();

        // Use one fixed reference point for every angular-momentum term.
        const Eigen::Vector3d angularReferenceLB = particles.at(movingParticle).center;

        auto accumulateMomentumChange = [&](const Eigen::Vector3d &positionLB,
                                            const Eigen::Vector3d &momentumChangeLB) {
            fluidMomentumChangeLB += momentumChangeLB;

            fluidAngularMomentumChangeLB += (positionLB - angularReferenceLB).cross(momentumChangeLB);
        };

        // Covered cells lose their liquid and its momentum.
        for (const auto &source : shifted.coveredLiquidCells) {
            accumulateMomentumChange(cellCenterEigen(*source.block, source.cell), -source.momentumLB);
        }

        // Newly uncovered cells acquire their proposed liquid momentum.
        for (const auto &state : uncoveredStates) {
            const double massLB = state.fill * state.density;

            accumulateMomentumChange(cellCenterEigen(*state.block, state.cell), massLB * state.velocityLB);
        }

        // Surviving interface cells gain or lose liquid at their
        // existing velocity. Their PDFs are not changed by this plan.
        for (std::size_t j = 0; j < interfaceRecipients.size(); ++j) {
            const double massDeltaLB = recipientMassDeltaLB[j];

            if (massDeltaLB == 0.0)
                continue;

            const auto &recipient = interfaceRecipients[j];

            const auto *pdfs = recipient.block->getData<PdfField_T>(pdfFieldID);

            const auto &cell = recipient.cell;

            const double density = double(pdfs->getDensity(cell.x(), cell.y(), cell.z()));

            if (!std::isfinite(density) || density <= 0.0) {
                throw std::runtime_error("Invalid recipient density in momentum plan");
            }

            const auto momentumDensity = pdfs->getMomentumDensity(cell.x(), cell.y(), cell.z());

            const Eigen::Vector3d velocityLB{double(momentumDensity[0]) / density,
                                             double(momentumDensity[1]) / density,
                                             double(momentumDensity[2]) / density};

            accumulateMomentumChange(recipient.positionLB, massDeltaLB * velocityLB);
        }

        const Eigen::Vector3d particleRemapImpulseLB = -fluidMomentumChangeLB;

        const Eigen::Vector3d particleRemapAngularImpulseLB = -fluidAngularMomentumChangeLB;

        struct PlannedFluidCell {
            IBlock *block;
            Cell cell;

            double fill;
            double density;
            Eigen::Vector3d velocityLB;

            bool initializePdfs;
        };

        std::vector<PlannedFluidCell> plannedFluidCells;

        // Newly uncovered cells need initialized PDFs.
        for (const auto &state : uncoveredStates) {
            plannedFluidCells.push_back(
                {state.block, state.cell, state.fill, state.density, state.velocityLB, true});
        }

        // Existing interface cells retain their PDFs.
        // Their liquid mass changes through their fill level.
        for (std::size_t j = 0; j < interfaceRecipients.size(); ++j) {
            const double massDeltaLB = recipientMassDeltaLB[j];

            if (massDeltaLB == 0.0)
                continue;

            const auto &recipient = interfaceRecipients[j];
            const auto &cell = recipient.cell;

            const auto *fill = recipient.block->getData<ScalarField_T>(fillFieldID);

            const auto *pdfs = recipient.block->getData<PdfField_T>(pdfFieldID);

            const double oldFill = double(fill->get(cell.x(), cell.y(), cell.z()));

            const double density = double(pdfs->getDensity(cell.x(), cell.y(), cell.z()));

            if (!std::isfinite(density) || density <= 0.0) {
                throw std::runtime_error("Invalid density while planning fluid cells");
            }

            const auto velocity = pdfs->getVelocity(cell.x(), cell.y(), cell.z());

            const Eigen::Vector3d velocityLB{static_cast<double>(velocity[0]),
                                             static_cast<double>(velocity[1]),
                                             static_cast<double>(velocity[2])};

            const double proposedFill = oldFill + massDeltaLB / density;

            plannedFluidCells.push_back({recipient.block, cell, proposedFill, density, velocityLB, false});
        }

        uint_t invalidPlannedCells = 0;

        for (const auto &state : plannedFluidCells) {

            const bool valid = capillary_coupling::validRemapFill(state.fill, state.initializePdfs) &&
                               std::isfinite(state.density) && state.density > 0.0 &&
                               state.velocityLB.allFinite();

            if (!valid) {
                ++invalidPlannedCells;
                continue;
            }
        }

        massPlanSupported = massPlanSupported && invalidPlannedCells == 0;

        for (const auto &uncovered : uncoveredStates) {
            // Positive-fill uncovered cells will become interface cells.
            // Exactly zero-fill uncovered cells will become gas.
            if (uncovered.fill != 0.0)
                continue;

            auto *block = uncovered.block;

            const auto *flags = block->getData<FlagField_T>(flagFieldID);

            const auto *fill = block->getData<ScalarField_T>(fillFieldID);

            const auto *pdfs = block->getData<PdfField_T>(pdfFieldID);

            for (auto d = LatticeModelStencil_T::beginNoCenter(); d != LatticeModelStencil_T::end(); ++d) {
                const Cell neighbor(uncovered.cell.x() + d.cx(), uncovered.cell.y() + d.cy(),
                                    uncovered.cell.z() + d.cz());

                if (!flags->xyzSize().contains(neighbor)) {
                    throw std::runtime_error("Topology planning currently requires "
                                             "changed cells away from block boundaries");
                }

                const auto flag = flags->get(neighbor.x(), neighbor.y(), neighbor.z());

                if (!flagInfo.isLiquid(flag))
                    continue;

                const Eigen::Vector3d positionLB = cellCenterEigen(*block, neighbor);

                // A cell covered by the new sphere becomes solid instead.
                if (isSolidAfterMove(positionLB, movingParticle, shiftedCenter)) {
                    continue;
                }

                // Avoid duplicates. Already-planned fluid cells use
                // their planned phase rather than their original flag.
                const bool alreadyPlanned = std::any_of(
                    plannedFluidCells.begin(), plannedFluidCells.end(), [&](const PlannedFluidCell &state) {
                        return state.block == block && state.cell == neighbor;
                    });

                if (alreadyPlanned)
                    continue;

                const double oldFill = double(fill->get(neighbor.x(), neighbor.y(), neighbor.z()));

                const double density = double(pdfs->getDensity(neighbor.x(), neighbor.y(), neighbor.z()));

                const auto velocity = pdfs->getVelocity(neighbor.x(), neighbor.y(), neighbor.z());

                const Eigen::Vector3d velocityLB{static_cast<double>(velocity[0]),
                                                 static_cast<double>(velocity[1]),
                                                 static_cast<double>(velocity[2])};

                plannedFluidCells.push_back({block, neighbor, oldFill, density, velocityLB, false});
            }
        }

        uint_t remapBlockCount = 0;

        for (auto it = blockForest->begin(); it != blockForest->end(); ++it) {
            ++remapBlockCount;
        }

        if (numProcesses != uint_c(1) || remapBlockCount != uint_c(1)) {
            throw std::runtime_error("The particle remapper currently requires "
                                     "one process and one block");
        }

        if (!massPlanSupported) {
            std::ostringstream message;
            message << "Unsupported particle remap: uncoveredWithoutNeighbors=" << unsupportedUncoveredCells
                    << " pendingMassLB=" << pendingTransfer.massLB
                    << " coveredMassLB=" << shifted.coveredLiquidMassLB
                    << " refillMassLB=" << requestedRefillMassLB
                    << " availableDonorMassLB=" << interfaceMassPlan.availablePositiveMass
                    << " eligibleRecipients=" << interfaceMassPlan.eligibleCount
                    << " invalidPlannedCells=" << invalidPlannedCells;
            throw std::runtime_error(message.str());
        }

        struct StoredExcessTransfer {
            IBlock *block;
            Cell source;
            Cell destination;
            double massLB;
        };

        std::vector<StoredExcessTransfer> storedExcessTransfers;

        // Validate all states, including the later topology promotions.
        for (const auto &state : plannedFluidCells) {
            if (!capillary_coupling::validRemapFill(state.fill, state.initializePdfs) ||
                !std::isfinite(state.density) || state.density <= 0.0 || !state.velocityLB.allFinite()) {
                throw std::runtime_error("Invalid fluid state before remapping");
            }

            if (!state.initializePdfs && remapExcessMassFieldID != ConstBlockDataID()) {
                const auto *storedExcess = state.block->getData<ScalarField_T>(remapExcessMassFieldID);

                const double extra =
                    double(storedExcess->get(state.cell.x(), state.cell.y(), state.cell.z()));

                if (!std::isfinite(extra)) {
                    throw std::runtime_error("Non-finite stored excess mass");
                }

                if (extra != 0.0) {
                    const auto *flags = state.block->getData<FlagField_T>(flagFieldID);

                    bool destinationFound = false;

                    for (auto d = LatticeModelStencil_T::beginNoCenter(); d != LatticeModelStencil_T::end();
                         ++d) {
                        const Cell neighbor(state.cell.x() + d.cx(), state.cell.y() + d.cy(),
                                            state.cell.z() + d.cz());

                        if (!flags->xyzSize().contains(neighbor))
                            continue;

                        const auto flag = flags->get(neighbor.x(), neighbor.y(), neighbor.z());

                        if (!flagInfo.isLiquid(flag))
                            continue;

                        const Eigen::Vector3d positionLB = cellCenterEigen(*state.block, neighbor);

                        // The destination must remain outside the particle.
                        if (isSolidAfterMove(positionLB, movingParticle, shiftedCenter)) {
                            continue;
                        }

                        // It must also remain an unchanged liquid cell.
                        const bool destinationChanges = std::any_of(
                            plannedFluidCells.begin(), plannedFluidCells.end(),
                            [&](const PlannedFluidCell &candidate) {
                                return candidate.block == state.block && candidate.cell == neighbor;
                            });

                        if (destinationChanges)
                            continue;

                        storedExcessTransfers.push_back({state.block, state.cell, neighbor, extra});

                        destinationFound = true;
                        break;
                    }

                    if (!destinationFound) {
                        throw std::runtime_error("No unchanged neighboring liquid cell can "
                                                 "receive stored excess mass");
                    }
                }
            }
        }

        auto measureRemapLiquidState = [&]() {
            double massLB = 0.0;
            Eigen::Vector3d momentumLB = Eigen::Vector3d::Zero();

            for (auto &block : *blockForest) {
                const auto *flags = block.getData<FlagField_T>(flagFieldID);

                const auto *fill = block.getData<ScalarField_T>(fillFieldID);

                const auto *pdfs = block.getData<PdfField_T>(pdfFieldID);

                const ScalarField_T *storedExcess = nullptr;

                if (remapExcessMassFieldID != ConstBlockDataID()) {
                    storedExcess = block.getData<ScalarField_T>(remapExcessMassFieldID);
                }

                for (const auto &cell : flags->xyzSize()) {
                    const auto flag = flags->get(cell.x(), cell.y(), cell.z());

                    if (!flagInfo.isLiquid(flag) && !flagInfo.isInterface(flag)) {
                        continue;
                    }

                    const double fraction = double(fill->get(cell.x(), cell.y(), cell.z()));

                    massLB += fraction * double(pdfs->getDensity(cell.x(), cell.y(), cell.z()));

                    if (storedExcess != nullptr) {
                        massLB += double(storedExcess->get(cell.x(), cell.y(), cell.z()));
                    }

                    const auto p = pdfs->getMomentumDensity(cell.x(), cell.y(), cell.z());

                    momentumLB +=
                        fraction * Eigen::Vector3d{static_cast<double>(p[0]), static_cast<double>(p[1]),
                                                   static_cast<double>(p[2])};
                }
            }

            return std::make_pair(massLB, momentumLB);
        };

        const auto beforeRemap = measureRemapLiquidState();

        // Newly covered cells become particle boundary cells.
        for (const auto &entry : shifted.coveredCells) {
            auto *handling =
                entry.block->getData<BoundaryHandling_T>(freeSurfaceBoundaryHandling->getHandlingID());

            auto *flags = handling->getFlagField();

            auto *fill = entry.block->getData<ScalarField_T>(fillFieldID);

            const auto &cell = entry.cell;

            if (remapExcessMassFieldID != ConstBlockDataID()) {
                const BlockDataID writableExcessID{static_cast<uint_t>(remapExcessMassFieldID)};

                auto *storedExcess = entry.block->getData<ScalarField_T>(writableExcessID);

                storedExcess->get(cell.x(), cell.y(), cell.z()) = real_c(0);
            }

            // Gas is not part of BoundaryHandling's fluid-domain mask,
            // so remove its flag explicitly.
            if (flagInfo.isGas(flags->get(cell.x(), cell.y(), cell.z()))) {
                handling->removeFlag(flagInfo.gasFlag, cell.x(), cell.y(), cell.z());
            }

            handling->forceBoundary(FreeSurfaceBoundaryHandling_T::ubbFlagID, cell.x(), cell.y(), cell.z(),
                                    UBB_T::Velocity(real_c(0), real_c(0), real_c(0)));

            // Match waLBerla's obstacle-fill initialization.
            // Solid cells are excluded from liquid-mass measurements.
            fill->get(cell.x(), cell.y(), cell.z()) = real_c(1);
        }

        // Apply uncovered-cell states, recipient changes,
        // and liquid-to-interface promotions.
        // Decide topology from the COMPLETE proposed state, not commit order.
        auto proposedCellIsGas = [&](IBlock *block, const Cell &cell) {
            if (isSolidAfterMove(cellCenterEigen(*block, cell), movingParticle, shiftedCenter))
                return false;
            const auto planned = std::find_if(plannedFluidCells.begin(), plannedFluidCells.end(),
                                              [&](const PlannedFluidCell &candidate) {
                                                  return candidate.block == block && candidate.cell == cell;
                                              });
            if (planned != plannedFluidCells.end())
                return planned->initializePdfs && planned->fill == 0.0;
            const auto *flags = block->getData<FlagField_T>(flagFieldID);
            return flagInfo.isGas(flags->get(cell.x(), cell.y(), cell.z()));
        };

        for (const auto &state : plannedFluidCells) {
            auto *handling =
                state.block->getData<BoundaryHandling_T>(freeSurfaceBoundaryHandling->getHandlingID());

            auto *fill = state.block->getData<ScalarField_T>(fillFieldID);

            auto *pdfs = state.block->getData<PdfField_T>(pdfFieldID);

            const auto &cell = state.cell;

            handling->clear(cell.x(), cell.y(), cell.z());

            const bool becomesGas = state.initializePdfs && state.fill == 0.0;

            bool becomesLiquid = state.initializePdfs && state.fill == 1.0;
            if (becomesLiquid) {
                for (auto d = LatticeModelStencil_T::beginNoCenter(); d != LatticeModelStencil_T::end();
                     ++d) {
                    const Cell neighbor(cell.x() + d.cx(), cell.y() + d.cy(), cell.z() + d.cz());
                    if (proposedCellIsGas(state.block, neighbor)) {
                        becomesLiquid = false;
                        break;
                    }
                }
            }

            if (becomesGas) {
                handling->setFlag(flagInfo.gasFlag, cell.x(), cell.y(), cell.z());
            } else {
                handling->setDomain(becomesLiquid ? flagInfo.liquidFlag : flagInfo.interfaceFlag, cell.x(),
                                    cell.y(), cell.z());
            }

            fill->get(cell.x(), cell.y(), cell.z()) = real_c(state.fill);

            if (state.initializePdfs) {

                if (remapExcessMassFieldID != ConstBlockDataID()) {
                    const BlockDataID writableExcessID{static_cast<uint_t>(remapExcessMassFieldID)};

                    auto *storedExcess = state.block->getData<ScalarField_T>(writableExcessID);

                    storedExcess->get(cell.x(), cell.y(), cell.z()) = real_c(0);
                }

                pdfs->setToEquilibrium(cell,
                                       Vector3<real_t>(real_c(state.velocityLB.x()),
                                                       real_c(state.velocityLB.y()),
                                                       real_c(state.velocityLB.z())),
                                       real_c(state.density));
            }
        }

        for (auto &block : *blockForest) {
            auto *handling = block.getData<BoundaryHandling_T>(freeSurfaceBoundaryHandling->getHandlingID());

            handling->refresh();
        }

        if (!storedExcessTransfers.empty()) {
            const BlockDataID writableExcessID{static_cast<uint_t>(remapExcessMassFieldID)};

            for (const auto &transfer : storedExcessTransfers) {
                auto *storedExcess = transfer.block->getData<ScalarField_T>(writableExcessID);

                storedExcess->get(transfer.destination.x(), transfer.destination.y(),
                                  transfer.destination.z()) += real_c(transfer.massLB);

                storedExcess->get(transfer.source.x(), transfer.source.y(), transfer.source.z()) = real_c(0);
            }
        }

        const auto afterRemap = measureRemapLiquidState();

        const Eigen::Vector3d momentumPredictionErrorLB =
            afterRemap.second - beforeRemap.second - fluidMomentumChangeLB;

        const auto geometryAfterRemap = detectParticleCellChanges(movingParticle, shiftedCenter);

        const double remapMassTolerance = 1e-10 * std::max(1.0, std::abs(beforeRemap.first));
        if (std::abs(afterRemap.first - beforeRemap.first) > remapMassTolerance ||
            momentumPredictionErrorLB.norm() > 1e-10 * std::max(1.0, beforeRemap.second.norm()) ||
            geometryAfterRemap.covered != 0 || geometryAfterRemap.uncovered != 0)
            throw std::runtime_error("Particle remap failed its mass, momentum, or geometry check");

        // Synchronize all geometric representations with the remapped sphere.
        auto &particle = particles.at(movingParticle);

        particle.center = shiftedCenter;
        particle.surface = makeSphereSurface(particle.center, double(solidRadius), particle.surfaceLevel);

        // Preserve the existing particle-0 aliases during the refactor.
        if (movingParticle == 0) {
            solidCenterEigen = particle.center;

            solidCenter = Vector3<real_t>(real_c(particle.center.x()), real_c(particle.center.y()),
                                          real_c(particle.center.z()));

            solidSurface = particle.surface;
        }

        return RemapResult{particleRemapImpulseLB, particleRemapAngularImpulseLB, angularReferenceLB};
    }; // End remapParticle.

    // communication after initialization
    Communication_T communication(blockForest, flagFieldID, fillFieldID, forceDensityFieldID);
    communication();

    PdfCommunication_T pdfCommunication(blockForest, pdfFieldID);
    pdfCommunication();

    // add bubble model
    std::shared_ptr<bubble_model::BubbleModelBase> bubbleModel = nullptr;
    if (enableBubbleModel) {
        const std::shared_ptr<bubble_model::BubbleModel<LatticeModelStencil_T>> bubbleModelDerived =
            std::make_shared<bubble_model::BubbleModel<LatticeModelStencil_T>>(blockForest,
                                                                               enableBubbleSplits);
        bubbleModelDerived->initFromFillLevelField(fillFieldID);
        bubbleModelDerived->setAtmosphere(
            Cell(domainSize[0] - uint_c(1), domainSize[1] - uint_c(1), domainSize[2] - uint_c(1)), real_c(1));

        bubbleModel = std::static_pointer_cast<bubble_model::BubbleModelBase>(bubbleModelDerived);
    } else {
        bubbleModel = std::make_shared<bubble_model::BubbleModelConstantPressure>(real_c(1));
    }

    // set density in non-liquid or non-interface cells to 1 (after initializing
    // with hydrostatic pressure)
    setDensityInNonFluidCellsToOne<FlagField_T, PdfField_T>(blockForest, flagInfo, flagFieldID, pdfFieldID);

    // create timeloop
    SweepTimeloop timeloop(blockForest, timesteps);

    Vector3<real_t> sphereForceLB(real_c(0), real_c(0), real_c(0));

    Eigen::Vector3d sphereTorqueLB = Eigen::Vector3d::Zero();

    uint_t sphereLinkCount = uint_c(0);

    auto updateSphereWallVelocity = [&]() {
        std::vector<Eigen::Vector3d> velocitiesLB(particles.size(), Eigen::Vector3d::Zero());

        std::vector<Eigen::Vector3d> angularVelocitiesLB(particles.size(), Eigen::Vector3d::Zero());

        for (std::size_t i = 0; i < particles.size(); ++i) {
            if (advanceDEM) {
                velocitiesLB[i] = dem.particleVelocitySI(i) * (units.dtSI / units.dxSI);

                angularVelocitiesLB[i] = dem.particleAngularVelocitySI(i) * units.dtSI;
            }

            if (!velocitiesLB[i].allFinite() || !angularVelocitiesLB[i].allFinite()) {
                throw std::runtime_error("DEM supplied non-finite particle velocity");
            }
        }

        for (auto &block : *blockForest) {
            auto *handling = block.getData<BoundaryHandling_T>(freeSurfaceBoundaryHandling->getHandlingID());

            auto *flags = handling->getFlagField();

            const auto particleMask = flags->getFlag(FreeSurfaceBoundaryHandling_T::ubbFlagID);

            auto &ubb =
                handling->template getBoundaryCondition<UBB_T>(FreeSurfaceBoundaryHandling_T::ubbBoundaryID);

            for (const auto &cell : flags->xyzSizeWithGhostLayer()) {
                if ((flags->get(cell.x(), cell.y(), cell.z()) & particleMask) == 0) {
                    continue;
                }

                const Eigen::Vector3d position = cellCenterEigen(block, cell);

                const std::size_t owner = particleOwnerAt(position);

                if (owner == particles.size()) {
                    throw std::runtime_error("Solid boundary cell has no particle owner");
                }

                const Eigen::Vector3d wallVelocity =
                    velocitiesLB[owner] +
                    angularVelocitiesLB[owner].cross(position - particles[owner].center);

                const UBB_T::Velocity configuration(real_c(wallVelocity.x()), real_c(wallVelocity.y()),
                                                    real_c(wallVelocity.z()));

                ubb.registerCell(particleMask, cell.x(), cell.y(), cell.z(), configuration);
            }
        }
    };

    double movingWallMassCorrectionLB = 0.0;

    std::vector<Eigen::Vector3d> particleWallForcesLB(particles.size(), Eigen::Vector3d::Zero());

    std::vector<Eigen::Vector3d> particleWallTorquesLB(particles.size(), Eigen::Vector3d::Zero());

    std::vector<uint_t> particleWallLinkCounts(particles.size(), uint_c(0));

    auto measureSphereForce = [&]() {
        sphereForceLB = Vector3<real_t>(real_c(0), real_c(0), real_c(0));
        sphereTorqueLB.setZero();
        sphereLinkCount = uint_c(0);
        movingWallMassCorrectionLB = 0.0;

        for (std::size_t i = 0; i < particles.size(); ++i) {
            particleWallForcesLB[i].setZero();
            particleWallTorquesLB[i].setZero();
            particleWallLinkCounts[i] = uint_c(0);
        }

        for (auto &block : *blockForest) {
            const auto *pdfs = block.getData<PdfField_T>(pdfFieldID);

            const auto *fill = block.getData<ScalarField_T>(fillFieldID);

            const auto *flags = block.getData<FlagField_T>(flagFieldID);

            if (remapExcessMassFieldID == ConstBlockDataID()) {
                throw std::runtime_error("Moving-wall mass correction requires "
                                         "the fallback-liquid excess-mass field");
            }

            const BlockDataID writableExcessID{static_cast<uint_t>(remapExcessMassFieldID)};

            auto *wallExcessMass = block.getData<ScalarField_T>(writableExcessID);

            auto *handling = block.getData<BoundaryHandling_T>(freeSurfaceBoundaryHandling->getHandlingID());

            const auto &ubb =
                handling->template getBoundaryCondition<UBB_T>(FreeSurfaceBoundaryHandling_T::ubbBoundaryID);

            const auto particleBoundaryMask = flags->getFlag(FreeSurfaceBoundaryHandling_T::ubbFlagID);

            // Visit owned fluid cells only, avoiding duplicate contributions
            // from ghost cells.
            for (const auto &cell : flags->xyzSize()) {
                const auto fluidFlag = flags->get(cell.x(), cell.y(), cell.z());

                if (!flagInfo.isLiquid(fluidFlag) && !flagInfo.isInterface(fluidFlag)) {
                    continue;
                }

                for (auto d = LatticeModelStencil_T::beginNoCenter(); d != LatticeModelStencil_T::end();
                     ++d) {
                    const Cell neighbor(cell.x() + d.cx(), cell.y() + d.cy(), cell.z() + d.cz());

                    const auto neighborFlag = flags->get(neighbor.x(), neighbor.y(), neighbor.z());

                    if ((neighborFlag & particleBoundaryMask) == 0)
                        continue;

                    // The floor also uses NoSlip. Count only solid cells
                    // belonging to our sphere.
                    const auto neighborPosition = blockForest->getBlockLocalCellCenter(block, neighbor);

                    const Eigen::Vector3d solidCellPosition{double(neighborPosition[0]),
                                                            double(neighborPosition[1]),
                                                            double(neighborPosition[2])};

                    const std::size_t owner = particleOwnerAt(solidCellPosition);

                    if (owner == particles.size()) {
                        throw std::runtime_error("Particle boundary link has no particle owner");
                    }

                    const real_t outgoingPDF = pdfs->get(cell.x(), cell.y(), cell.z(), d.toIdx());

                    const real_t density = pdfs->getDensity(cell.x(), cell.y(), cell.z());

                    // Read the exact velocity used by the boundary condition.
                    const auto &wallVelocity = ubb.getValue(neighbor.x(), neighbor.y(), neighbor.z());

                    const real_t directionDotVelocity = real_c(d.cx()) * wallVelocity[0] +
                                                        real_c(d.cy()) * wallVelocity[1] +
                                                        real_c(d.cz()) * wallVelocity[2];

                    const real_t weight = real_c(LatticeModel_T::w[d.toIdx()]);

                    // Same reflected PDF as waLBerla's compressible UBB.
                    const real_t reflectedPDF =
                        outgoingPDF - real_c(6) * density * weight * directionDotVelocity;

                    // Full-liquid cells acquire the moving-wall PDF mass change
                    // during streaming. Compensate it through excess-mass storage.
                    if (flagInfo.isLiquid(fluidFlag)) {
                        const real_t massCorrection = outgoingPDF - reflectedPDF;

                        wallExcessMass->get(cell.x(), cell.y(), cell.z()) += massCorrection;

                        movingWallMassCorrectionLB += double(massCorrection);
                    }

                    // Subtract the ambient reference pressure.
                    const real_t rhoAmbient = real_c(1);
                    const real_t momentumTransfer =
                        outgoingPDF + reflectedPDF - real_c(2) * weight * rhoAmbient;

                    const double fraction = double(fill->get(cell.x(), cell.y(), cell.z()));

                    const Eigen::Vector3d linkForce =
                        fraction * Eigen::Vector3d(double(momentumTransfer) * double(d.cx()),
                                                   double(momentumTransfer) * double(d.cy()),
                                                   double(momentumTransfer) * double(d.cz()));

                    particleWallForcesLB[owner] += linkForce;

                    const auto fluidPosition = blockForest->getBlockLocalCellCenter(block, cell);

                    const Eigen::Vector3d boundaryPosition(
                        0.5 * double(fluidPosition[0] + neighborPosition[0]),
                        0.5 * double(fluidPosition[1] + neighborPosition[1]),
                        0.5 * double(fluidPosition[2] + neighborPosition[2]));

                    particleWallTorquesLB[owner] +=
                        (boundaryPosition - particles[owner].center).cross(linkForce);

                    ++particleWallLinkCounts[owner];
                }
            }
        }

        // Temporary particle-0 aliases for code updated in the next step.
        sphereForceLB =
            Vector3<real_t>(real_c(particleWallForcesLB[0].x()), real_c(particleWallForcesLB[0].y()),
                            real_c(particleWallForcesLB[0].z()));

        sphereTorqueLB = particleWallTorquesLB[0];

        sphereLinkCount = uint_c(0);
        for (const auto count : particleWallLinkCounts) {
            sphereLinkCount += count;
        }
    };

    timeloop.addFuncBeforeTimeStep(
        [&]() {
            updateSphereWallVelocity();
            measureSphereForce();
        },
        "Apply DEM wall velocities and measure fluid loads");

    // Laplace pressure = 2 * surface tension * curvature; curvature computation
    // is not necessary with 0 surface tension
    bool computeCurvature = false;
    if (!realIsEqual(surfaceTension, real_c(0), real_c(1e-14))) {
        computeCurvature = true;
    }

    // add surface geometry handler
    const SurfaceGeometryHandler<LatticeModel_T, FlagField_T, ScalarField_T, VectorField_T> geometryHandler(
        blockForest, freeSurfaceBoundaryHandling, fillFieldID, curvatureModel, computeCurvature,
        enableWetting, contactAngle);

    double minimumResolvedFillLB = 1.0;
    double maximumResolvedFillLB = 0.0;
    auto measureLiquidState = [&]() {
        capillary_coupling::MomentumState state;
        minimumResolvedFillLB = 1.0;
        maximumResolvedFillLB = 0.0;
        for (auto &block : *blockForest) {
            const auto *flags = block.getData<FlagField_T>(flagFieldID);
            const auto *fill = block.getData<ScalarField_T>(fillFieldID);
            const auto *pdfs = block.getData<PdfField_T>(pdfFieldID);
            for (const auto &cell : flags->xyzSize()) {
                const auto flag = flags->get(cell.x(), cell.y(), cell.z());
                if (!flagInfo.isLiquid(flag) && !flagInfo.isInterface(flag))
                    continue;
                const double fraction = double(fill->get(cell.x(), cell.y(), cell.z()));
                // Native interface conversion has hysteresis; redistribution can
                // overshoot until a later conversion pass. Keep the signed value
                // in mass and momentum accounting rather than clipping it.
                if (!std::isfinite(fraction))
                    throw std::runtime_error("Non-finite resolved liquid fraction");
                minimumResolvedFillLB = std::min(minimumResolvedFillLB, fraction);
                maximumResolvedFillLB = std::max(maximumResolvedFillLB, fraction);
                const auto j = pdfs->getMomentumDensity(cell.x(), cell.y(), cell.z());
                const Eigen::Vector3d momentum =
                    fraction * Eigen::Vector3d{double(j[0]), double(j[1]), double(j[2])};
                const Eigen::Vector3d position = cellCenterEigen(block, cell);
                if (fraction > 1e-12 && (position.x() < 2.0 || position.y() < 2.0 || position.z() < 2.0 ||
                                         position.x() > double(domainSize[0]) - 2.0 ||
                                         position.y() > double(domainSize[1]) - 2.0 ||
                                         position.z() > double(domainSize[2]) - 2.0))
                    throw std::runtime_error("Isolated-particle projection requires "
                                             "liquid clear of all domain boundaries");
                if (!momentum.allFinite())
                    throw std::runtime_error("Non-finite resolved liquid momentum");
                state.linear += momentum;
                state.angular += (position - solidCenterEigen).cross(momentum);
            }
        }
        return state;
    };

    geometryHandler.addSweeps(timeloop);

    const ConstBlockDataID curvatureFieldID = geometryHandler.getConstCurvatureFieldID();
    const ConstBlockDataID normalFieldID = geometryHandler.getConstNormalFieldID();

    timeloop.addFuncBeforeTimeStep([&]() { interfacePlanes.clear(); }, "Reset interface-plane samples");

    timeloop.add() << Sweep(
        [&](IBlock *block) {
            using Vec = Eigen::Vector3d;

            const auto *flags = block->getData<FlagField_T>(flagFieldID);

            const auto *fill = block->getData<ScalarField_T>(fillFieldID);

            const auto *normals = block->getData<VectorField_T>(normalFieldID);

            for (const auto &cell : flags->xyzSize()) {
                const auto flag = flags->get(cell.x(), cell.y(), cell.z());

                if (!flagInfo.isInterface(flag))
                    continue;

                const auto cellCenter = blockForest->getBlockLocalCellCenter(*block, cell);

                const Vec center(cellCenter[0], cellCenter[1], cellCenter[2]);

                // Collect planes within three lattice spacings
                // of the analytic sphere surface.
                bool nearParticle = false;

                for (const auto &particle : particles) {
                    const double surfaceDistance =
                        std::abs((center - particle.center).norm() - double(solidRadius));

                    // Include cell centers whose reconstructed interface points
                    // can fall within the four-cell support radius.
                    if (surfaceDistance <= 4.5 + 0.5 * std::sqrt(3.0)) {
                        nearParticle = true;
                        break;
                    }
                }

                if (!nearParticle)
                    continue;

                const real_t fraction = fill->get(cell.x(), cell.y(), cell.z());

                if (!std::isfinite(fraction) || fraction <= real_c(1e-6) ||
                    fraction >= real_c(1) - real_c(1e-6)) {
                    continue;
                }

                auto normal = normals->get(cell.x(), cell.y(), cell.z());

                const Vec normalEigen(normal[0], normal[1], normal[2]);

                if (!normalEigen.allFinite() || normalEigen.squaredNorm() < 1e-14) {
                    continue;
                }

                normal = normal.getNormalizedOrZero();

                const auto localPoint = walberla::free_surface::getInterfacePoint(normal, fraction);

                // getInterfacePoint uses [0,1]^3 coordinates.
                // This demo has unit lattice spacing.
                const Vec globalPoint = center + Vec(localPoint[0] - real_c(0.5), localPoint[1] - real_c(0.5),
                                                     localPoint[2] - real_c(0.5));

                if (!globalPoint.allFinite())
                    continue;

                interfacePlanes.push_back({globalPoint, normalEigen.normalized()});
            }
        },
        "Collect interface planes near sphere");

    // add boundary handling for standard boundaries and free surface boundaries
    const SurfaceDynamicsHandler<LatticeModel_T, FlagField_T, ScalarField_T, VectorField_T> dynamicsHandler(
        blockForest, pdfFieldID, flagFieldID, fillFieldID, forceDensityFieldID, normalFieldID,
        curvatureFieldID, freeSurfaceBoundaryHandling, bubbleModel, pdfReconstructionModel, pdfRefillingModel,
        excessMassDistributionModel, relaxationRate, force, surfaceTension, useSimpleMassExchange,
        cellConversionThreshold, cellConversionForceThreshold);

    remapExcessMassFieldID = dynamicsHandler.getConstExcessMassFieldID();

    // Synchronize fields after the imposed startup displacement.
    communication();
    pdfCommunication();

    dynamicsHandler.addSweeps(timeloop);

    // All angular momenta below use the particle center held fixed during this LB
    // step.
    capillary_coupling::MomentumState stepBefore, rawResidual, correctedResidual;
    ContactLineLoad stepLoad;
    SurfaceLineResult stepLine;
    std::vector<ContactLineLoad> particleStepLoads(particles.size());

    std::vector<SurfaceLineResult> particleStepLines(particles.size());
    double maxAppliedCorrectionLB = 0.0;

    // Apply a specified fluid impulse, with angular impulse measured about
    // solidCenterEigen. Shared by the LB residual and moving-lever correction.
    auto applyFluidMomentumCorrection = [&](const Eigen::Vector3d &linearImpulseLB,
                                            const Eigen::Vector3d &angularImpulseLB) {
        double maxCorrection = 0.0;
        struct Recipient {
            PdfField_T *pdfs;
            Cell cell;
            Eigen::Vector3d position;
            double density;
        };
        std::vector<Recipient> recipients;
        std::vector<capillary_coupling::ProjectionCell> support;
        for (auto &block : *blockForest) {
            const auto *flags = block.getData<FlagField_T>(flagFieldID);
            const auto *fill = block.getData<ScalarField_T>(fillFieldID);
            auto *pdfs = block.getData<PdfField_T>(pdfFieldID);
            for (const auto &cell : flags->xyzSize()) {
                const auto flag = flags->get(cell.x(), cell.y(), cell.z());
                if (!flagInfo.isLiquid(flag) && !flagInfo.isInterface(flag))
                    continue;
                const double fraction = double(fill->get(cell.x(), cell.y(), cell.z()));
                if (!std::isfinite(fraction))
                    throw std::runtime_error("Non-finite liquid fraction before projection");
                // Correction support requires positive masses. Nonpositive
                // cells are unchanged, but their signed momentum remains
                // included in the measured residual.
                if (fraction <= 1e-12)
                    continue;
                const Eigen::Vector3d position = cellCenterEigen(block, cell);

                const double supportRadius = double(solidRadius) + 3.0;
                const double supportRadiusSquared = supportRadius * supportRadius;

                bool nearParticle = false;

                for (const auto &particle : particles) {
                    if ((position - particle.center).squaredNorm() <= supportRadiusSquared) {
                        nearParticle = true;
                        break;
                    }
                }

                if (!nearParticle) {
                    continue;
                }
                const double density = double(pdfs->getDensity(cell.x(), cell.y(), cell.z()));
                if (!std::isfinite(density) || density <= 0.0)
                    throw std::runtime_error("Invalid liquid density before projection");
                recipients.push_back({pdfs, cell, position, density});
                support.push_back({position, fraction * density});
            }
        }

        const auto projection = capillary_coupling::makeMomentumProjection(
            support, linearImpulseLB, angularImpulseLB, solidCenterEigen);

        // Preflight the entire correction before writing any PDFs.
        for (const auto &recipient : recipients) {
            const auto du = projection.velocityIncrement(recipient.position);
            maxCorrection = std::max(maxCorrection, du.norm());
            if (!du.allFinite() || du.norm() > maxCorrectionSpeedLB)
                throw std::runtime_error("Momentum correction exceeds maxCorrectionSpeedLB; reduce "
                                         "timestep or refine the grid");
            const auto oldVelocity =
                recipient.pdfs->getVelocity(recipient.cell.x(), recipient.cell.y(), recipient.cell.z());
            const Eigen::Vector3d newVelocity =
                du + Eigen::Vector3d{double(oldVelocity[0]), double(oldVelocity[1]), double(oldVelocity[2])};
            if (!newVelocity.allFinite() || newVelocity.norm() > 0.1)
                throw std::runtime_error("Corrected liquid speed exceeds 0.1 LB; reduce timestep");
            for (auto d = LatticeModelStencil_T::begin(); d != LatticeModelStencil_T::end(); ++d) {
                const Eigen::Vector3d direction{double(d.cx()), double(d.cy()), double(d.cz())};
                const double delta = capillary_coupling::momentumPdfIncrement(
                    double(LatticeModel_T::w[d.toIdx()]), recipient.density, direction, du);
                const double proposed = double(recipient.pdfs->get(recipient.cell.x(), recipient.cell.y(),
                                                                   recipient.cell.z(), d.toIdx())) +
                                        delta;
                if (!std::isfinite(proposed) || proposed < 0.0)
                    throw std::runtime_error("Momentum projection would create an invalid PDF");
            }
        }

        for (const auto &recipient : recipients) {
            const auto du = projection.velocityIncrement(recipient.position);
            for (auto d = LatticeModelStencil_T::begin(); d != LatticeModelStencil_T::end(); ++d) {
                const Eigen::Vector3d direction{double(d.cx()), double(d.cy()), double(d.cz())};
                recipient.pdfs->get(recipient.cell.x(), recipient.cell.y(), recipient.cell.z(), d.toIdx()) +=
                    real_c(capillary_coupling::momentumPdfIncrement(double(LatticeModel_T::w[d.toIdx()]),
                                                                    recipient.density, direction, du));
            }
        }
        pdfCommunication();
        return maxCorrection;
    };

    // A numerical conservation projection, not an alternative physical capillary
    // law. Keep waLBerla's Laplace-pressure treatment. Correct only the NET
    // residual left after that treatment, so the full line load is not applied
    // twice.
    timeloop.add() << Sweep(
        [&](IBlock *) {
            stepLoad = ContactLineLoad{};
            stepLine = SurfaceLineResult{};

            // Fluid angular momentum currently uses particle 0's center
            // as its reference point.
            const Eigen::Vector3d angularOrigin = solidCenterEigen;

            for (std::size_t i = 0; i < particles.size(); ++i) {
                auto &particle = particles[i];

                auto &line = particleStepLines[i];

                line = reconstructSurfaceLine(particle.surface, interfacePlanes, particle.center,
                                              double(surfaceTension), 4.5, double(contactAngle));

                // A small/rapidly moving contact patch may be under-resolved
                // on the geometric sphere even with valid fluid data. Refine
                // its actual triangles, then require a closed contour again.
                // Never close gaps with fabricated segments or reuse an old load.
                while ((line.invalidSegments != 0 || line.badEndpoints != 0) && particle.surfaceLevel < 6) {
                    ++particle.surfaceLevel;
                    particle.surface =
                        makeSphereSurface(particle.center, double(solidRadius), particle.surfaceLevel);
                    line = reconstructSurfaceLine(particle.surface, interfacePlanes, particle.center,
                                                  double(surfaceTension), 4.5, double(contactAngle));
                    WALBERLA_LOG_INFO_ON_ROOT("Refined particle surface: step="
                                              << timeloop.getCurrentTimeStep() << " particle=" << i
                                              << " level=" << particle.surfaceLevel << " segments="
                                              << line.segments.size() << " badEndpoints=" << line.badEndpoints
                                              << " invalid=" << line.invalidSegments);
                }
                if (line.invalidSegments != 0 || line.badEndpoints != 0) {
                    std::ofstream dump("contact_line_failure.txt");
                    dump << std::setprecision(17) << particle.center.transpose() << ' ' << solidRadius << ' '
                         << surfaceTension << ' ' << contactAngle << '\n';
                    for (const auto &plane : interfacePlanes)
                        dump << plane.point.transpose() << ' ' << plane.normal.transpose() << '\n';
                    dump.close();
                    std::ostringstream reason;
                    reason << "Invalid contact line: step=" << timeloop.getCurrentTimeStep()
                           << " particle=" << i << " segments=" << line.segments.size()
                           << " badEndpoints=" << line.badEndpoints << " invalid=" << line.invalidSegments
                           << " unsupported=" << line.unsupportedTriangles
                           << " degenerate=" << line.degenerateTriangles;
                    throw std::runtime_error(reason.str());
                }

                auto &load = particleStepLoads[i];

                load.force = particleWallForcesLB[i] + line.load.force;

                // This torque is about this particle's own center.
                load.torque = particleWallTorquesLB[i] + line.load.torque;

                if (!load.force.allFinite() || !load.torque.allFinite()) {
                    throw std::runtime_error("Non-finite load for particle " + std::to_string(i));
                }

                // Aggregate loads about the SAME reference point used
                // by the fluid angular-momentum measurement.
                stepLoad.force += load.force;

                stepLoad.torque += load.torque + (particle.center - angularOrigin).cross(load.force);

                // Combined contact-line geometry for existing output.
                stepLine.segments.insert(stepLine.segments.end(), line.segments.begin(), line.segments.end());

                stepLine.length += line.length;
                stepLine.unsupportedTriangles += line.unsupportedTriangles;
                stepLine.degenerateTriangles += line.degenerateTriangles;
                stepLine.invalidSegments += line.invalidSegments;
                stepLine.badEndpoints += line.badEndpoints;

                stepLine.load.force += line.load.force;

                stepLine.load.torque +=
                    line.load.torque + (particle.center - angularOrigin).cross(line.load.force);
            }

            const auto after = measureLiquidState();
            rawResidual.linear = after.linear - stepBefore.linear + stepLoad.force;
            rawResidual.angular = after.angular - stepBefore.angular + stepLoad.torque;
            maxAppliedCorrectionLB = 0.0;
            const double linearTolerance =
                5e-11 * std::max({1.0, stepBefore.linear.norm(), after.linear.norm(), stepLoad.force.norm()});
            const double angularTolerance = 5e-11 * std::max({1.0, stepBefore.angular.norm(),
                                                              after.angular.norm(), stepLoad.torque.norm()});

            if (rawResidual.linear.norm() > linearTolerance ||
                rawResidual.angular.norm() > angularTolerance) {
                maxAppliedCorrectionLB =
                    applyFluidMomentumCorrection(-rawResidual.linear, -rawResidual.angular);
            }

            const auto corrected = measureLiquidState();
            correctedResidual.linear = corrected.linear - stepBefore.linear + stepLoad.force;
            correctedResidual.angular = corrected.angular - stepBefore.angular + stepLoad.torque;
            if (correctedResidual.linear.norm() > linearTolerance ||
                correctedResidual.angular.norm() > angularTolerance)
                throw std::runtime_error("Fluid-particle momentum projection failed "
                                         "its conservation check");
        },
        "Coupling: conserve resolved linear and angular momentum");

    // add evaluator for total and excessive mass (mass that is currently
    // undistributed)
    const auto totalMass = std::make_shared<real_t>(real_c(0));

    const auto excessMass = std::make_shared<real_t>(real_c(0));

    TotalMassComputer<FreeSurfaceBoundaryHandling_T, PdfField_T, FlagField_T, ScalarField_T>
        totalMassComputer(blockForest, freeSurfaceBoundaryHandling, pdfFieldID, fillFieldID,
                          dynamicsHandler.getConstExcessMassFieldID(),
                          uint_c(1), // Update every LB step.
                          totalMass, excessMass);

    // Measure the initial state before advancing LB.
    totalMassComputer();

    const double initialLiquidMassLB = double(*totalMass);

    if (!std::isfinite(initialLiquidMassLB) || initialLiquidMassLB <= 0.0) {
        throw std::runtime_error("Initial liquid mass must be finite and positive");
    }

    timeloop.addFuncAfterTimeStep(totalMassComputer, "Evaluator: total mass");

    // add VTK output
    addVTKOutput<LatticeModel_T, FreeSurfaceBoundaryHandling_T, PdfField_T, FlagField_T, ScalarField_T,
                 VectorField_T>(blockForest, timeloop, walberlaEnv.config(), flagInfo, pdfFieldID,
                                flagFieldID, fillFieldID, forceDensityFieldID,
                                geometryHandler.getCurvatureFieldID(), geometryHandler.getNormalFieldID(),
                                geometryHandler.getObstNormalFieldID());

    // add triangle mesh output of free surface
    CapillarySurfaceMeshWriter<ScalarField_T, FlagField_T> surfaceMeshWriter(
        blockForest, fillFieldID, flagFieldID, flagIDs::liquidInterfaceGasFlagIDs, real_c(0),
        walberlaEnv.config());
    surfaceMeshWriter(); // write initial mesh
    timeloop.addFuncAfterTimeStep(surfaceMeshWriter, "Writer: surface mesh");

    const uint_t particleWriteFrequency = uint_c(dem.setup->particleFrequency);

    auto clearGeneratedFrames = [](const std::filesystem::path &folder, const std::string &prefix) {
        std::filesystem::create_directories(folder);
        for (const auto &entry : std::filesystem::directory_iterator(folder)) {
            const auto name = entry.path().filename().string();
            if (entry.is_regular_file() && name.starts_with(prefix) && entry.path().extension() == ".vtp")
                std::filesystem::remove(entry.path());
        }
    };
    if (particleWriteFrequency != uint_c(0)) {
        clearGeneratedFrames("particle-out", "particle_");
        clearGeneratedFrames("contact-line-out", "contact_line_");
        // Preserve the same time-series origin as the initial liquid/particle
        // frames. Geometry normals are not available yet, so this frame has no line
        // segments.
        writeContactLineVtp(stepLine,
                            std::filesystem::path("contact-line-out") / "contact_line_000000000000.vtp");
    }

    uint_t particleOutputCounter = 0;

    auto writeParticleMesh = [&]() {
        if (particleWriteFrequency == uint_c(0))
            return;

        const uint_t frame = particleOutputCounter++;

        if (frame % particleWriteFrequency != uint_c(0))
            return;

        const std::filesystem::path folder("particle-out");
        std::filesystem::create_directories(folder);

        std::ostringstream filename;
        filename << "particle_" << std::setfill('0') << std::setw(12) << frame << ".vtp";

        std::ofstream out(folder / filename.str());

        if (!out) {
            throw std::runtime_error("Cannot open particle visualization output");
        }

        out << std::setprecision(17);

        std::size_t triangleCount = 0;

        for (const auto &particle : particles) {
            triangleCount += particle.surface.size();
        }

        const std::size_t pointCount = 3 * triangleCount;

        out << "<?xml version=\"1.0\"?>\n"
            << "<VTKFile type=\"PolyData\" version=\"0.1\" "
               "byte_order=\"LittleEndian\">\n"
            << "<PolyData>\n"
            << "<Piece NumberOfPoints=\"" << pointCount
            << "\" NumberOfVerts=\"0\" NumberOfLines=\"0\" "
               "NumberOfStrips=\"0\" NumberOfPolys=\""
            << triangleCount << "\">\n";

        out << "<CellData Scalars=\"particle_id\">\n"
            << "<DataArray type=\"Int32\" Name=\"particle_id\" "
               "NumberOfComponents=\"1\" format=\"ascii\">\n";

        for (std::size_t i = 0; i < particles.size(); ++i) {
            for (std::size_t triangle = 0; triangle < particles[i].surface.size(); ++triangle) {
                out << i << '\n';
            }
        }

        out << "</DataArray>\n</CellData>\n";

        out << "<Points>\n"
            << "<DataArray type=\"Float64\" "
               "NumberOfComponents=\"3\" format=\"ascii\">\n";

        for (const auto &particle : particles) {
            for (const auto &triangle : particle.surface) {
                for (const auto &point : triangle) {
                    out << point.x() << ' ' << point.y() << ' ' << point.z() << '\n';
                }
            }
        }

        out << "</DataArray>\n</Points>\n"
            << "<Polys>\n"
            << "<DataArray type=\"Int64\" "
               "Name=\"connectivity\" format=\"ascii\">\n";

        for (std::size_t i = 0; i < triangleCount; ++i) {
            out << 3 * i << ' ' << 3 * i + 1 << ' ' << 3 * i + 2 << '\n';
        }

        out << "</DataArray>\n"
            << "<DataArray type=\"Int64\" "
               "Name=\"offsets\" format=\"ascii\">\n";

        for (std::size_t i = 0; i < triangleCount; ++i)
            out << 3 * (i + 1) << '\n';

        out << "</DataArray>\n</Polys>\n"
            << "</Piece>\n</PolyData>\n</VTKFile>\n";

        out.close();

        if (!out) {
            throw std::runtime_error("Failed writing particle visualization output");
        }
    };

    // Match the liquid mesh writer's initial frame and cadence.
    writeParticleMesh();

    timeloop.addFuncAfterTimeStep(writeParticleMesh, "Writer: particle mesh");

    uint_t contactLineOutputCounter = uint_c(0);
    timeloop.addFuncAfterTimeStep(
        [&]() {
            const uint_t frame = ++contactLineOutputCounter;
            if (particleWriteFrequency == uint_c(0) || frame % particleWriteFrequency != uint_c(0))
                return;
            const std::filesystem::path folder("contact-line-out");
            std::filesystem::create_directories(folder);
            std::ostringstream filename;
            filename << "contact_line_" << std::setfill('0') << std::setw(12) << frame << ".vtp";
            writeContactLineVtp(stepLine, folder / filename.str());
        },
        "Writer: applied contact line");

    WcTimingPool timingPool;

    DemLbParameters::Field forceSumSI(particles.size(), Eigen::Vector3d::Zero());

    DemLbParameters::Field torqueSumSI(particles.size(), Eigen::Vector3d::Zero());

    if (!dem.advanceDemStep || !dem.applyRemapImpulseSI || !dem.submitFluidLoadSI) {
        throw std::runtime_error("Missing DEM coupling callbacks");
    }

    // Verify the complete exchange, including the actual DEM integrator and
    // remapper.
    const double linearImpulseScaleSI = units.forceScaleSI * units.dtSI;
    const double angularImpulseScaleSI = units.torqueScaleSI * units.dtSI;
    auto measureSystemMomentum = [&]() -> capillary_coupling::MomentumState {
        const auto fluid = measureLiquidState();

        capillary_coupling::MomentumState systemState;

        systemState.linear = fluid.linear;

        // measureLiquidState() measures angular momentum about
        // solidCenterEigen. Convert it to the fixed LB origin.
        systemState.angular = fluid.angular + solidCenterEigen.cross(fluid.linear);

        for (std::size_t i = 0; i < particles.size(); ++i) {
            const Eigen::Vector3d particleMomentumLB =
                dem.particleMassSI * dem.particleVelocitySI(i) / linearImpulseScaleSI;

            const Eigen::Vector3d particleSpinLB =
                dem.particleInertiaSI * dem.particleAngularVelocitySI(i) / angularImpulseScaleSI;

            systemState.linear += particleMomentumLB;

            const Eigen::Vector3d actualCenterLB = demPositionToLb(dem.particlePositionSI(i));

            systemState.angular += actualCenterLB.cross(particleMomentumLB) + particleSpinLB;
        }

        return systemState;
    };

    const auto initialSystemState = measureSystemMomentum();

    const Eigen::Vector3d initialSystemMomentumLB = initialSystemState.linear;

    const Eigen::Vector3d initialSystemAngularLB = initialSystemState.angular;

    for (uint_t t = uint_c(0); t != timesteps; ++t) {
        stepBefore = measureLiquidState();

        timeloop.singleStep(timingPool, true);

        const Eigen::Vector3d totalFluidForceSI = units.forceScaleSI * stepLoad.force;
        const Eigen::Vector3d totalFluidTorqueSI = units.torqueScaleSI * stepLoad.torque;

        for (std::size_t i = 0; i < particles.size(); ++i) {
            forceSumSI[i] += units.forceScaleSI * particleStepLoads[i].force;

            torqueSumSI[i] += units.torqueScaleSI * particleStepLoads[i].torque;
        }

        if (diagnosticFrequency != uint_c(0) && t % diagnosticFrequency == uint_c(0)) {
            const double relativeMassChange =
                (double(*totalMass) - initialLiquidMassLB) / initialLiquidMassLB;
            WALBERLA_LOG_INFO_ON_ROOT("Coupling: step="
                                      << t << " forceSI_N=" << totalFluidForceSI.transpose()
                                      << " torqueSI_Nm=" << totalFluidTorqueSI.transpose()
                                      << " rawMomentumResidualLB=" << rawResidual.linear.norm()
                                      << " rawAngularResidualLB=" << rawResidual.angular.norm()
                                      << " momentumResidualLB=" << correctedResidual.linear.norm()
                                      << " angularResidualLB=" << correctedResidual.angular.norm()
                                      << " maxVelocityCorrectionLB=" << maxAppliedCorrectionLB
                                      << " relativeMassChange=" << relativeMassChange << " minFill="
                                      << minimumResolvedFillLB << " maxFill=" << maximumResolvedFillLB
                                      << " lineSegments=" << stepLine.segments.size()
                                      << " wallLinks=" << sphereLinkCount
                                      << " wallMassCorrectionLB=" << movingWallMassCorrectionLB);
        }

        if ((t + uint_c(1)) % units.lbStepsPerDemStep == uint_c(0)) {
            const double count = double(units.lbStepsPerDemStep);

            DemLbParameters::Field meanForcesSI(particles.size());
            DemLbParameters::Field meanTorquesSI(particles.size());

            for (std::size_t i = 0; i < particles.size(); ++i) {
                meanForcesSI[i] = forceSumSI[i] / count;
                meanTorquesSI[i] = torqueSumSI[i] / count;
            }

            dem.submitFluidLoadSI(meanForcesSI, meanTorquesSI);

            for (std::size_t i = 0; i < particles.size(); ++i) {
                forceSumSI[i].setZero();
                torqueSumSI[i].setZero();
            }

            if (advanceDEM)
                dem.advanceDemStep();

            if (advanceDEM) {
                // The force is held constant through one Verlet step. Its
                // orbital angular impulse is r_mid x J, while the LB sweep
                // used the fixed old center r_old. Transfer their difference
                // to the FLUID; do not invent a particle spin correction.
                Eigen::Vector3d leverTransportLB = Eigen::Vector3d::Zero();
                for (std::size_t i = 0; i < particles.size(); ++i) {
                    const Eigen::Vector3d displacementLB =
                        demPositionToLb(dem.particlePositionSI(i)) - particles[i].center;
                    const Eigen::Vector3d impulseLB =
                        meanForcesSI[i] * dem.demTimeStepSI / linearImpulseScaleSI;
                    leverTransportLB += 0.5 * displacementLB.cross(impulseLB);
                }
                if (!leverTransportLB.allFinite())
                    throw std::runtime_error("Non-finite moving force-lever impulse");
                if (leverTransportLB.squaredNorm() > 0.0) {
                    const auto before = measureLiquidState();
                    applyFluidMomentumCorrection(Eigen::Vector3d::Zero(), -leverTransportLB);
                    const auto after = measureLiquidState();
                    if ((after.linear - before.linear).norm() > 5e-11 ||
                        (after.angular - before.angular + leverTransportLB).norm() > 5e-10)
                        throw std::runtime_error("Moving force-lever transfer failed");

                }
            }

            // Validate every proposed position before modifying the fluid.
            std::vector<Eigen::Vector3d> newCentersLB(particles.size());

            for (std::size_t i = 0; i < particles.size(); ++i) {
                newCentersLB[i] = demPositionToLb(dem.particlePositionSI(i));

                const Eigen::Vector3d displacementLB = newCentersLB[i] - particles[i].center;

                if (!newCentersLB[i].allFinite() || displacementLB.norm() > 0.1) {
                    throw std::runtime_error("Particle " + std::to_string(i) +
                                             " moved more than 0.1 lattice cells; "
                                             "reduce the DEM timestep");
                }

                for (int axis = 0; axis < 3; ++axis) {
                    const double margin = double(solidRadius) + 2.0;

                    if (newCentersLB[i][axis] <= margin ||
                        newCentersLB[i][axis] >= double(domainSize[uint_c(axis)]) - margin) {
                        throw std::runtime_error("Particle " + std::to_string(i) +
                                                 " is too close to a domain boundary");
                    }
                }
            }

            // Remap sequentially. Each move preserves the other sphere
            // at its currently committed position.
            for (std::size_t i = 0; i < particles.size(); ++i) {
                const Eigen::Vector3d newCenterLB = newCentersLB[i];

                RemapResult result{Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero(), particles[i].center};

                const auto changes = detectParticleCellChanges(i, newCenterLB);

                if (changes.covered != 0 || changes.uncovered != 0) {
                    result = remapParticle(i, newCenterLB);

                    communication();
                    pdfCommunication();
                } else {
                    // Subcell motion without a change in solid occupancy.
                    particles[i].center = newCenterLB;

                    particles[i].surface =
                        makeSphereSurface(newCenterLB, double(solidRadius), particles[i].surfaceLevel);
                }

                // Convert angular impulse from the old center
                // to the particle's new center.
                const Eigen::Vector3d angularImpulseAtNewCenterLB =
                    result.angularImpulseAboutOldCenterLB -
                    (newCenterLB - result.oldCenterLB).cross(result.linearImpulseLB);

                const Eigen::Vector3d linearImpulseSI = linearImpulseScaleSI * result.linearImpulseLB;

                const Eigen::Vector3d angularImpulseSI = angularImpulseScaleSI * angularImpulseAtNewCenterLB;

                if (advanceDEM) {
                    dem.applyRemapImpulseSI(i, linearImpulseSI, angularImpulseSI);
                }
            }

            // Synchronize the remaining particle-0 aliases.
            solidCenterEigen = particles[0].center;

            solidCenter = Vector3<real_t>(real_c(solidCenterEigen.x()), real_c(solidCenterEigen.y()),
                                          real_c(solidCenterEigen.z()));

            solidSurface = particles[0].surface;

            // Ownership can change even when the union of occupied cells
            // stays unchanged, so refresh every solid cell's velocity.
            updateSphereWallVelocity();

            // Include any remapping in the current mass measurement.
            totalMassComputer();

            if (advanceDEM) {
                const auto systemState = measureSystemMomentum();

                const Eigen::Vector3d systemMomentumLB = systemState.linear;

                const Eigen::Vector3d systemAngularLB = systemState.angular;

                const double linearErrorLB = (systemMomentumLB - initialSystemMomentumLB).norm();
                const double angularErrorLB = (systemAngularLB - initialSystemAngularLB).norm();
                if (diagnosticFrequency != uint_c(0) &&
                    (t % diagnosticFrequency == uint_c(0) || t + uint_c(1) == timesteps)) {
                    WALBERLA_LOG_INFO_ON_ROOT("System momentum: step="
                                              << t << " linearErrorLB=" << linearErrorLB
                                              << " angularErrorLB=" << angularErrorLB << " combinedP_SI="
                                              << (linearImpulseScaleSI * systemMomentumLB).transpose());
                }
                if (!std::isfinite(linearErrorLB) || !std::isfinite(angularErrorLB) ||
                    linearErrorLB > 1e-9 * std::max(1.0, initialSystemMomentumLB.norm()) ||
                    angularErrorLB > 1e-8 * std::max(1.0, initialSystemAngularLB.norm()))
                    throw std::runtime_error("Full fluid-DEM momentum conservation failed after coupling");
            }
        }
        const double relativeMassChange = (double(*totalMass) - initialLiquidMassLB) / initialLiquidMassLB;
        if (!std::isfinite(relativeMassChange) || std::abs(relativeMassChange) > 1e-8)
            throw std::runtime_error("Liquid mass conservation failed during coupling");
    }
    WALBERLA_LOG_INFO_ON_ROOT("Completed " << timesteps << " LB steps; advanceDEM=" << advanceDEM);

    return EXIT_SUCCESS;
}
} // namespace DropWetting
} // namespace free_surface
} // namespace walberla

int runCapillaryLb(int argc, char **argv, const DemLbParameters &dem) {
    return walberla::free_surface::DropWetting::runWithDem(argc, argv, dem);
}
