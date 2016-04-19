//
// Created by Max Horn on 18/04/16.
//

#ifndef PROJECT_NAME_ENVIRONMENT_H
#define PROJECT_NAME_ENVIRONMENT_H

#define BORDER_SIZE 1
#define LAPLACIAN_SIZE 1 + 2 * BORDER_SIZE

#include <array>
#include <arrayfire.h>
#include <type_traits>

using namespace af;

template <typename RealType>
struct Ligand {
    std::string name;
    RealType initialConcentration;
    RealType globalProductionRate;
    RealType globalDegradationRate;
    RealType coliUptakeRate;
    RealType coliProductionRate;
    RealType diffusionCoefficient;
};

enum BoundaryCondition {
    BC_PERIODIC,
    BC_DIRICHELET,
    BC_NEUMANN
};

template <typename RealType>
struct EnvironmentSettings {
    // Definition of Size
    double resolution;
    std::vector<double> dimensions;

    // Definition of ligands
    std::vector<Ligand<RealType>> ligands;

    //Simulation parameters
    RealType dt;
    af_dtype dataType;
    BoundaryCondition boundaryCondition;
    // Visualization parameters
    Window *win = NULL;
};

template <typename RealType>
class Environment {
    // Only allow floating point types
    static_assert(std::is_floating_point<RealType>::value,
                  "class Environment can only be instantiated with floating point types");

    // Internal arrays
    af::array densities;
    af::array density_changes;
    af::array diffusion_filters;
    af::array globalProductionRates;

    //
    std::vector<Ligand<RealType>> ligands;

    // Simultation Parameters
    RealType dt;
    BoundaryCondition boundaryCondition;
    double resolution;

    Window *visualizationWin;

    // Boundary condition functions;
    static void _apply_neumann_2d(array input);

public:
    Environment(EnvironmentSettings<RealType> settings);
    static array getLaplacian(unsigned int dims);
    std::vector<Ligand<RealType>> getLigands();
    std::vector<dim_t> getSize();
    void printInternals();
    void simulate(double advanceTime);
    void simulateTimeStep();
    void test();
    void load_densitydistribution(unsigned int ligand, RealType *);
    std::function<void(void)> applyBoundaryCondition;
};

#endif //PROJECT_NAME_ENVIRONMENT_H
