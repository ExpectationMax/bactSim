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
    RealType uptakeRate;
    RealType productionRate;
    RealType diffusionCoefficient;
};

enum BorderCondition {
    BC_PERIODIC,
    BC_DIRICHELET,
    BC_NEUMANN
};

class BorderSetting {
    BorderCondition borderCondition;
    void applyBorderCondition(af::array);
};

template <typename RealType>
struct EnvironmentSettings {
    // Definition of Size
    double resolution;
    std::vector<double> dimensions;

    // Definition of ligands
    std::vector<Ligand<RealType>> ligands;

    //Simulation parameters
    double dt;
    af_dtype dataType;
};

template <typename RealType>
class Environment {
    // Only allow floating point types
    static_assert(std::is_floating_point<RealType>::value,
                  "class Environment can only be instantiated with floating point types");
    af::array densities;
    af::array density_changes;
    af::array diffusion_filters;
    af::array production_rates;
    af::array uptake_rates;
    std::vector<Ligand<RealType>> ligands;

public:
    Environment(EnvironmentSettings<RealType> settings);
    array getLaplacian(unsigned int dims);
    std::vector<Ligand<RealType>> getLigands();
    std::vector<dim_t> getSize();
    void printInternals();
    void simulate(double advanceTime);
    void simulateTimeStep();
    void test();
    void load_densitydistribution(unsigned int ligand, RealType *);
    void applyBoundaryCondition();
};

#endif //PROJECT_NAME_ENVIRONMENT_H
