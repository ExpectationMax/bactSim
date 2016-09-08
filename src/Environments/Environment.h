//
// Created by Max Horn on 18/04/16.
//

#ifndef PROJECT_NAME_ENVIRONMENT_H
#define PROJECT_NAME_ENVIRONMENT_H

#include "General/Types.h"
#define BORDER_SIZE 1
#define LAPLACIAN_SIZE 1 + 2 * BORDER_SIZE

#define LIGANDID 0
#define LIGANDINTERNAL 1

#include <functional>
#include <array>
#include <arrayfire.h>
#include <type_traits>
#include <map>
#include <memory>
#include <General/Ligand.h>
#include "Solvers/Solver.h"
#include "BoundaryCondition.h"

using namespace af;
using std::unique_ptr;
using std::shared_ptr;


struct EnvironmentSettings {
    // Definition of Size
    double resolution;
    std::vector<double> dimensions;

    // Definition of ligands
    std::vector<Ligand> ligands;

    //Simulation parameters
    GPU_REALTYPE dt;
    BoundaryCondition boundaryCondition;
};


class Environment {
private:
    void init(EnvironmentSettings settings, shared_ptr<Solver> odesolver);

protected:
    // Internal arrays
    af::array densities;
    af::array diffusion_filters;

    std::vector<Ligand> ligands;
    array ligandMapping;
    std::map<unsigned int, unsigned int> hostLigandMapping;

    // Simultation Parameters
    std::vector<dim_t> internal_dimensions;
#ifndef NO_GRAPHICS
    Window *visualizationWin;
    unsigned int numLigands;
    unsigned int rows;
#endif
    std::function<void(void)> applyBoundaryCondition;
    //std::function<void(void)> calculateTimeStep;
    void calculateTimeStep() {
        odesolver->solveStep(*diffusionEquation, densities, dt);
    }

    shared_ptr<Solver> odesolver;
    unique_ptr<DifferentialEquation> diffusionEquation;

    EnvironmentSettings settings;

    Environment(EnvironmentSettings settings, shared_ptr<Solver> odesolver);
    Environment(H5::Group);

    unique_ptr<H5::Group> storage;
public:
    GPU_REALTYPE dt;
    BoundaryCondition boundaryCondition;
    double resolution;

    void simulateTimeStep(void);
#ifndef NO_GRAPHICS
    void visualize(double normalizer);
#endif
    virtual std::vector<double> getSize() = 0;
    virtual array getDensity(int) = 0;
    virtual array getAllDensities() = 0;
    array getLigandMapping(std::vector<int> ligands);
    void setupVisualizationWindow(Window &win);
    BoundaryConditionType getBoundaryConditionType() { return boundaryCondition.type; }

    void setupStorage(unique_ptr<H5::Group> unique_ptr);
};

#endif //PROJECT_NAME_ENVIRONMENT_H
