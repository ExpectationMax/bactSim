//
// Created by Max Horn on 18/04/16.
//

#ifndef PROJECT_NAME_ENVIRONMENT_H
#define PROJECT_NAME_ENVIRONMENT_H

#define GPU_REALTYPE float
#define BORDER_SIZE 1
#define LAPLACIAN_SIZE 1 + 2 * BORDER_SIZE

#define LIGANDID 0
#define LIGANDINTERNAL 1

#include <functional>
#include <array>
#include <arrayfire.h>
#include <type_traits>
#include <map>
#include <General/Ligand.hpp>
#include "Solvers/Solver.h"
#include <memory>

using namespace af;
using std::unique_ptr;

enum BoundaryConditionType {
    BC_PERIODIC,
    BC_DIRICHELET,
    BC_NEUMANN
};

class BoundaryCondition {
public:
    // Default allways use Neumann boundaries
    BoundaryCondition() {
        this->type = BC_NEUMANN;
    };
    BoundaryCondition(BoundaryConditionType boundaryConditionType) {
        this->type = boundaryConditionType;
    };
    BoundaryCondition(BoundaryConditionType boundaryConditionType, double xpos, double xneg) :
            BoundaryCondition(boundaryConditionType) {
        this->xpos = xpos;
        this->xneg = xneg;
    };
    BoundaryCondition(BoundaryConditionType boundaryConditionType, double xpos, double xneg, double ypos, double yneg) :
            BoundaryCondition(boundaryConditionType, xpos, xneg) {
        this->ypos = ypos;
        this->yneg = yneg;
    };
    BoundaryCondition(BoundaryConditionType boundaryConditionType, double xpos, double xneg, double ypos, double yneg,
                      double zpos, double zneg) : BoundaryCondition(boundaryConditionType, xpos, xneg, ypos, yneg) {
        this->zpos = zpos;
        this->zneg = zneg;
    };
    BoundaryConditionType type;
    double xneg, xpos, yneg, ypos, zneg, zpos = 0;
};

enum ConvolutionType {
    CT_AFBATCH,
    CT_SERIAL
};


struct EnvironmentSettings {
    // Definition of Size
    double resolution;
    std::vector<double> dimensions;

    // Definition of ligands
    std::vector<Ligand> ligands;

    //Simulation parameters
    GPU_REALTYPE dt;
    af_dtype dataType;
    BoundaryCondition boundaryCondition;
    ConvolutionType convolutionType = CT_SERIAL;

    // Visualization parameters
};


class Environment {
protected:
    // Internal arrays
    af::array densities;
    af::array diffusion_filters;

    std::vector<Ligand> ligands;
    array ligandMapping;

    // Simultation Parameters
    GPU_REALTYPE dt;
    BoundaryCondition boundaryCondition;
    double resolution;
    std::vector<dim_t> internal_dimensions;
#ifndef NO_GRAPHICS
    Window *visualizationWin;
    unsigned int numLigands;
    unsigned int rows;
#endif
    std::function<void(void)> applyBoundaryCondition;
    //std::function<void(void)> calculateTimeStep;
    void calculateTimeStep() {
        densities = odesolver->solveStep(*diffusionEquation, densities, dt);
    }

    static void batchCalculateTimeStep(array densities, array densityChange, array diffusionFilters, double dt);

    Solver *odesolver;
    unique_ptr<DifferentialEquation> diffusionEquation = NULL;

    Environment(EnvironmentSettings settings, Solver &odesolver);
public:
    std::vector<Ligand> getLigands();
    void printInternals();
    void simulate(double advanceTime);
    void simulateTimeStep(void);
#ifndef NO_GRAPHICS
    void visualize(double normalizer);
#endif
    virtual std::vector<double> getSize() = 0;
    //virtual void simulateTimeStep() = 0;
    virtual array getDensity(int) = 0;
    virtual array getAllDensities() = 0;
    //void test();
    array getLigandMapping(std::vector<int> ligands);

    void setupVisualizationWindow(Window &win);

    BoundaryConditionType getBoundaryConditionType() { return boundaryCondition.type; }
};

#endif //PROJECT_NAME_ENVIRONMENT_H
