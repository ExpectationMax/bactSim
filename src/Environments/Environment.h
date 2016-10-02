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
    // Definition of Size and Boundary
    double resolution;
    std::vector<double> dimensions;
    BoundaryCondition boundaryCondition;

    // Definition of ligands
    std::vector<Ligand> ligands;
};


class Environment {
private:
    void init(EnvironmentSettings settings);

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
    EnvironmentSettings settings;

    Environment(EnvironmentSettings settings);
    Environment(H5::Group);

    unique_ptr<H5::Group> storage;
public:
    BoundaryCondition boundaryCondition;
    double resolution;
    virtual double getStabledt() = 0;
    virtual void simulateTimestep(double dt) = 0;
#ifndef NO_GRAPHICS
    void visualize(double normalizer);
#endif
    virtual std::vector<double> getSize() = 0;
    virtual array getDensity(int) = 0;
    virtual array getAllDensities() = 0;
    array getLigandMapping(std::vector<int> ligands);
    void setupVisualizationWindow(Window &win);
    virtual BoundaryConditionType getBoundaryConditionType() { return boundaryCondition.type; }
    virtual void save() = 0;
    virtual void setupStorage(unique_ptr<H5::Group> unique_ptr);
    virtual void closeStorage();
};

#endif //PROJECT_NAME_ENVIRONMENT_H
