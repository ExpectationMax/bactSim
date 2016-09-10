//
// Created by Max Horn on 19/04/16.
//

#ifndef CHEMOHYBRID_GPU_ENVIRONMENT2D_H
#define CHEMOHYBRID_GPU_ENVIRONMENT2D_H

#include "Environment.h"
#include <vector>
#include "Solvers/Solver.h"
#include <H5Cpp.h>

#define I_TOPLEFT 0
#define I_TOPRIGHT 1
#define I_BOTTOMLEFT 2
#define I_BOTTOMRIGHT 3
#define W_TOPLEFT 0
#define W_TOPRIGHT 1
#define W_BOTTOMLEFT 2
#define W_BOTTOMRIGHT 3


class Environment2D : public Environment {
    void init();
    // Boundary condition functions;
    static void applyNeumannBC(array *input, double resolution, BoundaryCondition *bc);
    static void applyDericheletBC(array *input,  BoundaryCondition *bc);
    static void applyPeriodicBC(array *input);
    static array getLaplacian();

    class Diffusion2D : public DifferentialEquation {
        Environment2D *parent;
    public:
        Diffusion2D(Environment2D *par): parent(par) {}
        virtual array rateofchange(array &input) override;
    };

    std::map<unsigned int, unique_ptr<H5::DataSet>> ligands_storage;

    array get_concentrations(array &indexes, array &ligands);

public:
    Environment2D(EnvironmentSettings settings, shared_ptr<Solver> solver);
    array getAllDensities() override;
    array getDensity(int) override;
    std::vector<double> getSize() override;

    void setInterpolatedPositions(array &xpos, array &ypos, array &pos, array &weights);

    void changeLigandConcentrationBy(array concDifferences, array positions, array weights, array ligands);

    array getLigandConcentrations(array positions, array weights, array ligands);

    void evalDensities();

    void closeStorage();

    void setupStorage(unique_ptr<H5::Group> storage);

    void save();

    Environment2D(H5::Group group);
};

#endif //CHEMOHYBRID_GPU_ENVIRONMENT2D_H
