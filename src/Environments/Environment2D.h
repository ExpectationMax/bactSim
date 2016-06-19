//
// Created by Max Horn on 19/04/16.
//

#ifndef CHEMOHYBRID_GPU_ENVIRONMENT2D_H
#define CHEMOHYBRID_GPU_ENVIRONMENT2D_H

#include "Environment.h"
#include "Solvers/Solver.h"

#define POS_LEFT 0
#define POS_RIGHT 1
#define POS_TOP 2
#define POS_BOTTOM 3
#define W_TOPLEFT 4
#define W_TOPRIGHT 5
#define W_BOTTOMLEFT 6
#define W_BOTTOMRIGHT 7


class Environment2D : public Environment {
    // Boundary condition functions;
    static void applyNeumannBC(array *input, double resolution, BoundaryCondition *bc);
    static void applyDericheletBC(array *input,  BoundaryCondition *bc);
    static void applyPeriodicBC(array *input);
    static void serialCalculateTimeStep(array *densities, array *diffusionFilters, double dt, std::vector<Ligand> *ligands);
    static void batchCalculateTimeStep(array *densities, array *diffusionFilters, double dt, std::vector<Ligand> *ligands);
    static array getLaplacian();

    class Diffusion2D : DifferentialEquation {
        Environment2D *parent;
    public:
        Diffusion2D(Environment2D *par): parent(par){}
        array rateofchange(array input) {
            array changes = array(parent->densities.dims());
            for (size_t i = 0; i < parent->densities.dims(2); i++) {
                changes(seq(1, end-1), seq(1, end-1), i) = convolve(input(span, span, i),
                                               parent->diffusion_filters(span, span, i))(seq(1, end-1), seq(1, end-1));

                changes(seq(1, end-1), seq(1, end-1), i) += parent->ligands[i].globalProductionRate
                                                            - parent->ligands[i].globalDegradationRate*input(seq(1, end-1), seq(1, end-1), i);
            }
            changes.eval();
            return changes;
        };
    };

    Diffusion2D diffequation;
public:
    Environment2D(EnvironmentSettings settings);
    array getAllDensities() override;
    array getDensity(int) override;
    std::vector<double> getSize() override;
    void test();

    array getInterpolatedPositions(array &xpos, array &ypos);

    void changeLigandConcentrationBy(array concDifferences, array posAndWeights, array ligands);

    array getLigandConcentrations(array posAndWeights, array ligands);

    void evalDensities();
};

#endif //CHEMOHYBRID_GPU_ENVIRONMENT2D_H
