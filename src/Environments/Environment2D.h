//
// Created by Max Horn on 19/04/16.
//

#ifndef CHEMOHYBRID_GPU_ENVIRONMENT2D_H
#define CHEMOHYBRID_GPU_ENVIRONMENT2D_H

#include "Environment.h"

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
public:
    Environment2D(EnvironmentSettings settings);
    array getAllDensities() override;
    array getDensity(int) override;
    dim4 getSize() override;
    void test();

    array getInterpolatedPositions(array &xpos, array &ypos);

    void changeLigandConcentrationBy(array concDifferences, array posAndWeights, array ligands);

    array getLigandConcentrations(array posAndWeights, array ligands);
};

#endif //CHEMOHYBRID_GPU_ENVIRONMENT2D_H
