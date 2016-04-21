//
// Created by Max Horn on 19/04/16.
//

#ifndef CHEMOHYBRID_GPU_ENVIRONMENT2D_H
#define CHEMOHYBRID_GPU_ENVIRONMENT2D_H

#include "Environment.h"

class Environment2D : public Environment {
    // Boundary condition functions;
    static void applyNeumannBC(array *input, double resolution, BoundaryCondition *bc);
    static void serialCalculateTimeStep(array *densities, array *diffusionFilters, double dt);
    static void batchCalculateTimeStep(array *densities, array *diffusionFilters, double dt);
    static array getLaplacian();
public:
    Environment2D(EnvironmentSettings settings);
    array getAllDensities() override;
    array getDensity(unsigned int ligand) override;
    dim4 getSize() override;
    void test();
};

#endif //CHEMOHYBRID_GPU_ENVIRONMENT2D_H
