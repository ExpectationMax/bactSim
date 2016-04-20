//
// Created by Max Horn on 19/04/16.
//

#ifndef CHEMOHYBRID_GPU_ENVIRONMENT1D_H
#define CHEMOHYBRID_GPU_ENVIRONMENT1D_H

#include "Environment.h"

class Environment1D : Environment {
    // Boundary condition functions;
    static void applyNeumannBC(array input, double resolution, BoundaryCondition bc);
    static array getLaplacian();
public:
    Environment1D(EnvironmentSettings settings);
};


#endif //CHEMOHYBRID_GPU_ENVIRONMENT1D_H
