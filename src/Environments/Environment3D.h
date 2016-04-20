//
// Created by Max Horn on 19/04/16.
//

#ifndef CHEMOHYBRID_GPU_ENVIRONMENT3D_H
#define CHEMOHYBRID_GPU_ENVIRONMENT3D_H

#include "Environment.h"

class Environment3D : Environment {
    // Boundary condition functions;
    static void applyNeumannBC(array input, double resolution, BoundaryCondition bc);
    static array getLaplacian();
public:
    Environment3D(EnvironmentSettings settings);
};


#endif //CHEMOHYBRID_GPU_ENVIRONMENT3D_H
