//
// Created by Max Horn on 19/04/16.
//

#include "Environment3D.h"

array Environment3D::getLaplacian() {
    GPU_REALTYPE data3[] =
            {0.0, 0.0, 0.0,
             0.0, 1.0, 0.0,
             0.0, 0.0, 0.0,

             0.0, 1.0, 0.0,
             1.0, -6.0, 1.0,
             0.0, 1.0, 0.0,

             0.0, 0.0, 0.0,
             0.0, 1.0, 0.0,
             0.0, 0.0, 0.0};
    return array(3,3,3, data3);
}

Environment3D::Environment3D(EnvironmentSettings settings) : Environment(settings) {
    densities = array(internal_dimensions[0], internal_dimensions[1],
                      internal_dimensions[2], internal_dimensions[3]);
    density_changes = constant(0.0, internal_dimensions[0], internal_dimensions[1],
                               internal_dimensions[2], internal_dimensions[3]);
    this->diffusion_filters = constant(0.0, LAPLACIAN_SIZE, LAPLACIAN_SIZE, LAPLACIAN_SIZE, (dim_t)settings.ligands.size());
    for(size_t i = 0; i < settings.ligands.size(); i++) {
        densities(span, span, span, i) = settings.ligands[i].initialConcentration;
        diffusion_filters(span, span, span, i) = Environment3D::getLaplacian();
        diffusion_filters(span, span, span, i) *= settings.ligands[i].diffusionCoefficient;
    }

    switch(this->boundaryCondition.type) {
        case BC_NEUMANN:
        default:
            this->applyBoundaryCondition = std::bind(Environment3D::applyNeumannBC, this->densities, this->resolution, this->boundaryCondition);
    }
}