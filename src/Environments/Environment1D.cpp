//
// Created by Max Horn on 19/04/16.
//

#include "Environment1D.h"

array Environment1D::getLaplacian() {
    GPU_REALTYPE data1[] = {1.0, -2.0, 1.0};
    return array(3, data1);
}

Environment1D::Environment1D(EnvironmentSettings settings) : Environment(settings) {
    densities = array(internal_dimensions[0], internal_dimensions[1], settings.dataType);
    density_changes = constant(0.0, internal_dimensions[0], internal_dimensions[1], settings.dataType);
    diffusion_filters = constant(0.0, LAPLACIAN_SIZE, (dim_t)settings.ligands.size());

    // Initialize ligand concentrations and diffusion filters
    for(size_t i = 0; i < settings.ligands.size(); i++) {
        densities(span, i) = settings.ligands[i].initialConcentration;
        diffusion_filters(span, i) = Environment1D::getLaplacian();
        diffusion_filters(span, i) *= settings.ligands[i].diffusionCoefficient;
    }

    switch(this->boundaryCondition.type) {
        case BC_NEUMANN:
        default:
            this->applyBoundaryCondition = std::bind(Environment1D::applyNeumannBC, this->densities, this->resolution, this->boundaryCondition);
    }
}
