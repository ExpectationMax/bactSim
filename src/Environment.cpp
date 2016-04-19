//
// Created by Max Horn on 18/04/16.
//

#include "Environment.h"
#include <exception>
#include "easylogging++.h"

template <typename RT>
array Environment<RT>::getLaplacian(unsigned dims) {
    array laplacian;
    switch(dims) {
        case 1: {
            RT data1[] = {1.0, -2.0, 1.0};
            laplacian = array(3, data1);
            return laplacian;
        }
        case 2: {
            RT data2 [] =
                    {0.0, 1.0, 0.0,
                     1.0, -4.0, 1.0,
                     0.0, 1.0, 0.0};
            laplacian = array(3, 3, data2);
            return laplacian;
        }
        case 3: {
            RT data3[] =
                                 {0.0, 0.0, 0.0,
                                  0.0, 1.0, 0.0,
                                  0.0, 0.0, 0.0,

                                  0.0, 1.0, 0.0,
                                  1.0, -6.0, 1.0,
                                  0.0, 1.0, 0.0,

                                  0.0, 0.0, 0.0,
                                  0.0, 1.0, 0.0,
                                  0.0, 0.0, 0.0};
            laplacian = array(3,3,3, data3);
            return laplacian;
        }
        default:
            throw exception("Unsupported number of dimensions in laplacian generation");
    }
}

template <typename RT>
Environment<RT>::Environment(EnvironmentSettings<RT> settings) {
    std::vector<dim_t> internal_dimensions(settings.dimensions.size() + 1);

    // Calculate dimensions of internal representations
    //
    for (auto i = 0; i < internal_dimensions.size(); i++) {
        internal_dimensions[i] = (dim_t) 2 * BORDER_SIZE + floor(settings.dimensions[i]/settings.resolution);
    }

    // Last dimension is the number of ligands
    internal_dimensions[internal_dimensions.size() - 1] = settings.ligands.size();

    // Initialize arrays and fill with initial concentrations
    switch(internal_dimensions.size()) {
        case 2:
            this->densities = array(internal_dimensions[0], internal_dimensions[1], settings.dataType);
            this->density_changes = constant(0.0, internal_dimensions[0], internal_dimensions[1], settings.dataType);
            this->diffusion_filters = constant(0.0, LAPLACIAN_SIZE, (dim_t)settings.ligands.size());
            // Initialize ligand concentrations and diffusion filters
            for(size_t i = 0; i < settings.ligands.size(); i++) {
                this->densities(span, i) = settings.ligands[i].initialConcentration;
                this->diffusion_filters(span, i) = this->getLaplacian(1);
                this->diffusion_filters(span, i) *= -settings.ligands[i].diffusionCoefficient;
            }
            break;
        case 3:
            this->densities = array(internal_dimensions[0], internal_dimensions[1], internal_dimensions[2], settings.dataType);
            this->density_changes =constant(0.0, internal_dimensions[0], internal_dimensions[1], internal_dimensions[2], settings.dataType);
            this->diffusion_filters = constant(0.0, LAPLACIAN_SIZE, LAPLACIAN_SIZE, (dim_t)settings.ligands.size());
            for(size_t i = 0; i < settings.ligands.size(); i++) {
                this->densities(span, span, i) = settings.ligands[i].initialConcentration;
                this->diffusion_filters(span, span, i) = this->getLaplacian(2);
                this->diffusion_filters(span, span, i) *= -settings.ligands[i].diffusionCoefficient;
            }
            break;
        case 4:
            this->densities = array(internal_dimensions[0], internal_dimensions[1],
                                    internal_dimensions[2], internal_dimensions[3]);
            this->density_changes = constant(0.0, internal_dimensions[0], internal_dimensions[1],
                                    internal_dimensions[2], internal_dimensions[3]);
            this->diffusion_filters = constant(0.0, LAPLACIAN_SIZE, LAPLACIAN_SIZE, LAPLACIAN_SIZE, (dim_t)settings.ligands.size());
            for(size_t i = 0; i < settings.ligands.size(); i++) {
                this->densities(span, span, span, i) = settings.ligands[i].initialConcentration;
                this->diffusion_filters(span, span, span, i) = this->getLaplacian(3);
                this->diffusion_filters(span, span, span, i) *= -settings.ligands[i].diffusionCoefficient;
            }
            break;
        default:
            throw new exception("Unsupported number of dimensions for grid creation.");
    }
}
template <typename RT>
void Environment<RT>::printInternals() {
    std::cout << "densities:" << std::endl;
    af_print(this->densities);
    std::cout << std::endl;

    std::cout << "diffusion filters:" << std::endl;
    af_print(this->diffusion_filters);
    std::cout << std::endl;
}

template <typename RT>
std::vector<Ligand<RT>> Environment<RT>::getLigands() {
    return this->ligands;
}

template <typename RT>
void Environment<RT>::load_densitydistribution(unsigned int ligand, RT *) {

}

template <typename RT>
void Environment<RT>::test() {
    this->densities(4,4,span) = 255;
}

template <typename RT>
void Environment<RT>::simulate(double advanceTime) {
    while()
}

template <typename RT>
void Environment::simulateTimeStep() {
    this->applyBoundaryCondition();
    this->density_changes = convolve(this->densities, this->diffusion_filters);
    this->densities += this->density_changes;
}

template <typename RT>
void Environment<RT>::applyBoundaryCondition() {
    //stub
}


template class Environment<float>;
template class Environment<double>;