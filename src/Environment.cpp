//
// Created by Max Horn on 18/04/16.
//

#include "Environment.h"
#include <exception>
#include "easylogging++.h"
#include <time.h>
#include "progress.h"

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
            densities = array(internal_dimensions[0], internal_dimensions[1], settings.dataType);
            density_changes = constant(0.0, internal_dimensions[0], internal_dimensions[1], settings.dataType);
            diffusion_filters = constant(0.0, LAPLACIAN_SIZE, (dim_t)settings.ligands.size());
            // Initialize ligand concentrations and diffusion filters
            for(size_t i = 0; i < settings.ligands.size(); i++) {
                densities(span, i) = settings.ligands[i].initialConcentration;
                diffusion_filters(span, i) = Environment<RT>::getLaplacian(1);
                diffusion_filters(span, i) *= settings.ligands[i].diffusionCoefficient;
            }
            break;
        case 3:
            densities = array(internal_dimensions[0], internal_dimensions[1], internal_dimensions[2], settings.dataType);
            density_changes =constant(0.0, internal_dimensions[0], internal_dimensions[1], internal_dimensions[2], settings.dataType);
            diffusion_filters = constant(0.0, LAPLACIAN_SIZE, LAPLACIAN_SIZE, (dim_t)settings.ligands.size());
            for(size_t i = 0; i < settings.ligands.size(); i++) {
                densities(span, span, i) = settings.ligands[i].initialConcentration;
                diffusion_filters(span, span, i) = Environment<RT>::getLaplacian(2);
                diffusion_filters(span, span, i) *= settings.ligands[i].diffusionCoefficient;
            }
            break;
        case 4:
            densities = array(internal_dimensions[0], internal_dimensions[1],
                                    internal_dimensions[2], internal_dimensions[3]);
            density_changes = constant(0.0, internal_dimensions[0], internal_dimensions[1],
                                    internal_dimensions[2], internal_dimensions[3]);
            this->diffusion_filters = constant(0.0, LAPLACIAN_SIZE, LAPLACIAN_SIZE, LAPLACIAN_SIZE, (dim_t)settings.ligands.size());
            for(size_t i = 0; i < settings.ligands.size(); i++) {
                densities(span, span, span, i) = settings.ligands[i].initialConcentration;
                diffusion_filters(span, span, span, i) = Environment<RT>::getLaplacian(3);
                diffusion_filters(span, span, span, i) *= settings.ligands[i].diffusionCoefficient;
            }
            break;
        default:
            throw new exception("Unsupported number of dimensions for grid creation.");
    }
    // Setup other properties
    this->dt = settings.dt;
    this->resolution = settings.resolution;
    this->visualizationWin = settings.win;
    this->boundaryCondition = settings.boundaryCondition;
    this->applyBoundaryCondition = std::bind(Environment<RT>::_apply_neumann_2d, this->densities);
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
    const unsigned Lx = this->densities.dims(0), nx = Lx + 1;
    const unsigned Ly = this->densities.dims(1), ny = Ly + 1;

    unsigned io = (unsigned)floor(Lx  / 5.0f),
            jo = (unsigned)floor(Ly / 5.0f),
            k = 20;
    array x = tile(moddims(seq(nx),nx,1), 1,ny);
    array y = tile(moddims(seq(ny),1,ny), nx,1);

    // Initial condition
    this->densities = 0.01f * exp((-((x - io) * (x - io) + (y - jo) * (y - jo))) / (k * k));

    this->simulate(20.0);
}

template <typename RT>
void Environment<RT>::simulate(double advanceTime) {

    RT normalizer = max<RT>(this->densities);

    timer t = timer::start();
    unsigned iter = 0;
    while (progress(iter, t, advanceTime)) {
        this->simulateTimeStep();
        if(this->visualizationWin != NULL){
            this->visualizationWin->image(this->densities/normalizer);
        }
        iter++;
    }
}

template <typename RT>
void Environment<RT>::simulateTimeStep() {
    this->applyBoundaryCondition();
    this->density_changes = convolve(this->densities, this->diffusion_filters);
    this->densities += this->density_changes*this->dt;
}

template <typename RT>
void Environment<RT>::_apply_neumann_2d(array input) {
    // X direction
    input(0, span, span) = input(1, span, span);
    input(end, span, span) = input(end-1, span, span);

    // Y direction
    input(span, 0, span) = input(span, 1, span);
    input(span, end, span) = input(span, end-1, span);

}


template class Environment<float>;
template class Environment<double>;