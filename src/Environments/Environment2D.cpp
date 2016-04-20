//
// Created by Max Horn on 19/04/16.
//

#include "Environment2D.h"

array Environment2D::getLaplacian() {
    GPU_REALTYPE data2 [] =
            {0.0, 1.0, 0.0,
             1.0, -4.0, 1.0,
             0.0, 1.0, 0.0};
    return array(3, 3, data2);
}

void Environment2D::applyNeumannBC(array input, double resolution, BoundaryCondition bc) {
    // X direction
    input(0, span, span) = input(1, span, span) - resolution*bc.xneg;
    input(end, span, span) = input(end-1, span, span) - resolution*bc.xpos;

    // Y direction
    input(span, 0, span) = input(span, 1, span) - resolution*bc.yneg;
    input(span, end, span) = input(span, end-1, span) - resolution*bc.ypos;
}

Environment2D::Environment2D(EnvironmentSettings settings) : Environment(settings) {
    densities = array(internal_dimensions[0], internal_dimensions[1], internal_dimensions[2], settings.dataType);
    density_changes =constant(0.0, internal_dimensions[0], internal_dimensions[1], internal_dimensions[2], settings.dataType);
    diffusion_filters = constant(0.0, LAPLACIAN_SIZE, LAPLACIAN_SIZE, (dim_t)settings.ligands.size());
    for(size_t i = 0; i < settings.ligands.size(); i++) {
        densities(span, span, i) = settings.ligands[i].initialConcentration;
        diffusion_filters(span, span, i) = Environment2D::getLaplacian();
        diffusion_filters(span, span, i) *= settings.ligands[i].diffusionCoefficient;
    }

    switch(this->boundaryCondition.type) {
        case BC_NEUMANN:
        default:
            this->applyBoundaryCondition = std::bind(Environment2D::applyNeumannBC, this->densities, this->resolution, this->boundaryCondition);
    }
}

void Environment2D::simulateTimeStep() {
    this->applyBoundaryCondition();
    this->density_changes = convolve(this->densities, this->diffusion_filters);
    this->densities += this->density_changes*this->dt;
    eval(this->density_changes, this->densities);
}

array Environment2D::getAllDensities() {
    return this->densities(seq(BORDER_SIZE,end-BORDER_SIZE), seq(BORDER_SIZE,end-BORDER_SIZE), span);
}

void Environment2D::test() {
    const unsigned Lx = this->densities.dims(0)-1, nx = Lx + 1;
    const unsigned Ly = this->densities.dims(1)-1, ny = Ly + 1;

    unsigned io = (unsigned)floor(Lx  / 5.0f),
            jo = (unsigned)floor(Ly / 5.0f),
            k = 20;
    array x = tile(moddims(seq(nx),nx,1), 1,ny);
    array y = tile(moddims(seq(ny),1,ny), nx,1);

    // Initial condition
    for (size_t i = 0; i < this->ligands.size(); i++){
        this->densities(span, span, i) = 20.0f * exp((-((x - io) * (x - io) + (y - jo) * (y - jo))) / (k * k));
    }
    this->simulate(60.0);
}

dim4 Environment2D::getSize() {
    dim4 dims = this->densities.dims();
    dims[0] = dims[0] - 2*BORDER_SIZE;
    dims[1] = dims[1] - 2*BORDER_SIZE;
    dims[2] = 1;
    return dims;
}

array Environment2D::getDensity(unsigned int ligand) {
    if(ligand < this->ligands.size())
        return this->densities(seq(1,end-1), seq(1,end-1), ligand);
    else
        throw new exception("Ligand index out of range");
}



