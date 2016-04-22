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

Environment2D::Environment2D(EnvironmentSettings settings) : Environment(settings) {
    densities = array(internal_dimensions[0], internal_dimensions[1], internal_dimensions[2], settings.dataType);
    diffusion_filters = constant(0.0, LAPLACIAN_SIZE, LAPLACIAN_SIZE, (dim_t)settings.ligands.size());
    for(size_t i = 0; i < settings.ligands.size(); i++) {
        densities(span, span, i) = settings.ligands[i].initialConcentration;
        diffusion_filters(span, span, i) = Environment2D::getLaplacian();
        diffusion_filters(span, span, i) *= settings.ligands[i].diffusionCoefficient;
    }

    switch (settings.convolutionType) {
        default:
        case CT_SERIAL:
            //density_changes =constant(0.0, internal_dimensions[0], internal_dimensions[1], settings.dataType);
            this->calculateTimeStep = std::bind(Environment2D::serialCalculateTimeStep, &this->densities,
                                                &this->diffusion_filters, this->dt, &this->ligands);
            break;
        case CT_AFBATCH:
            //density_changes =constant(0.0, internal_dimensions[0], internal_dimensions[1], internal_dimensions[2], settings.dataType);
            this->calculateTimeStep = std::bind(Environment2D::batchCalculateTimeStep, &this->densities,
                                                &this->diffusion_filters, this->dt, &this->ligands);
            break;
    }

    switch(this->boundaryCondition.type) {
        case BC_NEUMANN:
        default:
            this->applyBoundaryCondition = std::bind(Environment2D::applyNeumannBC, &this->densities, this->resolution, &this->boundaryCondition);
            break;
        case BC_DIRICHELET:
            this->applyBoundaryCondition = std::bind(Environment2D::applyDericheletBC, &this->densities, &this->boundaryCondition);
            break;
        case BC_PERIODIC:
            this->applyBoundaryCondition = std::bind(Environment2D::applyPeriodicBC, &this->densities);
            break;
    }
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

void Environment2D::applyNeumannBC(array *input, double resolution, BoundaryCondition *bc) {
    // Y direction
    input->operator()(0, span, span) = input->operator()(1, span, span) - resolution*bc->yneg;
    input->operator()(end, span, span) = input->operator()(end-1, span, span) - resolution*bc->ypos;
    input->eval();
    // X direction
    input->operator()(span, 0, span) = input->operator()(span, 1, span) - resolution*bc->xneg;
    input->operator()(span, end, span) = input->operator()(span, end-1, span) - resolution*bc->xpos;
    input->eval();
}

void Environment2D::applyDericheletBC(array *input, BoundaryCondition *bc) {
    // Y direction
    input->operator()(0, span, span) = input->operator()(1, span, span)*-1.0 + 2.0*bc->yneg;
    input->operator()(end, span, span) = input->operator()(end-1, span, span)*-1.0 + 2.0*bc->ypos;
    input->eval();
    // X direction
    input->operator()(span, 0, span) = input->operator()(span, 1, span)*-1.0 + 2.0*bc->xneg;
    input->operator()(span, end, span) = input->operator()(span, end-1, span)*-1.0 + 2.0*bc->xpos;
    input->eval();
}


void Environment2D::applyPeriodicBC(array *input) {

    // Y direction
    input->operator()(0, span, span) = input->operator()(end-1, span, span);
    input->operator()(end, span, span) = input->operator()(1, span, span);
    input->eval();
    // X direction
    input->operator()(span, 0, span) = input->operator()(span, end-1, span);
    input->operator()(span, end, span) = input->operator()(span, 1, span);
    input->eval();
}

void Environment2D::serialCalculateTimeStep(array *densities, array *diffusionFilters, double dt, std::vector<Ligand> *ligands) {
    for (size_t i = 0; i < densities->dims(2); i++) {
        array densityChange = convolve((*densities)(span, span, i),
                                       (*diffusionFilters)(span, span, i))(seq(1, end-1), seq(1, end-1));

        densities->operator()(seq(1, end-1), seq(1, end-1), i) += (densityChange
                + (*ligands)[i].globalProductionRate
                - (*ligands)[i].globalDegradationRate*(*densities)(seq(1, end-1), seq(1, end-1), i))*dt;

        densityChange.eval();
        densities->eval();
    }
}

void Environment2D::batchCalculateTimeStep(array *densities, array *diffusionFilters, double dt, std::vector<Ligand> *ligands) {
    array densityChanges = convolve(*densities, *diffusionFilters)(seq(1,end-1), seq(1,end-1), span);
    densities->operator()(seq(1,end-1), seq(1,end-1), span) += densityChanges*dt;
    densityChanges.eval();
    densities->eval();
}







