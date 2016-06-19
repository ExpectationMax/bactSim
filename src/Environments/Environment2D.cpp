//
// Created by Max Horn on 19/04/16.
//

#include "Environment2D.h"
#include <iostream>
#include <exception>

array Environment2D::getLaplacian() {
    GPU_REALTYPE data2 [] =
            {0.0, 1.0, 0.0,
             1.0, -4.0, 1.0,
             0.0, 1.0, 0.0};
    return array(3, 3, data2);
}

Environment2D::Environment2D(EnvironmentSettings settings) : Environment(settings), diffequation(Diffusion2D(this)) {
    densities = array(internal_dimensions[0], internal_dimensions[1], internal_dimensions[2], settings.dataType);
    diffusion_filters = constant(0.0, LAPLACIAN_SIZE, LAPLACIAN_SIZE, (dim_t)settings.ligands.size());

    for(size_t i = 0; i < ligands.size(); i++) {
        densities(span, span, i) = settings.ligands[i].initialConcentration;
        diffusion_filters(span, span, i) = Environment2D::getLaplacian();
        diffusion_filters(span, span, i) *= settings.ligands[i].diffusionCoefficient/pow(resolution, 2);
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
//    const unsigned Lx = this->densities.dims(0)-1, nx = Lx + 1;
//    const unsigned Ly = this->densities.dims(1)-1, ny = Ly + 1;
////
//    unsigned io = (unsigned)floor(Lx  / 5.0f),
//            jo = (unsigned)floor(Ly / 5.0f),
//            k = 20;
//    array x = tile(moddims(seq(nx),nx,1), 1,ny);
//    array y = tile(moddims(seq(ny),1,ny), nx,1);
//
//    // Initial condition
//    for (size_t i = 0; i < this->ligands.size(); i++){
//        this->densities(span, span, i) = 20.0f * exp((-((x - io) * (x - io) + (y - jo) * (y - jo))) / (k * k));
//    }

    GPU_REALTYPE changes[] = {50};

    array Vmax = array(1, changes);
    GPU_REALTYPE xpositions[] = {3.0, 4.0, 5.0};
    GPU_REALTYPE ypositions[] = {3.0, 4.0, 5.0};
    array xpos = array(3, xpositions);
    array ypos = array(3, ypositions);

    array interpolated = getInterpolatedPositions(xpos, ypos);
    int ligandids[] = {0};
    array ligandindexes = array(1, ligandids);
//    double normalizer = max<double>(this->densities);
    double normalizer = 10;
    for(int i =0; i < 20000; i++) {
        array concentrations = moddims(this->getLigandConcentrations(interpolated, ligandindexes), ligandindexes.dims(0));
        this->changeLigandConcentrationBy(dt*Vmax*concentrations/(concentrations+20.0), interpolated, ligandindexes);
        //af_print(getDensity(0));
        this->simulateTimeStep();
        std::cout << i*dt << std::endl;
    }
}

std::vector<double> Environment2D::getSize() {
    std::vector<double> size;

    dim4 dims = this->densities.dims();
    // x
    size.push_back((dims[1] - 3* BORDER_SIZE)*resolution);
    // y
    size.push_back((dims[0] - 3* BORDER_SIZE)*resolution);
    return size;
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

    // X direction
    input->operator()(span, 0, span) = input->operator()(span, 1, span)*-1.0 + 2.0*bc->xneg;
    input->operator()(span, end, span) = input->operator()(span, end-1, span)*-1.0 + 2.0*bc->xpos;
    input->eval();
}


void Environment2D::applyPeriodicBC(array *input) {

    // Y direction
    input->operator()(0, span, span) = input->operator()(end-1, span, span);
    input->operator()(end, span, span) = input->operator()(1, span, span);

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

array Environment2D::getDensity(int ligandId) {
    array pos = ligandMapping(span, LIGANDID) == ligandId;
    if (sum<int>(pos) == 0)
        throw exception("Could not find provided ligandId in Environment.");
    array index = ligandMapping(pos, LIGANDINTERNAL);
    return this->densities(seq(BORDER_SIZE, end-BORDER_SIZE), seq(BORDER_SIZE, end-BORDER_SIZE), index);
}

array Environment2D::getInterpolatedPositions(array &xpos, array &ypos) {
    array xindex = xpos/this->resolution + 1.5*BORDER_SIZE;
    array yindex = ypos/this->resolution + 1.5*BORDER_SIZE;

    xindex = moddims(xindex, 1, xindex.dims(0));
    yindex = moddims(yindex, 1, yindex.dims(0));
    //af_print(xindex);
    //af_print(yindex);
    array output = array(8, xpos.dims(0));
    output(POS_LEFT, span) = af::floor(xindex);  // Left
    output(POS_RIGHT, span) = af::ceil(xindex);   // Right
    output(POS_TOP, span) = af::floor(yindex);  // Top
    output(POS_BOTTOM, span) = af::ceil(yindex);   // Bottom
    //af_print(output(seq(POS_LEFT, POS_BOTTOM), span));

    output(W_TOPLEFT, span) = ((xindex - output(POS_LEFT, span)) * (yindex - output(POS_TOP, span)));
    output(W_TOPRIGHT, span) = ((output(POS_RIGHT, span) - xindex) * (yindex - output(POS_TOP, span)));
    output(W_BOTTOMLEFT, span) = ((xindex - output(POS_LEFT, span)) * (output(POS_BOTTOM, span) - yindex));
    output(W_BOTTOMRIGHT, span) = ((output(POS_RIGHT, span) - xindex) * (output(POS_BOTTOM, span) - yindex));
    //af_print(output(seq(W_TOPLEFT, W_BOTTOMRIGHT), span));
    double normalization = pow(this->resolution, -2);
    output(seq(W_TOPLEFT, W_BOTTOMRIGHT), span) *= normalization;
    //af_print(output);
    // eval(output);
    return output;
}

array Environment2D::getLigandConcentrations(array posAndWeights, array ligands) {

    // We need to tile the weight array along the ligand axis to allow element wise multiplication)
    dim4 tilingDims = {1, 1, 1, 1};
    tilingDims[2] = ligands.dims(0);
    array left = posAndWeights(POS_LEFT);
    array right = posAndWeights(POS_RIGHT);
    array top = posAndWeights(POS_TOP);
    array bottom = posAndWeights(POS_BOTTOM);
    
    array ligdensities =
            densities(top, left, ligands) * tile(posAndWeights(W_TOPLEFT), tilingDims) +
            densities(top, right, ligands) * tile(posAndWeights(W_TOPRIGHT), tilingDims) +
            densities(bottom, left, ligands) * tile(posAndWeights(W_BOTTOMLEFT), tilingDims) +
            densities(bottom, right, ligands) * tile(posAndWeights(W_BOTTOMRIGHT), tilingDims);

//    eval(ligdensities);
    return moddims(ligdensities, ligdensities.dims(2));
}

void Environment2D::changeLigandConcentrationBy(array concDifferences, array posAndWeights, array ligands) {

    if( concDifferences.dims(0) != ligands.dims(0))
        throw exception("The number of provided concentrations has to be equal to the number of ligands");
    dim_t concentrations = concDifferences.dims(0);
    dim4 targetdims = {1, 1, concentrations};
    array left = posAndWeights(POS_LEFT);
    array right = posAndWeights(POS_RIGHT);
    array top = posAndWeights(POS_TOP);
    array bottom = posAndWeights(POS_BOTTOM);
    // tile is required, to allow element wise calculation between two array datatypes (requires same size)
    // moddims changes the metadata to convert the result array of the calculation into a z vector (required as operation is in place)
    densities(top, left, ligands) += moddims(concDifferences*tile(posAndWeights(W_TOPLEFT), concentrations), targetdims);
    densities(top, right, ligands) += moddims(concDifferences*tile(posAndWeights(W_TOPRIGHT), concentrations), targetdims);
    densities(bottom, left, ligands) += moddims(concDifferences*tile(posAndWeights(W_BOTTOMLEFT), concentrations), targetdims);
    densities(bottom, right, ligands) += moddims(concDifferences*tile(posAndWeights(W_BOTTOMRIGHT), concentrations), targetdims);
    // eval(densities, left, right, top, bottom);
}

void Environment2D::evalDensities() {
    eval(densities);
}