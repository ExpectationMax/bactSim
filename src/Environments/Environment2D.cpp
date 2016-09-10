//
// Created by Max Horn on 19/04/16.
//

#include "Environment2D.h"
#include "General/StorageHelper.h"
#include "General/ArrayFireHelper.h"

array Environment2D::getLaplacian() {
    GPU_REALTYPE data2 [] =
            {0.0, 1.0, 0.0,
             1.0, -4.0, 1.0,
             0.0, 1.0, 0.0};
    return array(3, 3, data2);
}

Environment2D::Environment2D(EnvironmentSettings settings, shared_ptr<Solver> solver) : Environment(settings, solver) {
    init();
}

Environment2D::Environment2D(H5::Group group) : Environment(group) {
    init();
    hsize_t current_size[3], count[3], start[3];
    for(auto ligand: this->ligands) {
        H5::DataSet ligData = group.openDataSet(ligand.name);
        this->densities(seq(BORDER_SIZE, end-BORDER_SIZE), seq(BORDER_SIZE, end-BORDER_SIZE), this->hostLigandMapping[ligand.ligandId]) =
                StorageHelper::loadLastDataToGpu<GPU_REALTYPE>(ligData, HDF5_GPUTYPE, AF_GPUTYPE);

        this->ligands_storage[ligand.ligandId] = std::unique_ptr<H5::DataSet>(new H5::DataSet(ligData));
    }
}

void Environment2D::init() {
    densities = array(internal_dimensions[0], internal_dimensions[1], internal_dimensions[2], AF_GPUTYPE);
    diffusion_filters = constant(0.0, LAPLACIAN_SIZE, LAPLACIAN_SIZE, (dim_t)this->ligands.size());

    for(size_t i = 0; i < ligands.size(); i++) {
        densities(span, span, i) = this->ligands[i].initialConcentration;
        diffusion_filters(span, span, i) = Environment2D::getLaplacian();
        diffusion_filters(span, span, i) *= this->ligands[i].diffusionCoefficient/pow(resolution, 2);
    }
    diffusionEquation.reset(new Diffusion2D(this));

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

std::vector<double> Environment2D::getSize() {
    std::vector<double> size;

    dim4 dims = this->densities.dims();
    // x
    size.push_back((dims[1] - 2* BORDER_SIZE)*resolution);
    // y
    size.push_back((dims[0] - 2* BORDER_SIZE)*resolution);
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

array Environment2D::getDensity(int ligandId) {
    array pos = ligandMapping(span, LIGANDID) == ligandId;
    if (sum<int>(pos) == 0)
        throw exception("Could not find provided ligandId in Environment.");
    array index = ligandMapping(pos, LIGANDINTERNAL);
    return this->densities(seq(BORDER_SIZE, end-BORDER_SIZE), seq(BORDER_SIZE, end-BORDER_SIZE), index);
}

void Environment2D::setInterpolatedPositions(array &xpos, array &ypos, array &positions, array &weights) {
    array xindex = xpos/this->resolution + BORDER_SIZE;
    array yindex = ypos/this->resolution + BORDER_SIZE;
    array pos = array(xpos.elements(), 4);
//    xindex = moddims(xindex, 1, xindex.dims(0));
//    yindex = moddims(yindex, 1, yindex.dims(0));
    //af_print(xindex);
    //af_print(yindex);
//    array pos = array(xpos.dims(0), 4, af::dtype::u32);
    array left = af::floor(xindex);  // Left
    array right = af::ceil(xindex);   // Right
    array top = af::floor(yindex);  // Top
    array bottom = af::ceil(yindex);   // Bottom
    //af_print(pos(seq(POS_LEFT, POS_BOTTOM), span));
    
//    array weights = array(xpos.dims(0), 4, AF_GPUTYPE);
    weights(span, W_TOPLEFT) = (xindex - left) * (yindex - top);
    weights(span, W_TOPRIGHT) = (right - xindex) * (yindex - top);
    weights(span, W_BOTTOMLEFT) = (xindex - left) * (bottom - yindex);
    weights(span, W_BOTTOMRIGHT) = (right - xindex) * (bottom - yindex);
    
    positions(span, I_TOPLEFT) = ArrayFireHelper::coordinateIndexing(densities, top, left);
    positions(span, I_TOPRIGHT) = ArrayFireHelper::coordinateIndexing(densities, top, right);
    positions(span, I_BOTTOMLEFT) = ArrayFireHelper::coordinateIndexing(densities, bottom, left);
    positions(span, I_BOTTOMRIGHT) = ArrayFireHelper::coordinateIndexing(densities, bottom, right);

    eval(weights, positions);
}

array Environment2D::get_concentrations(array &indexes, array &ligands) {
    array index = ArrayFireHelper::indexZAxis(densities, indexes, ligands);
    return moddims(densities(index), indexes.dims(0), ligands.dims(0));
}

array Environment2D::getLigandConcentrations(array positions, array weights, array ligands) {
    // We need to tile the weight array along the ligand axis to allow element wise multiplication)
    dim_t nligands = ligands.dims(0);
    array topleft = positions(span, I_TOPLEFT);
    array topright = positions(span, I_TOPRIGHT);
    array bottomleft = positions(span, I_BOTTOMLEFT);
    array bottomright = positions(span, I_BOTTOMRIGHT);

    array ligdensities =
            get_concentrations(topleft, ligands) * tile(weights(span, W_TOPLEFT), 1, nligands) +
            get_concentrations(topright, ligands) * tile(weights(span, W_TOPRIGHT), 1, nligands) +
            get_concentrations(bottomleft, ligands) * tile(weights(span, W_BOTTOMLEFT), 1, nligands) +
            get_concentrations(bottomright, ligands) * tile(weights(span, W_BOTTOMRIGHT), 1, nligands);
    eval(ligdensities);
    return ligdensities;
}

void Environment2D::changeLigandConcentrationBy(array concDifferences, array positions, array weights, array ligands) {
    int nLigands = ligands.dims(0);
    if( concDifferences.dims(1) != nLigands)
        throw exception("The number of provided concentrations has to be equal to the number of ligands");

    array topleft = positions(span, I_TOPLEFT);
    array topright = positions(span, I_TOPRIGHT);
    array bottomleft = positions(span, I_BOTTOMLEFT);
    array bottomright = positions(span, I_BOTTOMRIGHT);

    array alltopleft = ArrayFireHelper::indexZAxis(densities, topleft, ligands);
    array alltopright = ArrayFireHelper::indexZAxis(densities, topright, ligands);
    array allbottomleft = ArrayFireHelper::indexZAxis(densities, bottomleft, ligands);
    array allbottomright = ArrayFireHelper::indexZAxis(densities, bottomright, ligands);
    densities(alltopleft)  += flat(concDifferences)*tile(weights(span, W_TOPLEFT), nLigands);
    densities(alltopright)  += flat(concDifferences)*tile(weights(span, W_TOPRIGHT), nLigands);
    densities(allbottomleft)  += flat(concDifferences)*tile(weights(span, W_BOTTOMLEFT), nLigands);
    densities(allbottomright)  += flat(concDifferences)*tile(weights(span, W_BOTTOMRIGHT), nLigands);
    eval(densities);
//    array tiled = concDifferences * reorder(tile(positions(W_TOPLEFT, span), concentrations), 1, 0);
//    af_print(tiled);
    // tile is required, to allow element wise calculation between two array datatypes (requires same size)
    // moddims changes the metadata to convert the result array of the calculation into a z vector (required as operation is in place)
//    densities(top, left, ligands) += moddims(concDifferences*tile(positions(W_TOPLEFT, span), concentrations), targetdims);
//    densities(top, right, ligands) += moddims(concDifferences*tile(positions(W_TOPRIGHT, span), concentrations), targetdims);
//    densities(bottom, left, ligands) += moddims(concDifferences*tile(positions(W_BOTTOMLEFT, span), concentrations), targetdims);
//    densities(bottom, right, ligands) += moddims(concDifferences*tile(positions(W_BOTTOMRIGHT, span), concentrations), targetdims);
////    eval(densities, left, right, top, bottom);
}

void Environment2D::evalDensities() {
    eval(densities);
}

void Environment2D::setupStorage(unique_ptr<H5::Group> storage)  {
    // Let parent init name, dt and boundary condition
    Environment::setupStorage(std::move(storage));

    // Define datatypes and dataspaces for strings
    H5::StrType varstrtype(0, H5T_VARIABLE);
    H5::DataSpace scalar(H5S_SCALAR);

    // Setup dimensions of datasets
    std::vector<hsize_t> dims;
    for(auto dim: this->internal_dimensions) {
        dims.push_back(static_cast<hsize_t>(dim)-2*BORDER_SIZE);
    }

    hsize_t initial_dims[3] = {0, dims[0], dims[1] };
    // Only time dimension (first) will be extended
    hsize_t max_dims[3] = {H5S_UNLIMITED, dims[0], dims[1] };
    H5::DataSpace dataSpace(3, initial_dims, max_dims);

    // Enable chunking for extendable storage, 4 time points into one chunk to decrease allocation / seeking time
    H5::DSetCreatPropList properties(H5::DSetCreatPropList::DEFAULT);
    hsize_t chunk_dims[3] = {4, dims[0], dims[1] };
    properties.setChunk(3, chunk_dims);


    for(auto ligand: this->ligands){
        H5::DataSet liganddataset = H5::DataSet(this->storage->createDataSet(ligand.name, H5::PredType::IEEE_F64LE, dataSpace, properties));
        liganddataset.createAttribute("Name", varstrtype, scalar).write(varstrtype, ligand.name);
        H5::Attribute properties = liganddataset.createAttribute("Properties", Ligand::getH5SaveType(), scalar);
        properties.write(Ligand::getH5ReadType(), &ligand);

        // Store a unique ptr to this dataset for fast storage later on
        this->ligands_storage[ligand.ligandId] = std::unique_ptr<H5::DataSet>(new H5::DataSet(liganddataset));
    }
}

void Environment2D::save() {
    if(!this->storage)
        return;

    // TODO: Maybe replace this with one copy operation of complete array and then writing via hyperslap
    for(auto ligand: this->ligands)
        StorageHelper::appendDataToDataSet<GPU_REALTYPE>(this->getDensity(ligand.ligandId), *this->ligands_storage[ligand.ligandId], HDF5_GPUTYPE);
}

void Environment2D::closeStorage() {
    // Calls destructor which also calls close
    this->ligands_storage.clear();
    this->storage.reset();
}

array Environment2D::Diffusion2D::rateofchange(array &input) {
    array changes = constant(0, parent->densities.dims());
    for (size_t i = 0; i < parent->densities.dims(2); i++) {
        changes(seq(1, end-1), seq(1, end-1), i) = convolve(input(span, span, i),
                                                            parent->diffusion_filters(span, span, i))(seq(1, end-1), seq(1, end-1));
                                                    + parent->ligands[i].globalProductionRate
                                                    - parent->ligands[i].globalDegradationRate*input(seq(1, end-1), seq(1, end-1), i);
    }
    changes.eval();
    return changes;
}