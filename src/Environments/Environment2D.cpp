//
// Created by Max Horn on 19/04/16.
//

#include "Environment2D.h"
#include <ios>


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

array Environment2D::getInterpolatedPositions(array &xpos, array &ypos) {
    array xindex = xpos/this->resolution + BORDER_SIZE;
    array yindex = ypos/this->resolution + BORDER_SIZE;

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

void Environment2D::setupStorage(unique_ptr<H5::Group> storage)  {
    // Store reference to group
    this->storage = std::move(storage);

    // Setup dimensions of datasets
    // TODO: this is the size in units and not the true dimension
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

    // Define datatypes and dataspaces for strings
    H5::StrType varstrtype(0, H5T_VARIABLE);
    H5::DataSpace scalar(H5S_SCALAR);
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

    hsize_t current_size[3], new_size[3];
    // Take a probe, aqnd calculate the new size after adding data
    this->ligands_storage.begin()->second->getSpace().getSimpleExtentDims(current_size);
    new_size[0] = current_size[0] + 1;
    new_size[1] = current_size[1];
    new_size[2] = current_size[2];

    // Parameters for hyperslap
    hsize_t count[3] = {1, new_size[1], new_size[2]};
    hsize_t start[3] = {current_size[0], 0, 0};

    // TODO: Maybe replace this with one copy operation of complete array and then writing via hyperslap
    for(auto ligand: this->ligands) {
        // Copy data from GPU
        GPU_REALTYPE *data = this->getDensity(ligand.ligandId).host<GPU_REALTYPE>();

        // Extend storage space
        this->ligands_storage[ligand.ligandId]->extend(new_size);
        // Define which data to write to where
        H5::DataSpace targetspace = this->ligands_storage[ligand.ligandId]->getSpace();
        targetspace.selectHyperslab(H5S_SELECT_SET, count, start);
        hsize_t dataLength = new_size[1]*new_size[2];

        H5::DataSpace sourceSpace(1, &dataLength);
        this->ligands_storage[ligand.ligandId]->write(data, HDF5_GPUTYPE, sourceSpace, targetspace);
        // Free allocated memory
        af::freeHost(data);
    }
}

void Environment2D::closeStorage() {
    // Calls destructor which also calls close
    this->ligands_storage.clear();
    this->storage.reset();
}

Environment2D::Environment2D(H5::Group group) : Environment(group) {
    init();
    int nLigands = group.getNumObjs();
    for(int i = 0; i < nLigands; i++) {
        // TODO
    }
}

void Environment2D::init() {
    densities = array(internal_dimensions[0], internal_dimensions[1], internal_dimensions[2], AF_GPUTYPE);
    diffusion_filters = constant(0.0, LAPLACIAN_SIZE, LAPLACIAN_SIZE, (dim_t)settings.ligands.size());

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


array Environment2D::Diffusion2D::rateofchange(array &input) {
    array changes = constant(0, parent->densities.dims());
    for (size_t i = 0; i < parent->densities.dims(2); i++) {
        changes(seq(1, end-1), seq(1, end-1), i) = convolve(input(span, span, i),
                                                            parent->diffusion_filters(span, span, i))(seq(1, end-1), seq(1, end-1));

        changes(seq(1, end-1), seq(1, end-1), i) += parent->ligands[i].globalProductionRate
                                                    - parent->ligands[i].globalDegradationRate*input(seq(1, end-1), seq(1, end-1), i);
    }
    changes.eval();
    return changes;
}