//
// Created by Max Horn on 10/08/16.
//

#include "ExamplePopulation.h"

ExamplePopulation::ExamplePopulation(shared_ptr<Environment2D> env, BacterialParameters params) : env(env), params(params) {
    std::vector<int> ligandIds;
    Koff = array((dim_t)params.interactions.size());
    Kon = array((dim_t)params.interactions.size());
    uptakeRates  = array((dim_t)params.interactions.size());
    productionRates = array((dim_t)params.interactions.size());

    for(size_t i = 0; i < params.interactions.size(); i++) {
        ligandIds.push_back(params.interactions[i].ligandId);
        Koff(i) = params.interactions[i].Koff;
        Kon(i) = params.interactions[i].Kon;
        uptakeRates(i) = params.interactions[i].uptakeRate;
        productionRates(i) = params.interactions[i].productionRate;
    }
    ligandmapping = env->getLigandMapping(ligandIds);
    std::vector<double> size = env->getSize();
    maxx = size[0];
    maxy = size[1];

    switch(env->getBoundaryConditionType()) {
        case BC_PERIODIC:
            validatePositions = std::bind(ExamplePopulation::applyPeriodicBoundary, maxx, maxy, std::ref(xpos), std::ref(ypos));
            break;
        default:
            validatePositions = std::bind(ExamplePopulation::applySolidBoundary, maxx, maxy, std::ref(xpos), std::ref(ypos), std::ref(tumbling));
    }
}

void ExamplePopulation::interactWithEnv(int individual) {
    array pos = interpolatedPositions(span, individual);
    interactWithEnvPos(pos);
}

void ExamplePopulation::interactWithEnv(array individuals) {
    array pos = interpolatedPositions(span, individuals);
    interactWithEnvPos(pos);
}

void ExamplePopulation::interactWithEnvPos(array pos) {
    array concentrations = env->getLigandConcentrations(pos, ligandmapping);
    array concentrationChange = -concentrations*uptakeRates*params.dt + productionRates*params.dt;
    env->changeLigandConcentrationBy(concentrationChange, pos, ligandmapping);
}


void ExamplePopulation::applyPeriodicBoundary(int maxx, int maxy, array &xpos, array &ypos) {
    // x axis
    xpos -= (xpos > maxx) * xpos;

    xpos += (xpos < 0) * (-xpos + maxx);

    // y axis
    ypos -= (ypos > maxy) * ypos;

    ypos += (ypos < 0) * (-ypos + maxy);
};

void ExamplePopulation::applySolidBoundary(int maxx, int maxy, array &xpos, array &ypos, array &tumbling) {
    // Shift bacteria that are out of range back into field
    // x axis
    array outofrange = array(xpos.dims(0), 4, b8);
    outofrange(span, 0) = xpos > maxx;
    xpos -= outofrange(span, 0) * (xpos-maxx);

    outofrange(span, 1) = xpos < 0;
    xpos -= outofrange(span, 1) * xpos;

    // y axis

    outofrange(span, 2) = ypos > maxy;
    ypos -= outofrange(span, 2) * (ypos-maxy);

    outofrange(span, 3) = ypos < 0;
    ypos -= outofrange(span, 3) * ypos;

    // set tumbling to true for out of range bacteria
    tumbling = tumbling || max(outofrange, 1);
};

void ExamplePopulation::liveTimestep() {
    simulate();
    move();
    validatePositions();
    updateInterpolatedPositions();
}

void ExamplePopulation::simulate() {
    array concentrations = env->getLigandConcentrations(interpolatedPositions, ligandmapping);
}

void ExamplePopulation::updateInterpolatedPositions() {
    interpolatedPositions = env->getInterpolatedPositions(xpos, ypos);
}


void ExamplePopulation::move() {
    xpos += !tumbling*af::cos(angle)*params.swimmSpeed*params.dt;
    ypos += !tumbling*af::sin(angle)*params.swimmSpeed*params.dt;
}

void ExamplePopulation::setupStorage(shared_ptr<H5::Group> storage) {
    // Create group with name of bacterial population
    this->storage.reset(new H5::Group(storage->createGroup(this->params.name)));

    // Store parameters of population
    H5::DataSpace scalar(H5S_SCALAR);
    H5::Attribute swimmSpeed = this->storage->createAttribute("Swimm speed", H5::PredType::IEEE_F64LE, scalar);
    swimmSpeed.write(HDF5_GPUTYPE, &this->params.swimmSpeed);

    H5::Attribute dt = this->storage->createAttribute("dt", H5::PredType::IEEE_F64LE, scalar);
    dt.write(HDF5_GPUTYPE, &this->params.dt);

    // Store interactions
    hsize_t interCount = this->params.interactions.size();
    H5::DataSpace interSpace(1, &interCount);
    H5::CompType interType = LigandInteraction::getH5type();
    H5::Attribute interactions = this->storage->createAttribute("Ligand interactions", interType, interSpace);
    interactions.write(interType, this->params.interactions.data());

    // Initialize DataSets for bacterial parameters
    hsize_t bactCount = this->size;
    hsize_t initDims[2] = {0, bactCount};
    hsize_t maxDims[2] = {H5S_UNLIMITED, bactCount};
    H5::DataSpace bactSpace(2, initDims, maxDims);
    H5::DSetCreatPropList properties(H5::DSetCreatPropList::DEFAULT);
    hsize_t chunk_dims[2] = {4, bactCount};
    properties.setChunk(2, chunk_dims);

    // Store as 64 bit double independent of architecture
    this->xposStorage.reset(
            new H5::DataSet(this->storage->createDataSet("xpos", H5::PredType::IEEE_F64LE, bactSpace, properties)));
    this->yposStorage.reset(
            new H5::DataSet(this->storage->createDataSet("ypos", H5::PredType::IEEE_F64LE, bactSpace, properties)));
    this->angleStorage.reset(
            new H5::DataSet(this->storage->createDataSet("angle", H5::PredType::IEEE_F64LE, bactSpace, properties)));
    this->tumblingStorage.reset(
            new H5::DataSet(this->storage->createDataSet("tumbling", H5::PredType::STD_I8LE, bactSpace, properties)));
}

void ExamplePopulation::save() {
    if(!this->storage)
        return;

    H5::DataSpace targetSpace;


    hsize_t current_size[2], new_size[2];
    // Take a probe, and calculate the new size after adding data, all dimensions are the same
    this->xposStorage->getSpace().getSimpleExtentDims(current_size);
    new_size[0] = current_size[0] + 1;
    new_size[1] = current_size[1];

    // Parameters for hyperslap
    hsize_t count[2] = {1, current_size[1]};
    hsize_t start[2] = {current_size[0], 0};

    // infer source dimensions
    H5::DataSpace sourceSpace(1, &current_size[1]);

    // Save xpos
    GPU_REALTYPE *data = this->xpos.host<GPU_REALTYPE>();
    this->xposStorage->extend(new_size);
    targetSpace = this->xposStorage->getSpace();
    targetSpace.selectHyperslab(H5S_SELECT_SET, count, start);
    this->xposStorage->write(data, HDF5_GPUTYPE, sourceSpace, targetSpace);
    af::freeHost(data);

    // Save ypos
    data = this->ypos.host<GPU_REALTYPE>();
    this->yposStorage->extend(new_size);
    targetSpace = this->yposStorage->getSpace();
    targetSpace.selectHyperslab(H5S_SELECT_SET, count, start);
    this->yposStorage->write(data, HDF5_GPUTYPE, sourceSpace, targetSpace);
    af::freeHost(data);

    // Save angle
    data = this->angle.host<GPU_REALTYPE>();
    this->angleStorage->extend(new_size);
    targetSpace = this->angleStorage->getSpace();
    targetSpace.selectHyperslab(H5S_SELECT_SET, count, start);
    this->angleStorage->write(data, HDF5_GPUTYPE, sourceSpace, targetSpace);
    af::freeHost(data);

    // Save tumbling
    char *bdata = this->tumbling.host<char>();
    this->tumblingStorage->extend(new_size);
    targetSpace = this->tumblingStorage->getSpace();
    targetSpace.selectHyperslab(H5S_SELECT_SET, count, start);
    this->tumblingStorage->write(bdata, H5::PredType::NATIVE_CHAR, sourceSpace, targetSpace);
    af::freeHost(bdata);
}


void ExamplePopulation::closeStorage() {
    // Call to reset also calls destructor
    this->xposStorage.reset();
    this->yposStorage.reset();
    this->angleStorage.reset();
    this->tumblingStorage.reset();

    // Finally close group
    this->storage.reset();
}

void ExamplePopulation::randomizeAngleAndTumbiling() {
    angle = 2 * af::Pi * randu(size);
    tumbling = (randu(size)>= 0.5).as(b8);
}

void ExamplePopulation::setPositions(array x, array y) {
    xpos = x;
    ypos = y;
    validatePositions();
    updateInterpolatedPositions();
}

ExamplePopulation::ExamplePopulation(shared_ptr<Environment2D> Env, BacterialParameters parameters,
                                     int nBacteria) :
        ExamplePopulation(Env, parameters) {
    size = nBacteria;

    randomizeAngleAndTumbiling();

    array randx = randu(size) * maxx;
    array randy = randu(size) * maxy;
    setPositions(randx, randy);
}

ExamplePopulation::ExamplePopulation(shared_ptr<Environment2D> Env, BacterialParameters parameters,
                                     int nBacteria, GPU_REALTYPE *initialx, GPU_REALTYPE *initialy) :
        ExamplePopulation(Env, parameters) {
    size = nBacteria;

    randomizeAngleAndTumbiling();

    setPositions(array(size, initialx), array(size, initialy));
}