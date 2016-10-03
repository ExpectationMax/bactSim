//
// Created by Max Horn on 10/08/16.
//

#include "SimplePopulation.h"
#include "General/StorageHelper.h"

SimplePopulation::SimplePopulation(std::string name, shared_ptr<Environment> env, SimplePopulationParameters params) : BacterialPopulation(params), env(env), params(params) {
    this->name = name;
    init();
}

void SimplePopulation::init() {
    std::vector<int> ligandIds;

    this->uptakeRates  = array(1, (dim_t)params.interactions.size());
    this->productionRates = array(1, (dim_t)params.interactions.size());

    for(size_t i = 0; i < params.interactions.size(); i++) {
        ligandIds.push_back(params.interactions[i].ligandId);
        uptakeRates(0, i) = params.interactions[i].uptakeRate;
        productionRates(0, i) = params.interactions[i].productionRate;
    }

    ligandmapping = env->getLigandMapping(ligandIds);
    std::vector<double> size = env->getSize();
    maxx = size[0];
    maxy = size[1];

    switch(env->getBoundaryConditionType()) {
        case BC_PERIODIC:
            validatePositions = std::bind(SimplePopulation::applyPeriodicBoundary, maxx, maxy, std::ref(xpos), std::ref(ypos));
            spaciallyLimitedEnv = false;
            break;
        default:
            validatePositions = std::bind(SimplePopulation::applySolidBoundary, maxx, maxy, std::ref(xpos), std::ref(ypos), std::ref(atborder));
            spaciallyLimitedEnv = true;
    }
}

void SimplePopulation::interactWithEnv(int individual, double dt) {
    interactWithEnvPos(interpolatedPositions(individual, span), weights(individual, span), individual, dt);
}

// This function could be used to interact with the environment if non-overlaping bacteria
// (no interactions with same grid points) are passed in the individuals array
void SimplePopulation::interactWithEnv(array individuals, double dt) {
    interactWithEnvPos(interpolatedPositions(individuals, span), weights(individuals, span), individuals, dt);
}

// this includes multiple integration steps, keep for future reference
//void SimplePopulation::interactWithEnvPos(array pos, int individual) {
//    array ligconcentrations = env->getLigandConcentrations(pos, ligandmapping);
//    array concentrationChange = constant(0, ligconcentrations.dims());
//    for(int i = 0; i < params.integrationMultiplyer; i++) {
//        array change = (-ligconcentrations*uptakeRates + productionRates)*(params.dt/params.integrationMultiplyer);
//        concentrationChange += change;
//        ligconcentrations += change;
//    }
////    af_print(concentrationChange);
////    af_print(ligconcentrations);
//    concentrations(individual, span) = ligconcentrations;
//    env->changeLigandConcentrationBy(concentrationChange, pos, ligandmapping);
//}

void SimplePopulation::interactWithEnvPos(array pos, array w, int individual, double dt) {
    array ligconcentrations = env->getLigandConcentrations(pos, w, ligandmapping);
    array concentrationChange = modelUptakeProductionRate(ligconcentrations);
    concentrations(individual, span) = ligconcentrations+concentrationChange*dt;
    env->changeLigandConcentrationBy(concentrationChange, pos, w, ligandmapping);
}

void SimplePopulation::interactWithEnvPos(array pos, array w, array individuals, double dt) {
    array ligconcentrations = env->getLigandConcentrations(pos, w, ligandmapping);
    array concentrationChange = modelUptakeProductionRate(ligconcentrations);
    concentrations(individuals, span) = ligconcentrations+concentrationChange*dt;
    env->changeLigandConcentrationBy(concentrationChange, pos, w, ligandmapping);
}

void SimplePopulation::applyPeriodicBoundary(double maxx, double maxy, array &xpos, array &ypos) {
    // x axis
    xpos += (xpos > maxx) * -xpos + (xpos < 0) * (-xpos + maxx);
    // y axis
    ypos = ypos - (ypos > maxy) * ypos + (ypos < 0) * (-ypos + maxy);

    eval(xpos, ypos);
}

void SimplePopulation::applySolidBoundary(double maxx, double maxy, array &xpos, array &ypos, array &atborder) {
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

    // Mark bacteria as at border
    atborder =  max(outofrange, 1);
    eval(xpos, ypos);
    eval(atborder);
}

void SimplePopulation::liveTimestep(double dt) {
    senseLigandConcentration();
    simulate(dt);
    move(dt);
    validatePositions();
    updateInterpolatedPositions();
}

void SimplePopulation::simulate(double dt) {
    concentrations = env->getLigandConcentrations(interpolatedPositions, weights, ligandmapping);
}

void SimplePopulation::updateInterpolatedPositions() {
    env->setInterpolatedPositions(xpos, ypos, interpolatedPositions, weights);
}

void SimplePopulation::move(double dt) {
    xpos += cos(angle)*params.swimmSpeed*dt;
    ypos += sin(angle)*params.swimmSpeed*dt;
}

void SimplePopulation::setupStorage(H5::Group storage) {
    // Name and type are initialised by base class
    BacterialPopulation::setupStorage(storage);

    H5::Attribute swimmSpeed = this->storage->createAttribute("Swimm speed", H5::PredType::IEEE_F64LE, StorageHelper::H5Scalar);
    swimmSpeed.write(HDF5_GPUTYPE, &this->params.swimmSpeed);

    // Store interactions
    hsize_t interCount = this->params.interactions.size();
    H5::DataSpace interSpace(1, &interCount);
    H5::CompType interType = LigandInteraction::getH5SaveType();
    H5::Attribute interactions = this->storage->createAttribute("Ligand interactions", interType, interSpace);
    interactions.write(interType, this->params.interactions.data());

    // Initialize DataSets for bacterial parameters
    hsize_t bactCount = this->size;
    hsize_t initDims[2] = {0, bactCount};
    hsize_t maxDims[2] = {H5S_UNLIMITED, bactCount};
    H5::DataSpace bactSpace(2, initDims, maxDims);
    this->storageSpace = bactSpace;
    H5::DSetCreatPropList properties(H5::DSetCreatPropList::DEFAULT);
    hsize_t chunk_dims[2] = {4, bactCount};
    properties.setChunk(2, chunk_dims);
    this->storageProperties = properties;

    // Store as 64 bit double independent of architecture
    this->xposStorage.reset(
            new H5::DataSet(this->storage->createDataSet("xpos", H5::PredType::IEEE_F64LE, this->storageSpace, this->storageProperties)));
    this->yposStorage.reset(
            new H5::DataSet(this->storage->createDataSet("ypos", H5::PredType::IEEE_F64LE, this->storageSpace, this->storageProperties)));
    this->angleStorage.reset(
            new H5::DataSet(this->storage->createDataSet("angle", H5::PredType::IEEE_F64LE, this->storageSpace, this->storageProperties)));
}

bool SimplePopulation::save() {
    if (!this->storage)
        return false;

    StorageHelper::appendDataToDataSet<GPU_REALTYPE>(xpos, *xposStorage, HDF5_GPUTYPE);
    StorageHelper::appendDataToDataSet<GPU_REALTYPE>(ypos, *yposStorage, HDF5_GPUTYPE);
    StorageHelper::appendDataToDataSet<GPU_REALTYPE>(angle, *angleStorage, HDF5_GPUTYPE);
    return true;
}

void SimplePopulation::closeStorage() {
    // Call to reset also calls destructor
    this->xposStorage.reset();
    this->yposStorage.reset();
    this->angleStorage.reset();

    // Finally close group
    this->storage.reset();
}

void SimplePopulation::randomizeAngle() {
    angle = (2 * af::Pi * randu(size)).as(AF_GPUTYPE);
}

void SimplePopulation::setPositions(array x, array y) {
    xpos = x.as(AF_GPUTYPE);
    ypos = y.as(AF_GPUTYPE);
    validatePositions();
    updateInterpolatedPositions();
}

SimplePopulation::SimplePopulation(std::string name, shared_ptr<Environment> Env, SimplePopulationParameters parameters,
                                     int nBacteria) : SimplePopulation(name, Env, parameters) {
    size = nBacteria;
    initializeArrays();
    randomizeAngle();

    array randx = randu(size) * maxx;
    array randy = randu(size) * maxy;
    setPositions(randx, randy);
    validatePositions();
    updateInterpolatedPositions();
    // Just get ligand concentrations
    simulate(0);
    senseLigandConcentration();
}

SimplePopulation::SimplePopulation(std::string name, shared_ptr<Environment> Env, SimplePopulationParameters parameters,
                                     int nBacteria, GPU_REALTYPE *initialx, GPU_REALTYPE *initialy) :
        SimplePopulation(name, Env, parameters) {
    size = nBacteria;
    initializeArrays();
    randomizeAngle();
    setPositions(array(size, initialx), array(size, initialy));
    validatePositions();
    updateInterpolatedPositions();
    // Just get ligand concentrations
    simulate(0);
    senseLigandConcentration();
}

SimplePopulation::SimplePopulation(shared_ptr<Environment> Env, H5::Group group) : BacterialPopulation(group) {
    // Name and storage initialized by base class
    SimplePopulationParameters parameters;

    // read parameters
    group.openAttribute("Swimm speed").read(HDF5_GPUTYPE, &parameters.swimmSpeed);

    // read ligand interactions
    H5::Attribute ligInteractions = group.openAttribute("Ligand interactions");
    hsize_t nInteractions;
    H5::DataSpace interationSpace = ligInteractions.getSpace();
    interationSpace.getSimpleExtentDims(&nInteractions);
    parameters.interactions.resize(nInteractions);
    ligInteractions.read(LigandInteraction::getH5ReadType(), parameters.interactions.data());

    this->env = Env;
    this->params = parameters;

    H5::DataSet xpos = group.openDataSet("xpos");
    hsize_t bactDims[2];
    xpos.getSpace().getSimpleExtentDims(bactDims);
    this->size = bactDims[1];
    this->init();
    initializeArrays();

    this->xpos = StorageHelper::loadLastDataToGpu<GPU_REALTYPE>(xpos, HDF5_GPUTYPE, AF_GPUTYPE);
    this->xposStorage.reset(new DataSet(xpos));

    H5::DataSet ypos = group.openDataSet("ypos");
    this->ypos = StorageHelper::loadLastDataToGpu<GPU_REALTYPE>(ypos, HDF5_GPUTYPE, AF_GPUTYPE);
    this->yposStorage.reset(new DataSet(ypos));

    H5::DataSet angle = group.openDataSet("angle");
    this->angle = StorageHelper::loadLastDataToGpu<GPU_REALTYPE>(angle, HDF5_GPUTYPE, AF_GPUTYPE);
    this->angleStorage.reset(new DataSet(angle));

    validatePositions();
    updateInterpolatedPositions();
    simulate(0);
    senseLigandConcentration();
}

REGISTER_DEF_TYPE(SimplePopulation);

void SimplePopulation::printInternals() {
    af_print(xpos);
    af_print(ypos);
    af_print(concentrations);
    af_print(angle);
    af_print(interpolatedPositions);
    af_print(uptakeRates);
    af_print(productionRates);
    af_print(ligandmapping);
}

void SimplePopulation::senseLigandConcentration() {
    sensedConcentration = sum(concentrations, 1);
    eval(sensedConcentration);
}

array SimplePopulation::modelUptakeProductionRate(array ligconc) {
    return (-ligconc*tile(uptakeRates, ligconc.dims(0)) + tile(productionRates, ligconc.dims(0)));
}

array SimplePopulation::getInterpolatedPositions() {
    return interpolatedPositions;
}

void SimplePopulation::initializeArrays() {
    concentrations = constant(0, size, params.interactions.size());
    interpolatedPositions = array(size, 4,  af::dtype::u32);
    weights = array(size, 4, AF_GPUTYPE);
}


