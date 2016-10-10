//
// Created by Max Horn on 05/09/16.
//

#include "Matthaeus2009Population.h"
#include "General/StorageHelper.h"
#include <General/ArrayFireHelper.h>

/**
 * Extract parameters from parameter struct and store as GPU arrays for faster access
 */
void Matthaeus2009Population::init() {
    this->Koff = array((dim_t)params.interactions.size());
    this->Kon = array((dim_t)params.interactions.size());
    this->swimming = constant(0, size, af::dtype::b8);
    for(size_t i = 0; i < params.interactions.size(); i++) {
        Koff(i) = params.interactions[i].Koff;
        Kon(i) = params.interactions[i].Kon;
    }

    equations.reserve(8);
    Tm.reserve(5);
    TmStorage.resize(5);
    Tma.reserve(5);
    TmaStorage.resize(5);

    for(int i = 0; i<5; i++) {
        // Methylation levels of rez are 0
        Tm.push_back(constant(params.T[i], size, AF_GPUTYPE));
        Tma.push_back(constant(0, size, AF_GPUTYPE));
        equations.push_back(std::make_tuple(unique_ptr<DifferentialEquation>(new dTm(this, i)), std::ref(Tm[i])));
    }
    Ap = constant(params.Ap, size, AF_GPUTYPE);
    equations.push_back(std::make_tuple(unique_ptr<DifferentialEquation>(new dAp(this)), std::ref(Ap)));
    Bp = constant(params.Bp, size, AF_GPUTYPE);
    equations.push_back(std::make_tuple(unique_ptr<DifferentialEquation>(new dBp(this)), std::ref(Bp)));
    Yp = constant(params.Yp, size, AF_GPUTYPE);
    equations.push_back(std::make_tuple(unique_ptr<DifferentialEquation>(new dYp(this)), std::ref(Yp)));
    tau = constant(0, size, AF_GPUTYPE);
}

void Matthaeus2009Population::setupStorage(H5::Group mystorage) {
    SimplePopulation::setupStorage(mystorage);
    // Store type of ode solver
    storage->createAttribute("Solver", StorageHelper::H5VariableString, StorageHelper::H5Scalar)
            .write(StorageHelper::H5VariableString, this->odesolver->getType());
    // Store additional required fields
    swimmingStorage.reset(
            new H5::DataSet(this->storage->createDataSet("swimming", H5::PredType::STD_I8LE, this->storageSpace, this->storageProperties)));
    YpStorage.reset(
            new H5::DataSet(this->storage->createDataSet("Yp", H5::PredType::IEEE_F64LE, this->storageSpace, this->storageProperties)));
    ApStorage.reset(
            new H5::DataSet(this->storage->createDataSet("Ap", H5::PredType::IEEE_F64LE, this->storageSpace, this->storageProperties)));
    BpStorage.reset(
            new H5::DataSet(this->storage->createDataSet("Bp", H5::PredType::IEEE_F64LE, this->storageSpace, this->storageProperties)));
    tauStorage.reset(
            new H5::DataSet(this->storage->createDataSet("tau", H5::PredType::IEEE_F64LE, this->storageSpace, this->storageProperties)));
    concentrationStorage.reset(
            new H5::DataSet(this->storage->createDataSet("Concentrations", H5::PredType::IEEE_F64LE, this->storageSpace, this->storageProperties)));

    // Rezeptor methylation stages and activities
    for(auto i = 0; i < 5; i++) {
        std::ostringstream TmStream, TmaStream;
        TmStream << "Tm[" << i << "]";
        TmaStream << "Tma[" << i << "]";
        TmStorage[i].reset(
                new H5::DataSet(this->storage->createDataSet(TmStream.str(), H5::PredType::IEEE_F64LE, this->storageSpace, this->storageProperties)));
        TmaStorage[i].reset(
                new H5::DataSet(this->storage->createDataSet(TmaStream.str(), H5::PredType::IEEE_F64LE, this->storageSpace, this->storageProperties)));
    }
}

bool Matthaeus2009Population::save() {
    if(SimplePopulation::save()) {
        StorageHelper::appendDataToDataSet<char>(swimming, *swimmingStorage, H5::PredType::NATIVE_CHAR);
        StorageHelper::appendDataToDataSet<GPU_REALTYPE>(Ap, *ApStorage, HDF5_GPUTYPE);
        StorageHelper::appendDataToDataSet<GPU_REALTYPE>(Bp, *BpStorage, HDF5_GPUTYPE);
        StorageHelper::appendDataToDataSet<GPU_REALTYPE>(Yp, *YpStorage, HDF5_GPUTYPE);
        StorageHelper::appendDataToDataSet<GPU_REALTYPE>(tau, *tauStorage, HDF5_GPUTYPE);
        StorageHelper::appendDataToDataSet<GPU_REALTYPE>(sensedConcentration, *concentrationStorage, HDF5_GPUTYPE);

        for(auto i = 0; i < 5; i++) {
            StorageHelper::appendDataToDataSet<GPU_REALTYPE>(Tm[i], *TmStorage[i], HDF5_GPUTYPE);
            StorageHelper::appendDataToDataSet<GPU_REALTYPE>(Tma[i], *TmaStorage[i], HDF5_GPUTYPE);
        }
        return true;
    }
    else
        return false;
}

void Matthaeus2009Population::closeStorage() {
    swimmingStorage.reset();
    ApStorage.reset();
    BpStorage.reset();
    YpStorage.reset();
    tauStorage.reset();
    concentrationStorage.reset();
    TmStorage.clear();
    TmaStorage.clear();
    SimplePopulation::closeStorage();
}

Matthaeus2009Population::Matthaeus2009Population(shared_ptr<Environment> Env, H5::Group group) : SimplePopulation(Env,
                                                                                                                  group) {
    this->params = Matthaeus2009Parameters(SimplePopulation::params);
    init();
    std::string solverName;
    group.openAttribute("Solver").read(StorageHelper::H5VariableString, solverName);
    this->odesolver = SolverFactory::createInstance(solverName);

    H5::DataSet swimming = group.openDataSet("swimming");
    this->swimming = StorageHelper::loadLastDataToGpu<char>(swimming, H5::PredType::NATIVE_CHAR, af::dtype::b8);
    this->swimmingStorage.reset(new DataSet(swimming));

    H5::DataSet Ap = group.openDataSet("Ap");
    this->Ap = StorageHelper::loadLastDataToGpu<GPU_REALTYPE>(Ap, HDF5_GPUTYPE, AF_GPUTYPE);
    this->ApStorage.reset(new DataSet(Ap));

    H5::DataSet Bp = group.openDataSet("Bp");
    this->Bp = StorageHelper::loadLastDataToGpu<GPU_REALTYPE>(Bp, HDF5_GPUTYPE, AF_GPUTYPE);
    this->BpStorage.reset(new DataSet(Bp));
    
    H5::DataSet Yp = group.openDataSet("Yp");
    this->Yp = StorageHelper::loadLastDataToGpu<GPU_REALTYPE>(Yp, HDF5_GPUTYPE, AF_GPUTYPE);
    this->YpStorage.reset(new DataSet(Yp));

    H5::DataSet tau = group.openDataSet("tau");
    this->tau = StorageHelper::loadLastDataToGpu<GPU_REALTYPE>(tau, HDF5_GPUTYPE, AF_GPUTYPE);
    this->tauStorage.reset(new DataSet(tau));

    H5::DataSet conc = group.openDataSet("Concentrations");
    this->sensedConcentration = StorageHelper::loadLastDataToGpu<GPU_REALTYPE>(conc, HDF5_GPUTYPE, AF_GPUTYPE);
    this->concentrationStorage.reset(new DataSet(conc));

    for(auto i = 0; i < 5; i++) {
        std::ostringstream TmStream, TmaStream;
        TmStream << "Tm[" << i << "]";
        H5::DataSet tm = this->storage->openDataSet(TmStream.str());
        Tm[i] = StorageHelper::loadLastDataToGpu<GPU_REALTYPE>(tm, HDF5_GPUTYPE, AF_GPUTYPE);
        TmStorage[i].reset(new DataSet(tm));

        TmaStream << "Tma[" << i << "]";
        H5::DataSet tma = this->storage->openDataSet(TmaStream.str());
        Tma[i] = StorageHelper::loadLastDataToGpu<GPU_REALTYPE>(tma, HDF5_GPUTYPE, AF_GPUTYPE);
        TmaStorage[i].reset(new DataSet(tma));
    }
}

Matthaeus2009Population::Matthaeus2009Population(std::string name, shared_ptr<Environment> Env,
                                               Matthaeus2009Parameters parameters, int nBacteria) : SimplePopulation(
        name, Env, parameters, nBacteria), params(parameters) {
    this->odesolver = params.odesolver;
    init();
};

REGISTER_DEF_TYPE(Matthaeus2009Population)

void Matthaeus2009Population::liveTimestep(double dt) {
    // Simulation
    senseLigandConcentration();
    updateTotalConc();
    calculateDividers();
    integrateEquations(dt);
//    af_print(Yp);
    // Movement
    updateSwimming(dt);
    move(dt);
    validatePositions();
    updateInterpolatedPositions();
    if(spaciallyLimitedEnv)
        setBorderBacteriaTumbling();
}

void Matthaeus2009Population::updateSwimming(double dt) {
    // Get new swimming candidates
    array subset = randu(size, AF_GPUTYPE) < (dt/params.pwDivider);
    tau = pow(Yp, params.H_c)/(pow(Yp, params.H_c) + pow(params.K_C, params.H_c));
//    tau = exp(params.H_c * log(Yp)) / ( exp(params.H_c * log(Yp)) + exp(params.H_c * log(params.K_C)) );

    array newswimming = !(randu(size, AF_GPUTYPE) < tau);
    // Change angle of swimming bacteria, directly calculate the angle in rad, 18.32 and -4.6 are the published values
    angle = !subset*angle + subset*(swimming*angle + !swimming*(angle + ArrayFireHelper::gammaSampler(size, 4, 18.32/360*2*Pi, -4.6/360*2*Pi)));

    // Update swimming
    swimming = !subset*swimming + subset*newswimming;
    eval(tau, angle);
    eval(swimming);
}

void Matthaeus2009Population::move(double dt) {
    xpos += swimming*cos(angle)*params.swimmSpeed*dt;
    ypos += swimming*sin(angle)*params.swimmSpeed*dt;
    eval(xpos,ypos);
}

void Matthaeus2009Population::updateTotalConc() {
    array T_tot = constant(0, size);
    array T_a = constant(0, size);
    for(int i = 0; i < 5; i++) {
        T_tot += Tm[i];
        Tma[i] = Tm[i] * params.T_V[i] * (1 - pow(sensedConcentration, params.T_H)/(pow(sensedConcentration, params.T_H) + pow(params.T_Km[i], params.T_H)));
        T_a += Tma[i];
    }
    Tt = T_tot;
    Ta = T_a;
    eval(Tt, T_a);
}

void Matthaeus2009Population::integrateEquations(double dt) {
    for(int i = 0; i < params.integrationMultiplyer; i++) {
        for(auto j = 0; j < equations.size(); j++) {
            odesolver->solveStep(*std::get<0>(equations[j]), std::get<1>(equations[j]),
                                 dt / params.integrationMultiplyer);
        }
    }
}

void Matthaeus2009Population::calculateDividers() {
    Ttdivider = 1/(params.K_R + Tt);
    Tadivider = 1/(params.K_B + Ta);
    eval(Ttdivider, Tadivider);
}

void Matthaeus2009Population::printInternals() {
    SimplePopulation::printInternals();
    af_print(sensedConcentration);
    for(int i = 0; i < 5; i++){
        af_print(Tm[i]);
        af_print(Tma[i]);
    }
    af_print(Ta);
    af_print(Tt);
    af_print(Ap);
    af_print(Bp);
    af_print(Yp);
}

void Matthaeus2009Population::setBorderBacteriaTumbling() {
    swimming = swimming && !atborder;
}

array Matthaeus2009Population::dTm::rateofchange(array &input) {
    array output =  - p->params.k_R*p->params.R_t*input*p->Ttdivider
                    - p->params.k_B*p->Bp * p->Tma[methylationLevel]*p->Tadivider;
    if(methylationLevel>0)
        output += p->params.k_R*p->params.R_t*p->Tm[methylationLevel-1]*p->Ttdivider;
    if(methylationLevel<4)
        output += p->params.k_B*p->Bp*(p->Tma[methylationLevel+1]*p->Tadivider);

    return output;

}

array Matthaeus2009Population::dAp::rateofchange(array &input) {
    array output = + p->params.k_A*(p->params.A_t - input)*p->Ta
                   - p->params.k_Y*input*(p->params.Y_t - p->Yp)
                   - p->params.kp_B*input*(p->params.B_t - p->Bp);
    return output;
}

array Matthaeus2009Population::dYp::rateofchange(array &input) {
    array output = + p->params.k_Y*p->Ap*(p->params.Y_t - input) - input* (p->params.k_Z*p->params.Z_t + p->params.g_Y);
    return output;
}

array Matthaeus2009Population::dBp::rateofchange(array &input) {
    array output = + p->params.kp_B*p->Ap*(p->params.B_t - input) - p->params.g_B*input;
    return output;
}

array Matthaeus2009Population::dR::rateofchange(array &input) {
    return array();
}
