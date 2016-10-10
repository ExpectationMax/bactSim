//
// Created by Max Horn on 05/09/16.
//

#ifndef BACTSIM_GPU_KOLLMANN2005POPULATION_H
#define BACTSIM_GPU_KOLLMANN2005POPULATION_H

#include "Environments/Environment.h"
#include "SimplePopulation.h"

using namespace af;

struct Matthaeus2009Parameters : SimplePopulationParameters {
    Matthaeus2009Parameters() : SimplePopulationParameters() {};
    Matthaeus2009Parameters(shared_ptr<Solver> odesolver, std::vector<LigandInteraction> interactions, GPU_REALTYPE swimmSpeed):
            SimplePopulationParameters(interactions, swimmSpeed), odesolver(odesolver) {};
    Matthaeus2009Parameters(SimplePopulationParameters paramsBase) {
        interactions = paramsBase.interactions;
        swimmSpeed = paramsBase.swimmSpeed;
    }

    shared_ptr<Solver> odesolver;

    // Warning, these following parameters are not actually saved
    unsigned int rezeptorMethylationLevels = 5;

    // Methylation dependent parameters
    std::vector<GPU_REALTYPE> T_Km {2.7, 20, 150, 1500, 60000}; // uM
    std::vector<GPU_REALTYPE> T_V {0, 0.25, 0.5, 0.75, 1};

    std::vector<GPU_REALTYPE> T {4.06, 1.08, 0.14, 0.013, 0.007};

    GPU_REALTYPE Tt = 5.3;
    GPU_REALTYPE T_H = 1.2;

    // Parameters for simulation of stochasticity
    GPU_REALTYPE beta = 0.0008;
    GPU_REALTYPE Rt_lower = 0;
    GPU_REALTYPE Rt_upper = 0.32;

    // MM constants
    GPU_REALTYPE K_R = 0.099; // uM
    GPU_REALTYPE K_B = 2.5; // uM
    GPU_REALTYPE K_C = 3.1; // uM

    // Hill coefficient
    GPU_REALTYPE H_c = 10.3;

    // initial concentrations
    GPU_REALTYPE A_t = 5.3; // uM
    GPU_REALTYPE Ap = 0.13; // uM
    GPU_REALTYPE B_t = 0.28; // uM
    GPU_REALTYPE Bp = 0.079; // uM

    GPU_REALTYPE R_t = 0.16; // uM
    GPU_REALTYPE Y_t = 9.7; // uM
    GPU_REALTYPE Yp = 2.7; // uM
    GPU_REALTYPE Z_t = 4.25; // uM

    GPU_REALTYPE pwDivider = (0.2 * 1.53);

    // Rate constants
    GPU_REALTYPE k_R = 0.39; // 1/s
    GPU_REALTYPE k_B = 6.3; // 1/s
    GPU_REALTYPE kp_B = 3; // 1/(uM*s)
    GPU_REALTYPE k_A = 50; // 1/(uM*s)
    GPU_REALTYPE k_Y = 100; // 1/(uM*s)
//    GPU_REALTYPE k_Z = 30/Z_t; // To be divided by [CheZ] 1/s * [CheZ]
    GPU_REALTYPE k_Z = 7.89; // To be divided by [CheZ] 1/s * [CheZ]
    GPU_REALTYPE g_B = 1; // 1/(uM*s)
    GPU_REALTYPE g_Y = 0.1; // 1/(uM*s)
};

class Matthaeus2009Population : public SimplePopulation {
public:
    Matthaeus2009Population(std::string name, shared_ptr<Environment> Env, Matthaeus2009Parameters parameters, int nBacteria);
    Matthaeus2009Population(std::string name, shared_ptr<Environment> env, Matthaeus2009Parameters params) : SimplePopulation(name, env, params), Matthaeus2009Parameters(params) {
        init();
    }

    Matthaeus2009Population(shared_ptr<Environment> Env, H5::Group group);

    // Simulation
    void liveTimestep(double dt) override;
    virtual double getStabledt() override {return params.integrationMultiplyer*0.002; };
    void printInternals() override;

    // Storage
    void setupStorage(H5::Group storage) override;
    bool save() override;
    void closeStorage() override;

    REGISTER_DEC_TYPE(Matthaeus2009Population);
protected:
    // Initialization
    void init();

    // Simulation

    void move(double dt) override;

    // Parameters
    Matthaeus2009Parameters params;

    // Additional Arrays
    array swimming;
    array Kon;
    array Koff;

    std::vector<array> Tm;
    std::vector<array> Tma;

    array Ta;
    array Tt;

    array Ap;
    array Yp;
    array Bp;

    array tau;

    // Values stored for accelleration
    array Ttdivider;
    array Tadivider;

    // Solver
    shared_ptr<Solver> odesolver;

    // Storage
    unique_ptr<H5::DataSet> swimmingStorage;
    std::vector<unique_ptr<H5::DataSet>> TmStorage;
    std::vector<unique_ptr<H5::DataSet>> TmaStorage;
    unique_ptr<H5::DataSet> tauStorage;
    unique_ptr<H5::DataSet> ApStorage;
    unique_ptr<H5::DataSet> YpStorage;
    unique_ptr<H5::DataSet> BpStorage;
    unique_ptr<H5::DataSet> concentrationStorage;

private:
    // Differential Equations
    class dTm : public DifferentialEquation {
        Matthaeus2009Population *p;
        int methylationLevel;
    public:
        dTm(Matthaeus2009Population *par, int methLevel): p(par), methylationLevel(methLevel) {}
        array rateofchange(array &input) override;
    };

    class dAp : public DifferentialEquation {
        Matthaeus2009Population *p;
    public:
        dAp(Matthaeus2009Population *par): p(par) {}
        array rateofchange(array &input) override;
    };

    class dYp : public DifferentialEquation {
        Matthaeus2009Population *p;
    public:
        dYp(Matthaeus2009Population *par): p(par) {}
        array rateofchange(array &input) override;
    };

    class dBp : public DifferentialEquation {
        Matthaeus2009Population *p;
    public:
        dBp(Matthaeus2009Population *par): p(par) {}
        array rateofchange(array &input) override;
    };

    // For simulation of stochasticity
    class dR : public DifferentialEquation {
        Matthaeus2009Population *p;
    public:
        dR(Matthaeus2009Population *par): p(par) {}
        array rateofchange(array &input) override;
    };

    std::vector<std::tuple<unique_ptr<DifferentialEquation>,array&>> equations;
    void updateSwimming(double dt);
    void calculateDividers();
    void integrateEquations(double dt);
    void updateTotalConc();
    void setBorderBacteriaTumbling();
};


#endif //BACTSIM_GPU_KOLLMANN2005_H
