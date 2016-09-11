//
// Created by Max Horn on 05/09/16.
//

#ifndef BACTSIM_GPU_KOLLMANN2005POPULATION_H
#define BACTSIM_GPU_KOLLMANN2005POPULATION_H

#include "Environments/Environment2D.h"
#include "SimplePopulation.h"

using namespace af;

struct Kollmann2005Parameters : SimplePopulationParameters {
    Kollmann2005Parameters() : SimplePopulationParameters() {};
    Kollmann2005Parameters(shared_ptr<Solver> odesolver, std::vector<LigandInteraction> interactions, GPU_REALTYPE dt, GPU_REALTYPE swimmSpeed):
            SimplePopulationParameters(interactions, swimmSpeed), odesolver(odesolver) {};
    Kollmann2005Parameters(SimplePopulationParameters paramsBase) {
        this->interactions = paramsBase.interactions;
        this->swimmSpeed = paramsBase.swimmSpeed;
    }

    shared_ptr<Solver> odesolver;

    // Warning, these following parameters are not actually saved
    unsigned int integrationMultiplyer = 5;
    unsigned int rezeptorMethylationLevels = 5;

    // Methylation dependent parameters
    std::vector<GPU_REALTYPE> T_Km {2.7, 200, 1500, 15000, 60000}; // uM
    std::vector<GPU_REALTYPE> T_V {0, 0.25, 0.5, 0.75, 1};
    GPU_REALTYPE T_H = 1.2;

    // Rate constants
    GPU_REALTYPE k_R = 0.39; // 1/s
    GPU_REALTYPE k_B = 6.3; // 1/s
    GPU_REALTYPE kp_B = 3; // 1/(uM*s)
    GPU_REALTYPE k_A = 50; // 1/(uM*s)
    GPU_REALTYPE k_Y = 100; // 1/(uM*s)
    GPU_REALTYPE k_Z = 30; // To be divided by [CheZ] 1/s * [CheZ]
    GPU_REALTYPE g_B = 1; // 1/(uM*s)
    GPU_REALTYPE g_Y = 0.1; // 1/(uM*s)

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
    GPU_REALTYPE B_t = 0.28; // uM
    GPU_REALTYPE R_t = 0.08; // uM
    GPU_REALTYPE Y_t = 9.7; // uM
    GPU_REALTYPE Z_t = 3.8; // uM

    GPU_REALTYPE pwDivider = (0.2 * 1.53);
};

class Kollmann2005Population : public SimplePopulation {
public:
    Kollmann2005Population(std::string name, shared_ptr<Environment2D> Env, Kollmann2005Parameters parameters, int nBacteria);
    Kollmann2005Population(std::string name, shared_ptr<Environment2D> env, Kollmann2005Parameters params) : SimplePopulation(name, env, params), params(params) {
        init();
    }

    Kollmann2005Population(shared_ptr<Environment2D> Env, H5::Group group);

    // Simulation
    void liveTimestep(double dt) override;
    virtual double getStabledt() override {return params.integrationMultiplyer*0.002; };
    void printInternals() override;

    // Storage
    void setupStorage(H5::Group storage) override;
    bool save() override;
    void closeStorage() override;

    REGISTER_DEC_TYPE(Kollmann2005Population);
protected:
    // Initialization
    void init();

    // Simulation

    void move(double dt) override;

    // Parameters
    Kollmann2005Parameters params;

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
        Kollmann2005Population *p;
        int methylationLevel;
    public:
        dTm(Kollmann2005Population *par, int methLevel): p(par), methylationLevel(methLevel) {}
        array rateofchange(array &input) override;
    };

    class dAp : public DifferentialEquation {
        Kollmann2005Population *p;
    public:
        dAp(Kollmann2005Population *par): p(par) {}
        array rateofchange(array &input) override;
    };

    class dYp : public DifferentialEquation {
        Kollmann2005Population *p;
    public:
        dYp(Kollmann2005Population *par): p(par) {}
        array rateofchange(array &input) override;
    };

    class dBp : public DifferentialEquation {
        Kollmann2005Population *p;
    public:
        dBp(Kollmann2005Population *par): p(par) {}
        array rateofchange(array &input) override;
    };

    // For simulation of stochasticity
    class dR : public DifferentialEquation {
        Kollmann2005Population *p;
    public:
        dR(Kollmann2005Population *par): p(par) {}
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
