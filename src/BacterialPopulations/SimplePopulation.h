//
// Created by Max Horn on 10/08/16.
//

#ifndef BACTSIM_GPU_SIMPLEPOPULATION_H
#define BACTSIM_GPU_SIMPLEPOPULATION_H

#include "BacterialPopulation.h"

struct SimplePopulationParameters : BacterialParameters {
    SimplePopulationParameters() {};
    SimplePopulationParameters(std::vector<LigandInteraction> interactions, GPU_REALTYPE dt, GPU_REALTYPE swimmSpeed):
            BacterialParameters(dt), interactions(interactions), swimmSpeed(swimmSpeed) {};
    std::vector<LigandInteraction> interactions;
    GPU_REALTYPE swimmSpeed;
};

class SimplePopulation : public BacterialPopulation {
public:
    SimplePopulation(std::string name, shared_ptr<Environment2D> Env, SimplePopulationParameters parameters, int nBacteria);
    SimplePopulation(std::string name, shared_ptr<Environment2D> Env, SimplePopulationParameters parameters, int nBacteria, GPU_REALTYPE *initialx, GPU_REALTYPE *initialy);
    SimplePopulation(shared_ptr<Environment2D> Env, H5::Group group);

    virtual void printInternals() override;

    void interactWithEnv(int individual) override ;
    void interactWithEnv(array individuals) override;

    int getSize() override { return size; }
    array getXpos() override { return xpos; }
    array getYpos() override { return ypos; }

    void liveTimestep() override;
    void senseLigandConcentration();

    void setupStorage(H5::Group storage) override;
    bool save() override;
    void closeStorage() override;

    virtual array getInterpolatedPositions() override;

    REGISTER_DEC_TYPE(SimplePopulation);

protected:
    // Initialization
    SimplePopulation(std::string name, shared_ptr<Environment2D> env, SimplePopulationParameters params);
    void randomizeAngle();
    void setPositions(array x, array y);
    void init();

    // Boundary
    static void applyPeriodicBoundary(int maxx, int maxy, array &xpos, array &ypos);
    static void applySolidBoundary(int maxx, int maxy, array &xpos, array &ypos, array &atborder);

    // Environment
    void updateInterpolatedPositions();
    std::function<void(void)> validatePositions;

    shared_ptr<Environment2D> env;
    array xpos;
    array ypos;
    array angle;
    array atborder;
    array interpolatedPositions;
    array weights;
    array sensedConcentration;
    int maxx;
    int maxy;
    int size;

    // Bacteria
    virtual void interactWithEnvPos(array pos, array w, int individual);
    virtual void interactWithEnvPos(array pos, array weights, array individuals);
    virtual array modelUptakeProduction(array ligconc);
    virtual void simulate();
    virtual void move();

    SimplePopulationParameters params;
    array uptakeRates;
    array productionRates;
    array ligandmapping;

    // Storage
    H5::DataSpace storageSpace;
    H5::DSetCreatPropList storageProperties;
    unique_ptr<H5::DataSet> xposStorage;
    unique_ptr<H5::DataSet> yposStorage;
    unique_ptr<H5::DataSet> angleStorage;
};


#endif //CHEMOHYBRID_GPU_EXAMPLEPOPULATION_H
