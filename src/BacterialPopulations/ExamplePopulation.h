//
// Created by Max Horn on 10/08/16.
//

#ifndef CHEMOHYBRID_GPU_EXAMPLEPOPULATION_H
#define CHEMOHYBRID_GPU_EXAMPLEPOPULATION_H

#include "BacterialPopulation.h"

struct ExampleParameters {
    std::vector<LigandInteraction> interactions;
    GPU_REALTYPE dt;
    GPU_REALTYPE swimmSpeed;
};

class ExamplePopulation : public BacterialPopulation {
protected:
    // Initialization
    ExamplePopulation(std::string name, shared_ptr<Environment2D> env, ExampleParameters params);
    void randomizeAngleAndTumbiling();
    void setPositions(array x, array y);

    // Boundary
    static void applyPeriodicBoundary(int maxx, int maxy, array &xpos, array &ypos);
    static void applySolidBoundary(int maxx, int maxy, array &xpos, array &ypos, array &tumbling);

    // Environment
    void updateInterpolatedPositions();
    std::function<void(void)> validatePositions;

    shared_ptr<Environment2D> env;
    array xpos;
    array ypos;
    array angle;
    array tumbling;
    array interpolatedPositions;
    int maxx;
    int maxy;
    int size;

    // Bacteria
    void interactWithEnvPos(array pos);
    void simulate();
    void move();

    ExampleParameters params;
    array Kon;
    array Koff;
    array uptakeRates;
    array productionRates;
    array ligandmapping;

    // Storage
    unique_ptr<H5::Group> storage;
    unique_ptr<H5::DataSet> xposStorage;
    unique_ptr<H5::DataSet> yposStorage;
    unique_ptr<H5::DataSet> angleStorage;
    unique_ptr<H5::DataSet> tumblingStorage;

public:
    ExamplePopulation(std::string name, shared_ptr<Environment2D> Env, ExampleParameters parameters, int nBacteria);
    ExamplePopulation(std::string name, shared_ptr<Environment2D> Env, ExampleParameters parameters, int nBacteria, GPU_REALTYPE *initialx, GPU_REALTYPE *initialy);
    ExamplePopulation(shared_ptr<Environment2D> Env, H5::Group group);


    void interactWithEnv(int individual) override ;
    void interactWithEnv(array individuals) override;

    int getSize() override { return size; }
    array getXpos() override { return xpos; }
    array getYpos() override { return ypos; }

    void liveTimestep() override;

    void setupStorage(H5::Group storage) override ;
    void save() override;
    void closeStorage() override;

    REGISTER_DEC_TYPE(ExamplePopulation);
};


#endif //CHEMOHYBRID_GPU_EXAMPLEPOPULATION_H
