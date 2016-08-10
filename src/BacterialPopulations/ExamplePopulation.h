//
// Created by Max Horn on 10/08/16.
//

#ifndef CHEMOHYBRID_GPU_EXAMPLEPOPULATION_H
#define CHEMOHYBRID_GPU_EXAMPLEPOPULATION_H

#include "BacterialPopulation.h"

struct BacterialParameters {
    std::string name;
    std::vector<LigandInteraction> interactions;
    GPU_REALTYPE dt;
    GPU_REALTYPE swimmSpeed;
};

class ExamplePopulation : public BacterialPopulation {
protected:
    // Initialization
    ExamplePopulation(shared_ptr<Environment2D> env, BacterialParameters params);
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

    BacterialParameters params;
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
    ExamplePopulation(shared_ptr<Environment2D> Env, BacterialParameters parameters, int nBacteria);
    ExamplePopulation(shared_ptr<Environment2D> Env, BacterialParameters parameters, int nBacteria, GPU_REALTYPE *initialx, GPU_REALTYPE *initialy);

    void interactWithEnv(int individual);
    void interactWithEnv(array individuals);

    int getSize() { return size; }
    array getXpos() { return xpos; }
    array getYpos() { return ypos; }

    void liveTimestep();

    void setupStorage(shared_ptr<H5::Group> storage);
    void save();
    void closeStorage();
};


#endif //CHEMOHYBRID_GPU_EXAMPLEPOPULATION_H
