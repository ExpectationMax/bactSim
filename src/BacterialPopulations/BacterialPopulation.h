//
// Created by Max Horn on 22/04/16.
//

#ifndef CHEMOHYBRID_GPU_BACTERIALPOPULATION_H
#define CHEMOHYBRID_GPU_BACTERIALPOPULATION_H

#import <arrayfire.h>
#import "Environments/Environment.h"
#include <Environments/Environment2D.h>
#import "General/Ligand.hpp"
#import <map>

struct BacterialParameters {
    std::string name;
    std::vector<LigandInteraction> interactions;
    GPU_REALTYPE dt;
    GPU_REALTYPE swimmSpeed;
};

class Bacterial2DPopulation {
protected:
    Environment2D *env;
    array interpolatedPositions;
    int maxx;
    int maxy;
    int size;
    BacterialParameters params;
    array xpos;
    array ypos;
    array angle;
    array tumbling;

    array Kds;
    array uptakeRates;
    array productionRates;

    array ligandmapping;

    std::function<void(void)> validatePositions;

    Bacterial2DPopulation(Environment2D *env, BacterialParameters params);;

    void initializeAngleAndTumbiling() {
        angle = 2 * af::Pi * randu(size);
        tumbling = (randu(size)>= 0.5).as(b8);
    }

    void setPositions(array x, array y) {
        xpos = x;
        ypos = y;
        validatePositions();
        updateInterpolatedPositions();
    }

    void updateInterpolatedPositions();

    unique_ptr<H5::Group> storage;
    unique_ptr<H5::DataSet> xposStorage;
    unique_ptr<H5::DataSet> yposStorage;
    unique_ptr<H5::DataSet> angleStorage;
    unique_ptr<H5::DataSet> tumblingStorage;

public:
    static void applyPeriodicBoundary(int maxx, int maxy, array &xpos, array &ypos);
    static void applySolidBoundary(int maxx, int maxy, array &xpos, array &ypos, array &tumbling);

    Bacterial2DPopulation(Environment2D *Env, BacterialParameters parameters, int nBacteria) :
            Bacterial2DPopulation(Env, parameters) {
        size = nBacteria;

        initializeAngleAndTumbiling();

        array randx = randu(size) * maxx;
        array randy = randu(size) * maxy;
        setPositions(randx, randy);
    }

    Bacterial2DPopulation(Environment2D *Env, BacterialParameters parameters, int nBacteria, GPU_REALTYPE *initialx, GPU_REALTYPE *initialy) :
            Bacterial2DPopulation(Env, parameters) {
        size = nBacteria;

        initializeAngleAndTumbiling();

        setPositions(array(size, initialx), array(size, initialy));
    }

    void move();

    void interactWithEnv(int individual);

    void interactWithEnv(array individuals);

    void interactWithEnvPos(array pos);

    void simulate();

    int getSize() { return size; }

    array getXpos() { return xpos; }

    array getYpos() { return ypos; }

    void liveTimestep();

    void setupStorage(shared_ptr<H5::Group> storage);

    void closeStorage();

    void save();
};


#endif //CHEMOHYBRID_GPU_BACTERIALPOPULATION_H
