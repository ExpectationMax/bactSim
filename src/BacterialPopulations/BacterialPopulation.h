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
    BacterialParameters params;
    array xpos;
    array ypos;
    array angle;
    array not_tumbling;
    array sensed_concentrations;

    array Kds;
    array uptakeRates;
    array productionRates;

    array ligandmapping;

    Bacterial2DPopulation(Environment2D *env, BacterialParameters params) : env(env), params(params) {
        std::vector<int> ligandIds;
        af::setSeed(time(NULL));
        Kds = array((dim_t)params.interactions.size());
        uptakeRates  = array((dim_t)params.interactions.size());
        productionRates = array((dim_t)params.interactions.size());

        for(size_t i = 0; i < params.interactions.size(); i++) {
            ligandIds.push_back(params.interactions[i].ligandId);
            Kds(i) = params.interactions[i].Kd;
            uptakeRates(i) = params.interactions[i].uptakeRate;
            productionRates(i) = params.interactions[i].productionRate;
        }
        ligandmapping = env->getLigandMapping(ligandIds);
        maxx = env->getSize()[0];
        maxy = env->getSize()[1];

    };

public:
    Bacterial2DPopulation(Environment2D *Env, BacterialParameters parameters, int nBacteria) :
            Bacterial2DPopulation(Env, parameters) {
        xpos = 2 + randu(nBacteria) * (maxx-4);
        ypos = 2 + randu(nBacteria) * (maxy-4);
        validatePositions();

        angle = 2 * af::Pi * randu(nBacteria);
        not_tumbling = (randu(nBacteria)>= 0.5).as(b8);
    }

    Bacterial2DPopulation(Environment2D *Env, BacterialParameters parameters, int nBacteria, double *initialx, double *initialy) :
    Bacterial2DPopulation(Env, parameters) {
        xpos = af::array(nBacteria, initialx);
        ypos = af::array(nBacteria, initialy);
        validatePositions();

        angle = 2 * af::Pi * randu(nBacteria);
        not_tumbling = (randu(nBacteria)>= 0.5).as(b8);
    }

    void move();

    void validatePositions();
    void updateInterpolatedPositions();

    int getSize();
    void interactWithEnv(array individuals);

    void simulateTimestep();

    void simulate();
};


#endif //CHEMOHYBRID_GPU_BACTERIALPOPULATION_H
