//
// Created by Max Horn on 22/04/16.
//

#ifndef CHEMOHYBRID_GPU_BACTERIALPOPULATION_H
#define CHEMOHYBRID_GPU_BACTERIALPOPULATION_H

#import <arrayfire.h>
#import "Environments/Environment.h"
#include <Environments/Environment2D.h>
#import "General/Ligands.hpp"
#import <map>

array randomDistribution(int dimension);

struct BacterialParameters {
    double swimmSpeed;
};

class Bacterial2DPopulation {
protected:
    Environment2D env;
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



    Bacterial2DPopulation(Environment2D env, BacterialParameters params, std::vector<LigandInteraction> interactions) : env(env), params(params) {
        std::vector<int> ligandIds;

        Kds = array((dim_t)interactions.size());
        uptakeRates  = array((dim_t)interactions.size());
        productionRates = array((dim_t)interactions.size());

        for(size_t i = 0; i < interactions.size(); i++) {
            ligandIds.push_back(interactions[i].ligandId);
            Kds(i) = interactions[i].Kd;
            uptakeRates(i) = interactions[i].uptakeRate;
            productionRates(i) = interactions[i].productionRate;
        }

        ligandmapping = env.getLigandMapping(ligandIds);
    };

public:
    Bacterial2DPopulation(Environment2D Env, BacterialParameters parameters, std::vector<LigandInteraction> interactions, int nBacteria) :
            Bacterial2DPopulation(Env, parameters, interactions) {
        dim4 size = Env.getSize();
        this->xpos = randomDistribution(nBacteria) * size[0];
        this->ypos = randomDistribution(nBacteria) * size[1];
        this->angle = 2 * af::Pi * randomDistribution(nBacteria);
        this->not_tumbling = randomDistribution(nBacteria).as(b8);
    }

    Bacterial2DPopulation(Environment2D Env, BacterialParameters parameters, std::vector<LigandInteraction> interactions, int nBacteria, double *initialx, double *initialy) :
    Bacterial2DPopulation(Env, parameters, interactions) {
        this->xpos = af::array(nBacteria, initialx);
        this->ypos = af::array(nBacteria, initialy);
        this->angle = 2 * af::Pi * af::randn(nBacteria);
    }

    void move(double dt) {
        this->xpos = af::cos(this->angle)*this->params.swimmSpeed*dt;
        this->ypos = af::sin(this->angle)*this->params.swimmSpeed*dt;
        
        array outofrange = this->xpos > env.getSize()[0];
        this->xpos(outofrange) = constant(env.getSize()[0], outofrange.dims(0));
        
        outofrange = this->xpos < 0;
        this->xpos(outofrange) = constant(env.getSize()[0], 0);

        outofrange = this->ypos > env.getSize()[1];
        this->ypos(outofrange) = constant(env.getSize()[0], outofrange.dims(0));

        outofrange = this->ypos < 0;
        this->ypos(outofrange) = constant(env.getSize()[1], 0);


    }

    void updateInterpolatedPositions();

    int getSize();
    void interactWithEnv(array individuals);
};


#endif //CHEMOHYBRID_GPU_BACTERIALPOPULATION_H
