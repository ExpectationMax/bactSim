//
// Created by Max Horn on 02/10/16.
//

#ifndef BACTSIM_GPU_CONSTANTENVIRONMENT_H
#define BACTSIM_GPU_CONSTANTENVIRONMENT_H

#include "Environment.h"

class ConstantEnvironment : public Environment {
public:
    ConstantEnvironment(std::map<unsigned int, GPU_REALTYPE *> constantConcentrations, EnvironmentSettings settings): Environment(settings){
        dim4 dims(internal_dimensions.dims[0]-2*BORDER_SIZE, internal_dimensions.dims[1]-2*BORDER_SIZE);

        for(auto& ligandDensity: constantConcentrations) {
            densities(seq(BORDER_SIZE, end-BORDER_SIZE), seq(BORDER_SIZE, end-BORDER_SIZE), hostLigandMapping[ligandDensity.first]) =
                    array(dims, ligandDensity.second).as(AF_GPUTYPE);
        }
        eval(densities);
    }

    // Don't change Environment
    virtual void changeLigandConcentrationBy(array concDifferences, array positions, array weights, array ligands) override {}
    virtual void simulateTimestep(double dt) override {}
    virtual double getStabledt() override {return 1.0;}
};

#endif //BACTSIM_GPU_CONSTANTENVIRONMENT_H
