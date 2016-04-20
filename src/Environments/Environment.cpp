//
// Created by Max Horn on 18/04/16.
//

#include "Environment.h"
#include <exception>
#include "../easylogging++.h"
#include <time.h>
#include "../progress.h"

Environment::Environment(EnvironmentSettings settings) {
    this->ligands = settings.ligands;

    this->dt = settings.dt;
    this->resolution = settings.resolution;
    this->boundaryCondition = settings.boundaryCondition;

    this->visualizationWin = settings.win;


    std::vector<dim_t> internalDim(settings.dimensions.size() + 1);

    // Calculate dimensions of internal representations
    //
    for (auto i = 0; i < internalDim.size(); i++) {
        internalDim[i] = (dim_t) 2 * BORDER_SIZE + floor(settings.dimensions[i]/settings.resolution);
    }

    // Last dimension is the number of ligands
    internalDim[internalDim.size() - 1] = settings.ligands.size();

    this->internal_dimensions = internalDim;

}

void Environment::printInternals() {
    af_print(this->densities);
    af_print(this->diffusion_filters);
}

std::vector<Ligand> Environment::getLigands() {
    return this->ligands;
}

void Environment::simulate(double advanceTime) {
    timer t = timer::start();
    unsigned iter = 0;
    double normalizer = max<double>(this->densities);
    while (progress(iter, t, advanceTime)) {
        this->simulateTimeStep();
        if(this->visualizationWin != NULL && (iter % 10) == 0){
            this->visualize(normalizer);
        }
        iter++;
    }
}

void Environment::visualize(double normalizer) {
    unsigned int numLigands = this->ligands.size();
    unsigned int rows = ceil(numLigands / 2);
    if(!this->WindowInitialized) {
        if(numLigands > 1)
            this->visualizationWin->grid(rows, 2);
        this->WindowInitialized = true;
    }
    //af_print(this->getDensity(0));
    //af_print(this->getDensity(1));
    /*
    this->visualizationWin->operator()(0, 0).image(this->getDensity(0), this->ligands[curLigand].name.c_str());
    this->visualizationWin->operator()(0, 1).image(this->getDensity(1), this->ligands[curLigand].name.c_str());
    this->visualizationWin->show();
    */
    //this->visualizationWin->grid(rows, 2);
    if(numLigands > 1) {
        for(size_t i = 0; i < rows; i++) {
            for (size_t j = 0; j < 2 && 2*i + j < numLigands; j++) {
                array ligandDensity = this->getDensity(2*i + j);
                this->visualizationWin->operator()(i, j).image(ligandDensity/normalizer, this->ligands[2*i + j].name.c_str());
            }
        }
    }
    else
        this->visualizationWin->image(this->densities/normalizer, this->ligands[0].name.c_str());



    this->visualizationWin->show();


}


