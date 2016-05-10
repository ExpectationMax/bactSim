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

    std::vector<dim_t> internalDim(settings.dimensions.size() + 1);

    // Calculate dimensions of internal representations
    for (auto i = 0; i < settings.dimensions.size(); i++) {
        // closest position should be 1.5 * border size to outer border
        internalDim[i] = (dim_t) 3 * BORDER_SIZE + ceil(settings.dimensions[i]/settings.resolution);
    }

    // Last dimension is the number of ligands
    internalDim[internalDim.size() - 1] = settings.ligands.size();

    this->internal_dimensions = internalDim;

    ligandMapping = constant(0, ligands.size(), 2, u16);
    for(size_t i = 0; i < ligands.size(); i++) {
        ligandMapping(i, LIGANDID) = ligands[i].ligandId;
        ligandMapping(i, LIGANDINTERNAL) = i;
    }
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
        iter++;
    }
}
#ifndef NO_GRAPHICS

void Environment::setupVisualizationWindow(Window &win) {
    visualizationWin = &win;
    numLigands = this->ligands.size();
    rows = ceil(numLigands / 2);
    if(numLigands > 1)
        visualizationWin->grid(rows, 2);
}

void Environment::visualize(double normalizer) {
    if(numLigands > 1) {
        for(size_t i = 0; i < rows; i++) {
            for (size_t j = 0; j < 2 && 2*i + j < numLigands; j++) {
                array ligandDensity = this->getDensity(this->ligands[2*i + j].ligandId);
                visualizationWin->operator()(i, j).image(ligandDensity/normalizer, this->ligands[2*i + j].name.c_str());
            }
        }
    }
    else
        visualizationWin->image(this->densities/normalizer, this->ligands[0].name.c_str());

    visualizationWin->show();
}
#endif

void Environment::simulateTimeStep(void) {
    this->applyBoundaryCondition();
    this->calculateTimeStep();
}

array Environment::getLigandMapping(std::vector<int> ligandIds) {
    array mapping = constant(0, ligandIds.size(), u16);

    for(size_t i = 0; i < ligandIds.size(); i++) {
        mapping(i) = ligandMapping(ligandMapping(span, LIGANDID) == ligandIds[i], LIGANDINTERNAL);
    }
    return mapping;
}




