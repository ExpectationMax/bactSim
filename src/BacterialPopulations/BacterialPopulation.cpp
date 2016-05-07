//
// Created by Max Horn on 22/04/16.
//

#include "BacterialPopulation.h"

array randomDistribution(int dimension) {
    return (randn(dimension) + 1)/2;
}



void Bacterial2DPopulation::interactWithEnv(array individuals) {
    array concentrations = env.getLigandConcentrations(interpolatedPositions(span, individuals), ligandmapping);
//    for(auto interaction: this->ligandInteractions) {
//      //  if()
//    }
}

void Bacterial2DPopulation::updateInterpolatedPositions() {
    interpolatedPositions = env.getInterpolatedPositions(xpos, ypos);
}



