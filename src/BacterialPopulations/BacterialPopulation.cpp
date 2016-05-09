//
// Created by Max Horn on 22/04/16.
//

#include "BacterialPopulation.h"


void Bacterial2DPopulation::interactWithEnv(array individuals) {
    array pos = interpolatedPositions(span, individuals);
    // af_print(pos);
    array concentrations = env->getLigandConcentrations(pos, ligandmapping);
//    af_print(concentrations);
    array concentrationChange = -concentrations*uptakeRates*params.dt + productionRates*params.dt;
    //af_print(concentrationChange);
    // eval(concentrationChange);
    env->changeLigandConcentrationBy(concentrationChange, interpolatedPositions(pos, individuals), ligandmapping);
}

void Bacterial2DPopulation::validatePositions() {

    // Shift bacteria that are out of range back into field
    // x axis
    array outofrange = array(xpos.dims(0), 4, b8);
    outofrange(span, 0) = xpos > maxx;
    xpos -= outofrange(span, 0) * (xpos-maxx);

    outofrange(span, 1) = xpos < 0;
    xpos -= outofrange(span, 1) * xpos;

    // y axis

    outofrange(span, 2) = ypos > maxy;
    ypos -= outofrange(span, 2) * (ypos-maxy);

    outofrange(span, 3) = ypos < 0;
    ypos -= outofrange(span, 3) * ypos;

    // set tumbling to true for out of range bacteria
    tumbling = tumbling || max(outofrange, 1);
    // TODO: Check if speed improvement without this...
    // eval(xpos, ypos, tumbling);
}

void Bacterial2DPopulation::simulateTimestep() {
    validatePositions();
    updateInterpolatedPositions();
    array randomizer = af::randu(xpos.dims(0));
    array sorted_randomizer, shuffled_indices;
    sort(sorted_randomizer, shuffled_indices, randomizer, 0);
    for(auto i = 0; i < xpos.dims(0); i++) {
        array index = shuffled_indices(i);
        interactWithEnv(index);
    }
    //simulate();
    //move();
}

void Bacterial2DPopulation::simulate() {
    array concentrations = env->getLigandConcentrations(interpolatedPositions, ligandmapping);
}

void Bacterial2DPopulation::updateInterpolatedPositions() {
    interpolatedPositions = env->getInterpolatedPositions(xpos, ypos);
}


void Bacterial2DPopulation::move() {
    xpos += !tumbling*af::cos(angle)*params.swimmSpeed*params.dt;
    ypos += !tumbling*af::sin(angle)*params.swimmSpeed*params.dt;
}




