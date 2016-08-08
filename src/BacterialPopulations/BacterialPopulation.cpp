//
// Created by Max Horn on 22/04/16.
//

#include "BacterialPopulation.h"


void Bacterial2DPopulation::interactWithEnv(int individual) {
    array pos = interpolatedPositions(span, individual);
    interactWithEnvPos(pos);
}

void Bacterial2DPopulation::interactWithEnv(array individuals) {
    array pos = interpolatedPositions(span, individuals);
    interactWithEnvPos(pos);
}

void Bacterial2DPopulation::interactWithEnvPos(array pos) {
    array concentrations = env->getLigandConcentrations(pos, ligandmapping);
    array concentrationChange = -concentrations*uptakeRates*params.dt + productionRates*params.dt;
    env->changeLigandConcentrationBy(concentrationChange, pos, ligandmapping);
}


void Bacterial2DPopulation::applyPeriodicBoundary(int maxx, int maxy, array &xpos, array &ypos) {
    // x axis
    xpos -= (xpos > maxx) * xpos;

    xpos += (xpos < 0) * (-xpos + maxx);

    // y axis
    ypos -= (ypos > maxy) * ypos;

    ypos += (ypos < 0) * (-ypos + maxy);
};

void Bacterial2DPopulation::applySolidBoundary(int maxx, int maxy, array &xpos, array &ypos, array &tumbling) {
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
};

void Bacterial2DPopulation::liveTimestep() {
    simulate();
    move();
    validatePositions();
    updateInterpolatedPositions();
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




