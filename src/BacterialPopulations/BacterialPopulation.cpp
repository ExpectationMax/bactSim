//
// Created by Max Horn on 22/04/16.
//

#include "BacterialPopulation.h"


void Bacterial2DPopulation::interactWithEnv(array individuals) {
    array pos = interpolatedPositions(span, individuals);
    af_print(xpos(individuals));
    af_print(ypos(individuals));
    af_print(pos);
    array concentrations = env->getLigandConcentrations(pos, ligandmapping);
    af_print(concentrations);
    array concentrationChange = -concentrations*uptakeRates*params.dt + productionRates*params.dt;
    //af_print(concentrationChange);
    env->changeLigandConcentrationBy(concentrationChange, interpolatedPositions(pos, individuals), ligandmapping);
}

void Bacterial2DPopulation::validatePositions() {
    unsigned int *host_outofrange;

    array outofrange = array(xpos.dims(0), 4, b8);
//    af_print(outofrange);
//    af_print(this->xpos > maxx);
    outofrange(span, 0) = xpos > maxx;
    outofrange(span, 1) = xpos < 0;
    outofrange(span, 2) = ypos > maxy;
    outofrange(span, 3) = ypos < 0;

    host_outofrange = af::sum(outofrange, 0).host<unsigned int>();

    if(host_outofrange[0] + host_outofrange[1] +  host_outofrange[2] + host_outofrange[3]){
        if(host_outofrange[0]) {
            array tmp = outofrange(span, 0);
            xpos(tmp) = constant(maxx, host_outofrange[0], xpos.type());
            not_tumbling(tmp) = false;
        }
        if(host_outofrange[1]) {
            array tmp = outofrange(span, 1);
            xpos(tmp) = constant(0, host_outofrange[1], xpos.type());
            not_tumbling(tmp) = false;
        }
        if(host_outofrange[2]) {
            array tmp = outofrange(span, 2);
            ypos(tmp) = constant(maxy, host_outofrange[2], ypos.type());
            not_tumbling(tmp) = false;
        }
        if(host_outofrange[3]) {
            array tmp = outofrange(span, 3);
            ypos(tmp) = constant(0, host_outofrange[3], ypos.type());
            not_tumbling(tmp) = false;
        }
        eval(xpos, ypos, not_tumbling);
    }
    af::freeHost(host_outofrange);

}

void Bacterial2DPopulation::simulateTimestep() {

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
    // TODO: sometime results in errors if all acteria are swimming
    af_print(not_tumbling);
    xpos(not_tumbling) = af::cos(angle(not_tumbling))*params.swimmSpeed*params.dt;
    ypos(not_tumbling) = af::sin(angle(not_tumbling))*params.swimmSpeed*params.dt;

    validatePositions();
}




