//
// Created by Max Horn on 20/06/16.
//

#ifndef CHEMOHYBRID_GPU_SOLVER_H
#define CHEMOHYBRID_GPU_SOLVER_H

#import <arrayfire.h>

using namespace af;

class DifferentialEquation {
public:
    virtual array rateofchange(array &input) = 0;
};

class Solver {
protected:
    Solver() {};
public:
    virtual array solveStep(DifferentialEquation &eq, array &inital_state, double stepsize) = 0;
};


#endif //CHEMOHYBRID_GPU_SOLVER_H
