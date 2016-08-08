//
// Created by Max Horn on 20/06/16.
//

#include "ForwardEulerSolver.h"

array ForwardEulerSolver::solveStep(DifferentialEquation &eq, array &inital_state, double stepsize) const {
    return inital_state + eq.rateofchange(inital_state)*stepsize;
}

