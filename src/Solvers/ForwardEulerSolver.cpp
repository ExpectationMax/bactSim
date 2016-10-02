//
// Created by Max Horn on 20/06/16.
//

#include "ForwardEulerSolver.h"

void ForwardEulerSolver::solveStep(DifferentialEquation &eq, array &inital_state, GPU_REALTYPE stepsize) const {
    inital_state += eq.rateofchange(inital_state)*stepsize;
    eval(inital_state);
}

REGISTER_DEF_SOLVER(ForwardEulerSolver);
