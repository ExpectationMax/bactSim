//
// Created by Max Horn on 20/06/16.
//

#include "RungeKuttaSolver.h"

array RungeKuttaSolver::solveStep(DifferentialEquation &eq, array &initial_state, double stepsize) {
    array dxk = eq.rateofchange(initial_state);
    array xa = initial_state + 0.5*stepsize*dxk;
    array dxa = eq.rateofchange(xa);
    array xb = initial_state + 0.5*stepsize*dxa;
    array dxb = eq.rateofchange(xb);
    array xc = initial_state + stepsize*dxb;
    array dxc = eq.rateofchange(xc);
    array result = initial_state + stepsize/6 * (dxk + 2*(dxa + dxb) + dxc);
    result.eval();
    return result;
}