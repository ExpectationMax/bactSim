//
// Created by Max Horn on 20/06/16.
//

#ifndef CHEMOHYBRID_GPU_RUNGEKUTTASOLVER_H
#define CHEMOHYBRID_GPU_RUNGEKUTTASOLVER_H

#include "Solver.h"

class RungeKuttaSolver : Solver{
    array solveStep(DifferentialEquation, double stepsize) override;
};


#endif //CHEMOHYBRID_GPU_RUNGEKUTTASOLVER_H
