//
// Created by Max Horn on 20/06/16.
//

#ifndef CHEMOHYBRID_GPU_RUNGEKUTTASOLVER_H
#define CHEMOHYBRID_GPU_RUNGEKUTTASOLVER_H

#include "Solver.h"
#include "Environments/EnvironmentBase.h"

class RungeKuttaSolver : public Solver{
public:
    RungeKuttaSolver() {};
    virtual void solveStep(DifferentialEquation &eq, array &initial_state, GPU_REALTYPE stepsize) const override;
    REGISTER_DEC_SOLVER(RungeKuttaSolver);
};


#endif //CHEMOHYBRID_GPU_RUNGEKUTTASOLVER_H
