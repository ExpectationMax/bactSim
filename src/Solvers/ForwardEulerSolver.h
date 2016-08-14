//
// Created by Max Horn on 20/06/16.
//

#ifndef CHEMOHYBRID_GPU_FORWARDEULERSOLVER_H
#define CHEMOHYBRID_GPU_FORWARDEULERSOLVER_H

#include "Solver.h"

class ForwardEulerSolver : public Solver{
public:
    virtual array solveStep(DifferentialEquation &eq, array &inital_state, double stepsize) const override;
    REGISTER_DEC_SOLVER(ForwardEulerSolver);
};


#endif //CHEMOHYBRID_GPU_FORWARDEULERSOLVER_H
