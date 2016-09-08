//
// Created by Max Horn on 20/06/16.
//

#ifndef CHEMOHYBRID_GPU_SOLVER_H
#define CHEMOHYBRID_GPU_SOLVER_H

#include <arrayfire.h>
#include <memory>
#include <map>
#include "General/Types.h"

#define REGISTER_DEC_SOLVER(NAME) \
    static RegisterSolver<NAME> reg; \
    static const std::string type; \
    std::string getType() override

#define REGISTER_DEF_SOLVER(NAME) \
    RegisterSolver<NAME> NAME::reg(#NAME); \
    const std::string NAME::type = #NAME; \
    std::string NAME::getType() {return NAME::type;}

using namespace af;
using std::shared_ptr;

class DifferentialEquation {
public:
    virtual array rateofchange(array &input) = 0;
};

class Solver {
protected:
    Solver() {};
public:
    virtual void solveStep(DifferentialEquation &eq, array &inital_state, GPU_REALTYPE stepsize) const = 0;
    virtual std::string getType() = 0;
};

// This allows solvers to be registered and later initialized based on a type string
template<typename T> shared_ptr<Solver> createSolver() {
    return shared_ptr<Solver>(
            static_cast<Solver *>(new T()));
}

class SolverFactory {
    typedef std::map<std::string, shared_ptr<Solver>(*)()> map_type;
public:
    static shared_ptr<Solver>
    createInstance(std::string const &s);
protected:
    static map_type * getMap();
    static map_type * map;
};

template<typename T>
class RegisterSolver : SolverFactory {
public:
    RegisterSolver(std::string const& s) {
        getMap()->insert(std::make_pair(s, &createSolver<T>));
    }
};


#endif //CHEMOHYBRID_GPU_SOLVER_H
