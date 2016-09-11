//
// Created by Max Horn on 10/08/16.
//
#include "Solver.h"

SolverFactory::map_type *SolverFactory::map = NULL;

shared_ptr<Solver>
SolverFactory::createInstance(std::string const &s) {
    map_type::iterator it = getMap()->find(s);
    if(it == getMap()->end())
        return shared_ptr<Solver>();
    return it->second();
}

SolverFactory::map_type * SolverFactory::getMap() {
    // never delete'ed. (exist until program termination)
    // because we can't guarantee correct destruction order
    if(!map) { map = new map_type; }
    return map;
}
