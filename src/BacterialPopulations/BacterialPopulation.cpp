//
// Created by Max Horn on 10/08/16.
//

#include "BacterialPopulation.h"

BacteriaFactory::map_type *BacteriaFactory::map = NULL;

shared_ptr<BacterialPopulation>
BacteriaFactory::createInstance(std::string const &s, shared_ptr<Environment2D> env, H5::Group group) {
    map_type::iterator it = getMap()->find(s);
    if(it == getMap()->end())
        return shared_ptr<BacterialPopulation>();
    return it->second(env, group);
}

BacteriaFactory::map_type * BacteriaFactory::getMap() {
    // never delete'ed. (exist until program termination)
    // because we can't guarantee correct destruction order
    if(!map) { map = new map_type; }
    return map;
}

shared_ptr<BacterialPopulation> BacterialPopulation::createFromGroup(shared_ptr<Environment2D> env, H5::Group group) {
    H5::Attribute popType = group.openAttribute("Type");
    std::string type;
    H5::StrType varstrtype(0, H5T_VARIABLE);
    popType.read(varstrtype, type);
    return BacteriaFactory::createInstance(type, env, group);
}
