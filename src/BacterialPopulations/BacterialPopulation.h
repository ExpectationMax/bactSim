//
// Created by Max Horn on 22/04/16.
//

#ifndef CHEMOHYBRID_GPU_BACTERIALPOPULATION_H
#define CHEMOHYBRID_GPU_BACTERIALPOPULATION_H

#import <arrayfire.h>
#import "Environments/Environment.h"
#include <Environments/Environment2D.h>
#import "General/Ligand.hpp"
#import <map>

#define REGISTER_DEC_TYPE(NAME) \
    static DerivedRegister<NAME> reg

#define REGISTER_DEF_TYPE(NAME) \
    DerivedRegister<NAME> NAME::reg(#NAME)


class BacterialPopulation {
public:
    virtual void interactWithEnv(int individual) = 0;
    virtual void interactWithEnv(array individuals) = 0;

    virtual int getSize() = 0;
    virtual array getXpos() = 0;
    virtual array getYpos() = 0;

    virtual void liveTimestep() = 0;

    virtual void setupStorage(shared_ptr<H5::Group> storage) = 0;
    virtual void save() = 0;
    virtual void closeStorage() = 0;
};


// This allows bacterial populations to be registered and later initialized based on a type string
template<typename T> BacterialPopulation * createT(H5::Group group) { return new T(group); }

struct BacteriaFactory {
    typedef std::map<std::string, BacterialPopulation*(*)()> map_type;

    static BacterialPopulation * createInstance(std::string const& s) {
        map_type::iterator it = getMap()->find(s);
        if(it == getMap()->end())
            return 0;
        return it->second();
    }

protected:
    static map_type * getMap() {
        // never delete'ed. (exist until program termination)
        // because we can't guarantee correct destruction order
        if(!map) { map = new map_type; }
        return map;
    }

private:
    static map_type * map;
};

template<typename T>
struct DerivedRegister : BacteriaFactory {
    DerivedRegister(std::string const& s) {
        getMap()->insert(std::make_pair(s, &createT<T>));
    }
};

#endif //CHEMOHYBRID_GPU_BACTERIALPOPULATION_H
