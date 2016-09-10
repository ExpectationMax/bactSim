//
// Created by Max Horn on 22/04/16.
//

#ifndef BACSIM_GPU_BACTERIALPOPULATION_H
#define BACSIM_GPU_BACTERIALPOPULATION_H

#include <arrayfire.h>
#include "Environments/Environment.h"
#include <Environments/Environment2D.h>
#include "General/Ligand.h"
#include <map>

#define REGISTER_DEC_TYPE(NAME) \
    static DerivedRegister<NAME> reg; \
    static const std::string type; \
    std::string getType() override


#define REGISTER_DEF_TYPE(NAME) \
    DerivedRegister<NAME> NAME::reg(#NAME); \
    const std::string NAME::type = #NAME; \
    std::string NAME::getType() {return NAME::type;}

struct BacterialParameters {
    BacterialParameters() {};
    BacterialParameters(GPU_REALTYPE dt) : dt(dt) {};
    GPU_REALTYPE dt;
};


class BacterialPopulation {
public:
    static shared_ptr<BacterialPopulation> createFromGroup(shared_ptr<Environment2D> env, H5::Group group);
    std::string name;
    virtual void interactWithEnv(int individual) = 0;
    virtual void interactWithEnv(array individuals) = 0;

    virtual std::string getType() = 0;
    virtual int getSize() = 0;
    virtual array getXpos() = 0;
    virtual array getYpos() = 0;
    virtual GPU_REALTYPE getdt();

    virtual void liveTimestep() = 0;

    virtual void setupStorage(H5::Group storage);
    virtual bool save() = 0;
    virtual void closeStorage() = 0;

    array concentrations;
    virtual array getInterpolatedPositions() = 0;
    virtual void printInternals() = 0;
protected:
    BacterialPopulation(H5::Group group);
    BacterialPopulation(BacterialParameters params): params(params) {};
    void setupBaseStorage(H5::Group storage);

    BacterialParameters params;
    shared_ptr<H5::Group> storage;
};

// This allows bacterial populations to be registered and later initialized based on a type string and a H5Group
template<typename T> shared_ptr<BacterialPopulation> createT(shared_ptr<Environment2D> env, H5::Group group) {
    return shared_ptr<BacterialPopulation>(
            static_cast<BacterialPopulation *>(new T(env, group)));
}

class BacteriaFactory {
    typedef std::map<std::string, shared_ptr<BacterialPopulation>(*)(shared_ptr<Environment2D> env, H5::Group group)> map_type;
public:
    static shared_ptr<BacterialPopulation>
    createInstance(std::string const &s, shared_ptr<Environment2D> env, H5::Group group);
protected:
    static map_type * getMap();
    static map_type * map;
};

template<typename T>
class DerivedRegister : BacteriaFactory {
public:
    DerivedRegister(std::string const& s) {
        getMap()->insert(std::make_pair(s, &createT<T>));
    }
};

#endif //CHEMOHYBRID_GPU_BACTERIALPOPULATION_H
