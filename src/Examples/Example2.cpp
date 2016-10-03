//
// Created by Max Horn on 17/09/16.
//
#include <H5Cpp.h>
#include "BacterialPopulations/SimplePopulation.h"
#include "General/StorageHelper.h"
using namespace af;

class ExamplePopulation: public SimplePopulation {
public:
    ExamplePopulation(std::string name, shared_ptr<Environment> Env, SimplePopulationParameters parameters, int nBacteria, int newParameter): SimplePopulation(name, Env, parameters, nBacteria) {
        theNewParameter = constant(newParameter, size);
    }

    ExamplePopulation(shared_ptr<Environment> Env, H5::Group group) : SimplePopulation(Env, group) {
        auto newParameter = group.openDataSet("newParameter");
        theNewParameter = StorageHelper::loadLastDataToGpu<int>(newParameter, H5::PredType::NATIVE_INT, af::dtype::s32);
    }

    void setupStorage(H5::Group storage) override {
        SimplePopulation::setupStorage(storage);
        newStorage.reset(new H5::DataSet(this->storage->createDataSet("newParameter", H5::PredType::INTEL_I32, this->storageSpace, this->storageProperties)));
    }

    void closeStorage() override {
        newStorage.reset();
        SimplePopulation::closeStorage();
    }

    bool save() override {
        if(SimplePopulation::save())
            StorageHelper::appendDataToDataSet<int>(theNewParameter, *newStorage, H5::PredType::NATIVE_INT);
    }

    REGISTER_DEC_TYPE(ExamplePopulation);
private:
    array theNewParameter;
    std::unique_ptr<H5::DataSet> newStorage;
};

REGISTER_DEF_TYPE(ExamplePopulation);