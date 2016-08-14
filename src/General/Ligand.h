//
// Created by Max Horn on 05/05/16.
//

#ifndef CHEMOHYBRID_GPU_LIGANDS_H
#define CHEMOHYBRID_GPU_LIGANDS_H

#ifndef __H5Cpp_H
#include <H5Cpp.h>
#endif

class Ligand {
public:
    std::string name;
    unsigned int ligandId;
    double initialConcentration;
    double globalProductionRate;
    double globalDegradationRate;
    double diffusionCoefficient;
    static H5::CompType getH5SaveType();
    static H5::CompType getH5ReadType();
};

class LigandInteraction {
public:
    unsigned int ligandId;
    double uptakeRate;
    double productionRate;
    double Kon;
    double Koff;
    static H5::CompType getH5SaveType();
    static H5::CompType getH5ReadType();
};

#endif //CHEMOHYBRID_GPU_LIGANDS_H
