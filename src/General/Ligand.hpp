//
// Created by Max Horn on 05/05/16.
//

#ifndef CHEMOHYBRID_GPU_LIGANDS_H
#define CHEMOHYBRID_GPU_LIGANDS_H

#ifndef __H5Cpp_H
#include <H5Cpp.h>
#endif

struct Ligand {
    std::string name;
    uint8_t ligandId;
    double initialConcentration;
    double globalProductionRate;
    double globalDegradationRate;
    double diffusionCoefficient;
};

class LigandInteraction {
public:
    int ligandId;
    double uptakeRate;
    double productionRate;
    double Kd;
    static H5::CompType getH5type();
};

#endif //CHEMOHYBRID_GPU_LIGANDS_H
