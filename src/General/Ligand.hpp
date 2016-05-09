//
// Created by Max Horn on 05/05/16.
//

#ifndef CHEMOHYBRID_GPU_LIGANDS_H
#define CHEMOHYBRID_GPU_LIGANDS_H

struct Ligand {
    std::string name;
    int ligandId;
    double initialConcentration;
    double globalProductionRate;
    double globalDegradationRate;
    double diffusionCoefficient;
};

struct LigandInteraction {
    int ligandId;
    double uptakeRate;
    double productionRate;
    double Kd;
};

#endif //CHEMOHYBRID_GPU_LIGANDS_H
