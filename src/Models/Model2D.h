//
// Created by Max Horn on 10/05/16.
//

#ifndef CHEMOHYBRID_GPU_MODEL2D_H
#define CHEMOHYBRID_GPU_MODEL2D_H

#include <array>
#include "Environments/Environment2D.h"
#include "BacterialPopulations/BacterialPopulation.h"

struct bacteriumRef {
    Bacterial2DPopulation *population;
    int individual;
};

class Model2D {
    Environment2D *env = NULL;
#ifndef NO_GRAPHICS
    Window *populationsWin;
#endif
    std::vector<Bacterial2DPopulation *> bacterialPopulations;
    int totalBacteria = 0;
    bacteriumRef *allBacteria = NULL;
    unsigned int *callOrder = NULL;

public:
    Model2D(Environment2D *environment, std::vector<Bacterial2DPopulation *> populations);

    ~Model2D() { delete[] allBacteria; delete[] callOrder; }

    void simulateTimestep();

    void setupVisualizationWindows(Window &winDiff, Window &winPop);

    void visualize();
};


#endif //CHEMOHYBRID_GPU_MODEL2D_H
