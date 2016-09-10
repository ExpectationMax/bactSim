//
// Created by Max Horn on 10/05/16.
//

#ifndef CHEMOHYBRID_GPU_MODEL2D_H
#define CHEMOHYBRID_GPU_MODEL2D_H

#include <array>
#include "Environments/Environment2D.h"
#include "BacterialPopulations/BacterialPopulation.h"

struct bacteriumRef {
    shared_ptr<BacterialPopulation> population;
    int individual;
};

class Model2D {
#ifndef NO_GRAPHICS
    Window *populationsWin;
#endif
    int totalBacteria = 0;
    unique_ptr<H5::H5File> storage;
public:
    Model2D(shared_ptr<Environment2D> environment, std::vector<shared_ptr<BacterialPopulation>> populations);

    Model2D(H5::H5File &input);

    shared_ptr<Environment2D> env;

    std::vector<shared_ptr<BacterialPopulation>> bacterialPopulations;

    void simulateTimestep();

    GPU_REALTYPE simulateFor(GPU_REALTYPE t);

    void setupVisualizationWindows(Window &winDiff, Window &winPop);

    void setupStorage(std::string path, int savestep);

    void visualize();

    void setupStorage(H5::H5File &output, int savestep);

    void closeStorage();

    void save();
private:
    void init();
    int simulationsSinceLastSave = 0;
    int savestep = 1;
    array processBacteriaParallel();
    void processOverlappingBacteria(array &overlaping);
};



#endif //CHEMOHYBRID_GPU_MODEL2D_H
