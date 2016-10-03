//
// Created by Max Horn on 10/05/16.
//

#ifndef CHEMOHYBRID_GPU_MODEL2D_H
#define CHEMOHYBRID_GPU_MODEL2D_H

#include <array>
#include "Environments/Environment.h"
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
    Model2D(shared_ptr<Environment> environment, std::vector<shared_ptr<BacterialPopulation>> populations, double dt);

    Model2D(H5::H5File &input);

    shared_ptr<Environment> env;

    std::vector<shared_ptr<BacterialPopulation>> bacterialPopulations;

    void simulateTimestep();

    GPU_REALTYPE simulateFor(GPU_REALTYPE t, bool *continueSim);

#ifndef NO_GRAPHICS
    void setupVisualizationWindows(Window &winDiff, Window &winPop);

    void visualize();
#endif
    void setupStorage(std::string path, int savestep);



    void setupStorage(H5::H5File &output, int savestep);

    void closeStorage();

    void save();
private:
    double EnvironmentDt;
//    double PopulationDt;
    double Modeldt;
    int simulationsSinceLastSave = 0;
    int savestep = 1;

    void init();
    array processBacteriaParallel(double dt);
    void processOverlappingBacteria(array &overlaping, double dt);

    void processAllBacteriaParallel(double dt);
};



#endif //CHEMOHYBRID_GPU_MODEL2D_H
