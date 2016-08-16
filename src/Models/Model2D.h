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
    std::vector<shared_ptr<BacterialPopulation>> bacterialPopulations;
    int totalBacteria = 0;
    bacteriumRef *allBacteria = NULL;
    unsigned int *callOrder = NULL;
    unique_ptr<H5::H5File> storage;
public:
    Model2D(shared_ptr<Environment2D> environment, std::vector<shared_ptr<BacterialPopulation>> populations);

    Model2D(H5::H5File &input);

    ~Model2D() { delete[] allBacteria; delete[] callOrder; }

    shared_ptr<Environment2D> env;

    void simulateTimestep();

    void setupVisualizationWindows(Window &winDiff, Window &winPop);

    void setupStorage(std::string path);

    void visualize();

    void setupStorage(H5::H5File &output);

    void closeStorage();

    void save();
private:
    void init();
};



#endif //CHEMOHYBRID_GPU_MODEL2D_H
