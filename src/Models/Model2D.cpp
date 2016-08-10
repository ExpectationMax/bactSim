//
// Created by Max Horn on 10/05/16.
//

#include "Model2D.h"
#include <random>
#include <H5Cpp.h>

Model2D::Model2D(Environment2D *environment, std::vector<Bacterial2DPopulation *> populations):
        env(environment), bacterialPopulations(populations) {
    for(auto population: bacterialPopulations) {
        totalBacteria += population->getSize();
    }

    allBacteria = new bacteriumRef[totalBacteria];

    // This array should contain all indexes that can be called on allBacteria
    callOrder = new unsigned int[totalBacteria];
    int currentBacterium = 0;

    for(auto population: bacterialPopulations) {
        int currentSize = population->getSize();
        for(size_t i = 0; i < currentSize; i++) {
            allBacteria[currentBacterium].individual = i;
            allBacteria[currentBacterium].population = population;
            callOrder[currentBacterium] = currentBacterium;
            currentBacterium++;
        }
    }
}

void Model2D::simulateTimestep() {
    std::random_shuffle(&callOrder[0], &callOrder[totalBacteria-1]);
    for(size_t i = 0; i < totalBacteria; i++) {
        bacteriumRef curBacterium = allBacteria[i];
        curBacterium.population->interactWithEnv(curBacterium.individual);
    }

    for(auto population: bacterialPopulations) {
        population->liveTimestep();
    }
    env->simulateTimeStep();

}

void Model2D::setupVisualizationWindows(Window &winDiff, Window &winPop) {
    winDiff.setTitle("Ligand diffusion");
    env->setupVisualizationWindow(winDiff);

    populationsWin = &winPop;
    populationsWin->setTitle("Bacterial positions");
    if(bacterialPopulations.size() > 1)
        populationsWin->grid(1, bacterialPopulations.size());
}

void Model2D::visualize() {
    double normalizer = max<double>(env->getAllDensities());
    env->visualize(normalizer);
    if(bacterialPopulations.size() > 1) {
        for(size_t i = 0; i < bacterialPopulations.size(); i++)
            populationsWin->operator()(0, i).scatter(bacterialPopulations[i]->getXpos(), -1*bacterialPopulations[i]->getYpos());
    } else {
        populationsWin->scatter(bacterialPopulations[0]->getXpos(), -1*bacterialPopulations[0]->getYpos());
    }

    populationsWin->show();
}

void Model2D::setupStorage(H5::CommonFG &output) {

    // Let environment initialize its group
    unique_ptr<H5::Group> envGroup(new H5::Group(output.createGroup("Environment")));
    this->env->setupStorage(std::move(envGroup));

    // Let the bacterial populations create their groups below the populations group
    // As these functions should not really store the the Group but only use it in the setupStorage function we must not bother about ownership
    shared_ptr<H5::Group> popGroup(new H5::Group(output.createGroup("Populations")));
    for(auto population: this->bacterialPopulations) {
        population->setupStorage(popGroup);
    }
}

void Model2D::closeStorage() {
    for(auto population: this->bacterialPopulations) {
        population->closeStorage();
    }
    this->env->closeStorage();

    if(this->storage)
        this->storage->close();
}

void Model2D::setupStorage(std::string path) {
    this->storage.reset(new H5::H5File(path, H5F_ACC_TRUNC));
    this->setupStorage(*this->storage.get());
}

void Model2D::save() {
    if(!this->storage)
        return;

    this->env->save();
    for(auto population: this->bacterialPopulations) {
        population->save();
    }
}