//
// Created by Max Horn on 10/05/16.
//

#include "Model2D.h"

#include <random>

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