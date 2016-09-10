//
// Created by Max Horn on 10/05/16.
//

#include "Model2D.h"
#include <random>
#include <H5Cpp.h>
#include <General/StorageHelper.h>
#include <General/ArrayFireHelper.h>

Model2D::Model2D(shared_ptr<Environment2D> environment, std::vector<shared_ptr<BacterialPopulation>> populations):
        env(environment), bacterialPopulations(populations) {
    init();
}

void Model2D::init() {
    for(auto population: this->bacterialPopulations) {
        totalBacteria += population->getSize();
    }
}

void Model2D::simulateTimestep() {
    // Let bacteria interact with environment
    array successful = processBacteriaParallel();
    array overlappingBacteria = where(!successful);
    processOverlappingBacteria(overlappingBacteria);

    // Simulate bacteria
    for(auto population: bacterialPopulations) {
        population->liveTimestep();
    }

    for(double ddt = 0; ddt < bacterialPopulations[0]->getdt(); ddt += env->dt){
        // Simulate environment
        env->evalDensities();
        env->simulateTimeStep();
    }

    simulationsSinceLastSave++;
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

void Model2D::setupStorage(H5::H5File &output, int saveStepsize) {
    this->storage.reset(new H5::H5File(output));
    // Let environment initialize its group
    unique_ptr<H5::Group> envGroup(new H5::Group(output.createGroup("Environment")));
    this->env->setupStorage(std::move(envGroup));

    // Let the bacterial populations create their groups below the populations group
    // As these functions should not really store the the Group but only use it in the setupStorage function we must not bother about ownership
    H5::Group popGroup(output.createGroup("Populations"));
    for(auto population: this->bacterialPopulations) {
        H5::Group curPopulation = popGroup.createGroup(population->name);
        population->setupStorage(curPopulation);
    }
    savestep = saveStepsize;
    this->storage->createAttribute("saveStep", PredType::INTEL_I32, StorageHelper::H5Scalar).write(PredType::NATIVE_INT, &savestep);

}

void Model2D::closeStorage() {
    for(auto population: this->bacterialPopulations) {
        population->closeStorage();
    }
    this->env->closeStorage();

    if(this->storage)
        this->storage.reset();
}

void Model2D::setupStorage(std::string path, int saveStepsize) {
    H5::H5File thefile (path, H5F_ACC_TRUNC);
    this->setupStorage(thefile, saveStepsize);
}

void Model2D::save() {
    if (!this->storage)
        return;
    if(simulationsSinceLastSave % savestep == 0) {
        this->env->save();
        for (auto population: this->bacterialPopulations) {
            population->save();
        }
        simulationsSinceLastSave = 0;
    }
}

Model2D::Model2D(H5::H5File &input) {
    shared_ptr<Environment2D> environment(new Environment2D(input.openGroup("Environment")));
    H5::Group populations = input.openGroup("Populations");
    int nPopulations = populations.getNumObjs();
    bacterialPopulations.reserve(nPopulations);
    for(int i = 0; i < populations.getNumObjs(); i++) {
        std::string name = populations.getObjnameByIdx(i);
        H5::Group popGroup = populations.openGroup(name);
        this->bacterialPopulations.push_back(BacterialPopulation::createFromGroup(environment, popGroup));
    }

    this->env = environment;
    this->bacterialPopulations = bacterialPopulations;
    init();
    this->storage = unique_ptr<H5::H5File>(new H5::H5File(input));
    this->storage->openAttribute("saveStep").read(H5::PredType::NATIVE_INT, &savestep);
}

GPU_REALTYPE Model2D::simulateFor(GPU_REALTYPE t) {
    // TODO: implement
}

array Model2D::processBacteriaParallel() {
    // Accumulate all interacting grid points from all populations
    array allPositions(totalBacteria, 4, af::dtype::u32);
    array indexes(totalBacteria, af::dtype::u32);
    std::vector<unsigned int> spaces;
    spaces.resize(bacterialPopulations.size());
    auto curSpace = 0;
    for(auto i =0; i<bacterialPopulations.size(); i++){
        spaces[i] = curSpace;
        auto curSize = bacterialPopulations[i]->getSize();
        allPositions(seq(curSpace, curSpace+curSize-1), span) = bacterialPopulations[i]->getInterpolatedPositions();
        indexes(seq(curSpace, curSpace+curSize-1)) = range(curSize);
        curSpace += curSize;
    }

    // Find those bacteria that uniquely interact with grid points --> can be calculated in parallel
    array filter = ArrayFireHelper::isUnique(allPositions(span, I_TOPLEFT))
                   && ArrayFireHelper::isUnique(allPositions(span, I_TOPRIGHT))
                   && ArrayFireHelper::isUnique(allPositions(span, I_BOTTOMLEFT))
                   && ArrayFireHelper::isUnique(allPositions(span, I_TOPRIGHT));


    // Perform parallel calculations
    for(auto i =0; i<bacterialPopulations.size(); i++){
        auto popRange = seq(spaces[i], spaces[i]+bacterialPopulations[i]->getSize()-1);
        array popIndexes = allPositions(popRange);
        array filtered = filter(popRange);
        array selector = popIndexes(filtered);
        bacterialPopulations[i]->interactWithEnv(selector);
    }
    return filter;
}

void Model2D::processOverlappingBacteria(array &overlappingBacteria) {
    int overlappingCount = overlappingBacteria.dims(0);
    // randomly process overlapping individuals
    if(overlappingCount) {
        std::cout << overlappingCount << " overlaping bacteria." <<std::endl;
        unsigned int *missed = overlappingBacteria.host<unsigned int>();
        std::vector<bacteriumRef> missingBacteria;
        missingBacteria.reserve(overlappingCount);
        for(auto i = 0; i<overlappingCount; i++) {
            int missingIndex = missed[i];
            int startIndex = 0;
            for(auto j = 0; j<bacterialPopulations.size(); j++) {
                if(startIndex <= missingIndex && missingIndex < startIndex+bacterialPopulations[j]->getSize()) {
                    bacteriumRef ref;
                    ref.population = bacterialPopulations[j];
                    ref.individual = missingIndex - startIndex;
                    missingBacteria.push_back(ref);
                    break;
                } else
                    startIndex += bacterialPopulations[j]->getSize();
            }
        }

        std::vector<int> seq(overlappingCount);
        for(auto i = 0; i< overlappingCount; i++)
            seq[i] = i;

        std::random_shuffle(seq.begin(), seq.end());
        for(size_t i = 0; i < overlappingCount; i++) {
            bacteriumRef curBacterium = missingBacteria[seq[i]];
            curBacterium.population->interactWithEnv(curBacterium.individual);
        }
    }
}

