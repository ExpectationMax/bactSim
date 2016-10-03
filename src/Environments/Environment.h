//
// Created by Max Horn on 19/04/16.
//

#ifndef CHEMOHYBRID_GPU_ENVIRONMENT2D_H
#define CHEMOHYBRID_GPU_ENVIRONMENT2D_H

#include "EnvironmentBase.h"
#include <vector>
#include "Solvers/Solver.h"
#include "General/CoordinateIndexer.h"
#include <H5Cpp.h>

#define I_TOPLEFT 0
#define I_TOPRIGHT 1
#define I_BOTTOMLEFT 2
#define I_BOTTOMRIGHT 3
#define W_TOPLEFT 0
#define W_TOPRIGHT 1
#define W_BOTTOMLEFT 2
#define W_BOTTOMRIGHT 3


class Environment : public EnvironmentBase {
    void init();
    // Boundary condition functions;
    std::function<void(void)> applyBoundaryCondition;
    static void applyNeumannBC(array &input, double resolution, BoundaryCondition &bc);
    static void applyDericheletBC(array &input,  BoundaryCondition &bc);
    static void applyPeriodicBC(array &input);
    static array getLaplacian();

//    array degradationRates;
//    array productionRates;

    array get_concentrations(array &indexes, array &ligands);

    std::map<unsigned int, unique_ptr<H5::DataSet>> ligands_storage;

    CoordinateIndexer densityIndexer;

public:
    Environment(EnvironmentSettings settings);
    array getAllDensities() override;
    array getDensity(int) override;
    std::vector<double> getSize() override;

    virtual double getStabledt() override;

    void setInterpolatedPositions(array &xpos, array &ypos, array &pos, array &weights);

    void changeLigandConcentrationBy(array concDifferences, array positions, array weights, array ligands);

    array getLigandConcentrations(array positions, array weights, array ligands);

    void evalDensities();

    virtual void simulateTimestep(double dt) override;

    virtual void closeStorage() override;

    virtual void setupStorage(unique_ptr<H5::Group> storage) override;

    virtual void save() override;

    Environment(H5::Group group);
};

#endif //CHEMOHYBRID_GPU_ENVIRONMENT2D_H
