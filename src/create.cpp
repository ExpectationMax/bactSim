#include <arrayfire.h>
#include "Environments/Environment2D.h"
#include <random>
#include "BacterialPopulations/ExamplePopulation.h"
#include "Models/Model2D.h"
#include "Solvers/ForwardEulerSolver.h"
#include "Solvers/RungeKuttaSolver.h"
#include <memory>

int main(int argc, char** argv)
{
	char t_device_name[64] = {0};
	char t_device_platform[64] = {0};
	char t_device_toolkit[64] = {0};
	char t_device_compute[64] = {0};
	af::deviceInfo(t_device_name, t_device_platform, t_device_toolkit, t_device_compute);

    printf("Device name: %s\n", t_device_name);
    printf("Platform name: %s\n", t_device_platform);
    printf("Toolkit: %s\n", t_device_toolkit);
    printf("Compute version: %s\n", t_device_compute);

    // Setup Environment
    EnvironmentSettings ESettings;
    ESettings.resolution = 0.5;
    ESettings.dimensions = std::vector<double> {10, 10};
    ESettings.dt = 0.008;

    BoundaryCondition boundaryCondition(BC_PERIODIC);
    boundaryCondition.xpos = 0;
    boundaryCondition.ypos = 0;
    ESettings.boundaryCondition = boundaryCondition;

    // Setup Ligands

    Ligand ligand1 = {"Ligand1", 0, 10, 0.0, 0.0, 10.0};
    ESettings.ligands.push_back(ligand1);
    Ligand ligand2 = {"Ligand2", 1, 10, 0.0, 0.0, 10.0};
    ESettings.ligands.push_back(ligand2);
//    Ligand ligand3 = {"Ligand3", 2, 10.0, 0.0, 0.0, 15.0};
//    ESettings.ligands.push_back(ligand3);
//    Ligand ligand4 = {"Ligand4", 3, 10.0, 0.0, 0.05,20.0};
//    ESettings.ligands.push_back(ligand4);

    shared_ptr<Solver> solver(static_cast<Solver *>(new RungeKuttaSolver));
    shared_ptr<Environment2D> simEnv(new Environment2D(ESettings, solver));

    // Update randomness
    af::setSeed(time(NULL));

    std::vector<shared_ptr<BacterialPopulation>> populations;

    // Setup population 1
    std::vector<LigandInteraction> ligandInteractions1;

    LigandInteraction interaction11 = {0, 0, 0, 0, 0};
    ligandInteractions1.push_back(interaction11);
    LigandInteraction interaction12 = {1, 5, 0, 0, 0};
    ligandInteractions1.push_back(interaction12);

    ExampleParameters bactParams = {ligandInteractions1, ESettings.dt, 10};
    populations.push_back(shared_ptr<BacterialPopulation>(static_cast<BacterialPopulation *>(new ExamplePopulation("Ligand 1 eater", simEnv, bactParams, 20))));

    // Setup population 2
    std::vector<LigandInteraction> ligandInteractions2;
//    LigandInteraction interaction21 = {0, 0, 0, 0};
//    ligandInteractions2.push_back(interaction21);
    LigandInteraction interaction22 = {1, 0, 0, 0, 0};
    ligandInteractions2.push_back(interaction22);


    ExampleParameters bactParams2 = {ligandInteractions2, ESettings.dt, 10};
    populations.push_back(shared_ptr<BacterialPopulation>(new ExamplePopulation("Ligand2 eater", simEnv, bactParams2, 20)));

    // Setup model
    Model2D mymodel(simEnv, populations);
    mymodel.setupStorage("test.h5");
    H5::Group testgr;

    // Enable visualization
#ifndef NO_GRAPHICS
    Window diffusionwindow(1024, 512,"Simple Diffusion simulation");
    diffusionwindow.setColorMap(AF_COLORMAP_HEAT);
    Window populationwindow(1024,512, "Populations");
    mymodel.setupVisualizationWindows(diffusionwindow, populationwindow);
#endif

    time_t start= time(NULL);
    mymodel.save();
    for(int i =0; i < 1/ESettings.dt; i++) {
        mymodel.simulateTimestep();
        mymodel.save();
#ifndef NO_GRAPHICS
        mymodel.visualize();
#endif
    }
    mymodel.closeStorage();

    return 0;
}
