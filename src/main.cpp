#include <arrayfire.h>
#include "easylogging++.h"

#include "Environments/Environment2D.h"
#include <random>
#include "BacterialPopulations/BacterialPopulation.h"
#include "Models/Model2D.h"
#include "Solvers/ForwardEulerSolver.h"
#include "Solvers/RungeKuttaSolver.h"

INITIALIZE_EASYLOGGINGPP

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
    ESettings.dataType = f32;
    ESettings.convolutionType = CT_SERIAL;

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

    RungeKuttaSolver solver;
    Environment2D simEnv(ESettings, solver);

    // Update randomness
    af::setSeed(time(NULL));

    std::vector<Bacterial2DPopulation *> populations;

    // Setup population 1
    std::vector<LigandInteraction> ligandInteractions1;

    LigandInteraction interaction11 = {0, 5, 0, 0};
    ligandInteractions1.push_back(interaction11);
//    LigandInteraction interaction12 = {1, 5, 0, 0};
//    ligandInteractions1.push_back(interaction12);

    std::vector<LigandInteraction> ligandInteractions = {interaction11};//, interaction12};
    BacterialParameters bactParams = {ligandInteractions, ESettings.dt, 10};
    Bacterial2DPopulation bacteria(&simEnv, bactParams, 20);
    populations.push_back(&bacteria);

    // Setup population 2
    std::vector<LigandInteraction> ligandInteractions2;
//    LigandInteraction interaction21 = {0, 0, 0, 0};
//    ligandInteractions2.push_back(interaction21);
    LigandInteraction interaction22 = {1, 5, 0, 0};
    ligandInteractions2.push_back(interaction22);

    BacterialParameters bactParams2 = {ligandInteractions2, ESettings.dt, 10};
    Bacterial2DPopulation bacteria2(&simEnv, bactParams2, 20);
    populations.push_back(&bacteria2);

    // Setup model


    Model2D mymodel(&simEnv, populations);

    // Enable visualization
#ifndef NO_GRAPHICS
    Window diffusionwindow(1024, 512,"Simple Diffusion simulation");
    diffusionwindow.setColorMap(AF_COLORMAP_HEAT);
    Window populationwindow(1024,512, "Populations");
    mymodel.setupVisualizationWindows(diffusionwindow, populationwindow);
#endif


    time_t start= time(NULL);
//    simEnv.test();

    for(int i =0; i < 20/ESettings.dt; i++) {
        mymodel.simulateTimestep();
#ifndef NO_GRAPHICS
        mymodel.visualize();
#endif
    }

//    Window bacteriaWin(512, 512,"Simulation of bacterial population");
//    std::cout << time(NULL)-start;
//    for(int i =0; i < 20/ESettings.dt; i++) {
//        double normalizer = max<double>(simEnv.getAllDensities());
//        bacteria.simulateTimestep();
//        bacteriaWin.scatter(bacteria.getXpos(), bacteria.getYpos());
//        simEnv.evalDensities();
//
//        simEnv.simulateTimeStep();
//        simEnv.visualize(normalizer);
//
//        if(i % 200 == 0){
//            std::cout << pow((double)(time(NULL)-start)/200.0, -1) << " iterations per second" << std::endl;
//            start = time(NULL);
//        }
//
//    }

    // simEnv.test();
    // simEnv.printInternals();

    // array test = (randomDistribution(1000000) >= 0.5).as(b8);

    // printf("%f\n", sum<float>(test)/1000000);

    // simEnv.test();

    return 0;
}
