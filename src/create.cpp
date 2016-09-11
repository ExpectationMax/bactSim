#include <csignal>
#include <random>
#include <arrayfire.h>
#include "Environments/Environment2D.h"
#include "BacterialPopulations/Kollmann2005Population.h"
#include "Models/Model2D.h"
#include "Solvers/ForwardEulerSolver.h"
#include "Solvers/RungeKuttaSolver.h"

#define NO_GRAPHICS

namespace
{
    bool continueSimulation = true;
}

void signal_handler(int signal)
{
    std::cout << "Got signal, stoping simulation..." << std::endl;
    continueSimulation = false;
}

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

    ESettings.resolution = 2;
    ESettings.dimensions = std::vector<double> {200, 200};

    BoundaryCondition boundaryCondition(BC_NEUMANN);
    boundaryCondition.xpos = -20;
    boundaryCondition.ypos = 0;
    ESettings.boundaryCondition = boundaryCondition;

    // Setup Ligands
    Ligand ligand1 = {"Ligand1", 0, 10, 0.0, 0.02, 500.0};
    ESettings.ligands.push_back(ligand1);
//    Ligand ligand2 = {"Ligand2", 1,  5, 0.0, 1.0, 500.0};
//    ESettings.ligands.push_back(ligand2);

    shared_ptr<Environment2D> simEnv(new Environment2D(ESettings));

    GPU_REALTYPE bactdt = 0.009;

    // Update randomness
    af::setSeed(time(NULL));

    std::vector<shared_ptr<BacterialPopulation>> populations;

    shared_ptr<Solver> BactSolver(static_cast<Solver *>(new RungeKuttaSolver));

    // Setup population 1
    std::vector<LigandInteraction> ligandInteractions1;
    LigandInteraction interaction11 = {0, 5, 0, 0, 0};
    ligandInteractions1.push_back(interaction11);

    Kollmann2005Parameters bactParams = {BactSolver, ligandInteractions1, bactdt, 20};
    populations.push_back(shared_ptr<BacterialPopulation>(static_cast<BacterialPopulation *>(new Kollmann2005Population("Population 1", simEnv, bactParams, 10))));


//    std::vector<LigandInteraction> ligandInteractions2;
//    LigandInteraction interaction21 = {1,0.8, 0, 0, 0};
//    ligandInteractions2.push_back(interaction21);

//    Kollmann2005Parameters bactParams2 = {BactSolver, ligandInteractions2, bactdt, 5};
//    populations.push_back(shared_ptr<BacterialPopulation>(static_cast<BacterialPopulation *>(new Kollmann2005Population("Population 2", simEnv, bactParams2, 20))));

    // Setup model
    Model2D mymodel(simEnv, populations, bactdt);
    mymodel.setupStorage("test.h5", 50);
    mymodel.save();

#ifndef NO_GRAPHICS
    Window diffusionwindow(1024, 512,"Simple Diffusion simulation");
    diffusionwindow.setColorMap(AF_COLORMAP_HEAT);
    Window populationwindow(1024,512, "Populations");
    mymodel.setupVisualizationWindows(diffusionwindow, populationwindow);
#endif
    std::signal(SIGINT, signal_handler);
    std::signal(SIGTERM, signal_handler);
    double simulationTime = 100;
    double simulatedTime = mymodel.simulateFor(simulationTime, &continueSimulation);
    std::cout << "Simulated for " << simulatedTime << " of " << simulationTime << std::endl;
    mymodel.closeStorage();
    return 0;
}
