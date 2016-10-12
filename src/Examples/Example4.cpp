/**
 * Example4: Simulation of Matthäus 2009 population and Simple population within the same environment
 */

#include <csignal>
#include <random>
#include <arrayfire.h>
#include "Environments/ConstantEnvironment.h"
#include "BacterialPopulations/Matthaeus2009Population.h"
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
    ESettings.resolution = 1;
    ESettings.dimensions = std::vector<double> {2000, 2000};

    BoundaryCondition boundaryCondition(BC_PERIODIC);
    boundaryCondition.xpos = 0;
    boundaryCondition.ypos = 0;
    ESettings.boundaryCondition = boundaryCondition;

    // Dummy ligand properties, not relevant as constant environment does not regard them
    Ligand ligand1 = {"Attractor", 0, 0.0, 0.0, 0.05, 500.0};
    ESettings.ligands.push_back(ligand1);
    shared_ptr<Environment> simEnv(new Environment(ESettings));
    GPU_REALTYPE bactdt = 0.01;

    af::setSeed(time(NULL));
    std::vector<shared_ptr<BacterialPopulation>> populations;

    shared_ptr<Solver> BactSolver(static_cast<Solver *>(new ForwardEulerSolver));

    // Setup population 1
    std::vector<LigandInteraction> ligandInteractions1;
    // Dummy interaction, currently populations without interactions are not supported
    LigandInteraction interaction11 = {0, 50, 0, 5.0, 0, 0};
    ligandInteractions1.push_back(interaction11);

    Matthaeus2009Parameters bactParams = {BactSolver, ligandInteractions1, 30};
    populations.push_back(shared_ptr<BacterialPopulation>(static_cast<BacterialPopulation *>(new Matthaeus2009Population("Population 1", simEnv, bactParams, 5000))));

    // Setup model
    Model2D mymodel(simEnv, populations, bactdt);
    mymodel.setupStorage("Example4.h5", 200);
    mymodel.save();

#ifndef NO_GRAPHICS
    Window diffusionwindow(1024, 512,"Simple Diffusion simulation");
    diffusionwindow.setColorMap(AF_COLORMAP_HEAT);
    Window populationwindow(1024,512, "Populations");
    mymodel.setupVisualizationWindows(diffusionwindow, populationwindow);
#endif
    std::signal(SIGINT, signal_handler);
    std::signal(SIGTERM, signal_handler);
    double simulationTime = 2400;
    double simulatedTime = mymodel.simulateFor(simulationTime, &continueSimulation);
    std::cout << "Simulated for " << simulatedTime << " of " << simulationTime << std::endl;
    mymodel.closeStorage();
    return 0;
}
