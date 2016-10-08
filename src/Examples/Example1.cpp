//
// Created by Max Horn on 12/09/16.
//

#include <csignal>
#include "Environments/Environment.h"
#include "BacterialPopulations/Kollmann2005Population.h"
#include "Solvers/RungeKuttaSolver.h"
#include "Solvers/ForwardEulerSolver.h"
#include "Models/Model2D.h"
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

int main(int argc, char *argv[]) {
    // Initialize random number generator
    af::setSeed(time(NULL));

    // Setup Environment
    // =================

    EnvironmentSettings ESettings;
    // It should be $200 \mu m$ x $200 \mu m$ and resolve $0.8 \mu m$ as one grid point
    ESettings.resolution = 1;
    ESettings.dimensions = std::vector<double> {400, 400};
    ESettings.boundaryCondition = BoundaryCondition(BC_PERIODIC);

    // Create and Setup Ligands

    // Food, should be degraded slowly ($0.01 \frac{\mu M}{s}) and diffuse with $500 \frac{\u m^2}{s}
    Ligand ligand1;
    ligand1.name = "Food";
    ligand1.ligandId = 0;
    ligand1.initialConcentration = 10;
    ligand1.globalProductionRate = 0.1;
    ligand1.globalDegradationRate = 0;
    ligand1.diffusionCoefficient = 200.0;
    ESettings.ligands.push_back(ligand1);

    // Attractor () should be degraded at a higher rate ($0.5 \frac{\mu M}{s}) and diffuse with same speed
    Ligand ligand2 = {"Attractor", 1,  0, 0.0, 0.5, 200.0};
    ESettings.ligands.push_back(ligand2);

    // Setup the Environment by passing the settings
    shared_ptr<Environment> simEnv(new Environment(ESettings));

    // Setup populations
    // =================

    std::vector<shared_ptr<BacterialPopulation>> populations;

    // Setup Population 1 which should degrade "Food" and produce "Attractor"
    std::vector<LigandInteraction> ligandInteractions1;
    LigandInteraction foodInteraction;
    foodInteraction.ligandId = 0;
    foodInteraction.uptakeRate = 0.2;
    foodInteraction.productionRate = 0;
    // Ignore Kon and Koff as they are currently unused
    ligandInteractions1.push_back(foodInteraction);

    LigandInteraction attractorInteraction;
    attractorInteraction.ligandId = 1;
    attractorInteraction.uptakeRate = 0;
    attractorInteraction.productionRate = 0.5;
//    ligandInteractions1.push_back(attractorInteraction);

    // Setup Solver for differential equations
    shared_ptr<Solver> BactSolver(static_cast<Solver *>(new RungeKuttaSolver));

    // Initialize parameters: integrate using ForwardEulerSolver, swimmSpeed = $20 \frac{\mu m}{s}$
    Kollmann2005Parameters bactParams;
    bactParams.odesolver = BactSolver;
    bactParams.interactions = ligandInteractions1;
    bactParams.swimmSpeed = 10;

    populations.push_back(shared_ptr<BacterialPopulation>(
            new Kollmann2005Population("Population 1", simEnv, bactParams, 1000)
    ));


    // Setup Model
    // ===========

    // Setup model using environment and populations, simulate with dt=0.01
    Model2D mymodel(simEnv, populations, 0.05);

    // Save every 50th calculation step to the file Example1.h5
    mymodel.setupStorage("Example1.h5", 10);
    // Save the initial state!
    mymodel.save();

    // Simulation
    // ==========
    // Stop simulation in case of abort signal
    std::signal(SIGINT, signal_handler);
    std::signal(SIGTERM, signal_handler);
#ifndef NO_GRAPHICS
    // Setup visualisation windows
    Window diffusionwindow(1024, 512,"Simple Diffusion simulation");
    diffusionwindow.setColorMap(AF_COLORMAP_HEAT);
    Window populationwindow(1024,512, "Populations");
    mymodel.setupVisualizationWindows(diffusionwindow, populationwindow);
#endif
    // Run Simulation for 300s, second parameter is to halt the process in case of an abort signal
    mymodel.simulateFor(20, &continueSimulation);

    // Close storage
    mymodel.closeStorage();
}
