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

    ESettings.resolution = 10;
    ESettings.dimensions = std::vector<double> {10000, 10000};

    BoundaryCondition boundaryCondition(BC_PERIODIC);
    boundaryCondition.xpos = 0;
    boundaryCondition.ypos = 0;
    ESettings.boundaryCondition = boundaryCondition;

    int maxi = int(ceil(ESettings.dimensions[0]/ESettings.resolution));
    int maxj = int(ceil(ESettings.dimensions[1]/ESettings.resolution));
    double size = 2000;
    std::cout << maxj << ' ' << maxj << std::endl;
//     Initialize with gaussian
    GPU_REALTYPE *initialValues = new GPU_REALTYPE[maxi*maxj];
    for(int i = 0; i < maxi; i++) {
        for(int j = 0; j < maxj; j++) {
//            initialValues[i*maxi+j] = 0;
            initialValues[i*maxi+j] = 1000 *
                    (
                     exp(-( (pow(i - maxi/2.0, 2)/(2*size/ESettings.resolution)) + (pow(j - maxj/2.0, 2)/(2*size/ESettings.resolution)) ))
                    + exp(-( (pow(i - maxi*0.25, 2)/(2*size/ESettings.resolution)) + (pow(j - maxj*0.25, 2)/(2*size/ESettings.resolution)) ))
                    + exp(-( (pow(i - maxi*0.75, 2)/(2*size/ESettings.resolution)) + (pow(j - maxj*0.75, 2)/(2*size/ESettings.resolution)) ))
                    + exp(-( (pow(i - maxi*0.25, 2)/(2*size/ESettings.resolution)) + (pow(j - maxj*0.75, 2)/(2*size/ESettings.resolution)) ))
                    + exp(-( (pow(i - maxi*0.75, 2)/(2*size/ESettings.resolution)) + (pow(j - maxj*0.25, 2)/(2*size/ESettings.resolution)) ))
                    );
        }
    }

    std::map<unsigned int, GPU_REALTYPE *> initLigands = {{0, initialValues}};


    // Setup Ligands
    Ligand ligand1 = {"Food", 0, 0.0, 0.0, 0.0, 500.0};
    ESettings.ligands.push_back(ligand1);
//    Ligand ligand2 = {"Attractor", 1,  5, 0.0, 0.5, 500.0};
//    ESettings.ligands.push_back(ligand2);


//    shared_ptr<Environment> simEnv(new ConstantEnvironment(initLigands, ESettings));
    shared_ptr<Environment> simEnv(new ConstantEnvironment(initLigands, ESettings));
    delete[] initialValues;
    GPU_REALTYPE bactdt = 0.01;

    // Update randomness
    af::setSeed(time(NULL));

    std::vector<shared_ptr<BacterialPopulation>> populations;

    shared_ptr<Solver> BactSolver(static_cast<Solver *>(new ForwardEulerSolver));

    // Setup population 1
    std::vector<LigandInteraction> ligandInteractions1;
    LigandInteraction interaction11 = {0, 0, 0, 5.0, 0, 0};
//    LigandInteraction interaction12 = {1,   0, 0, 0, 0};
    ligandInteractions1.push_back(interaction11);
//    ligandInteractions1.push_back(interaction12);

    Matthaeus2009Parameters bactParams = {BactSolver, ligandInteractions1, 30};
    populations.push_back(shared_ptr<BacterialPopulation>(static_cast<BacterialPopulation *>(new Matthaeus2009Population("Population 1", simEnv, bactParams, 500))));

//    std::vector<LigandInteraction> ligandInteractions2;
//    LigandInteraction interaction21 = {1,0.8, 0, 0, 0};
//    ligandInteractions2.push_back(interaction21);

//    Matthaeus2009Parameters bactParams2 = {BactSolver, ligandInteractions2, bactdt, 5};
//    populations.push_back(shared_ptr<BacterialPopulation>(static_cast<BacterialPopulation *>(new Matthaeus2009Population("Population 2", simEnv, bactParams2, 20))));

    // Setup model
    Model2D mymodel(simEnv, populations, bactdt);
    mymodel.setupStorage("gaussian_highres.h5", 50);
    mymodel.save();

#ifndef NO_GRAPHICS
    Window diffusionwindow(1024, 512,"Simple Diffusion simulation");
    diffusionwindow.setColorMap(AF_COLORMAP_HEAT);
    Window populationwindow(1024,512, "Populations");
    mymodel.setupVisualizationWindows(diffusionwindow, populationwindow);
#endif
    std::signal(SIGINT, signal_handler);
    std::signal(SIGTERM, signal_handler);
    double simulationTime = 500;
    double simulatedTime = mymodel.simulateFor(simulationTime, &continueSimulation);
    std::cout << "Simulated for " << simulatedTime << " of " << simulationTime << std::endl;
    mymodel.closeStorage();
    return 0;
}
