#include <arrayfire.h>
#include "easylogging++.h"

#include "Environments/Environment2D.h"
#include <random>
#include "BacterialPopulations/BacterialPopulation.h"
#include "Models/Model2D.h"

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
    ESettings.dt = 0.0005;
    ESettings.dataType = f32;
    ESettings.convolutionType = CT_SERIAL;

    BoundaryCondition boundaryCondition(BC_PERIODIC);
    boundaryCondition.xpos = 0;
    boundaryCondition.ypos = 0;
    ESettings.boundaryCondition = boundaryCondition;

    // Setup Ligands
    Ligand ligand1 = {"Ligand1", 0, 10.0, 0.0, 0.0, 5.0};
    Ligand ligand2 = {"Ligand2", 1, 10.0, 0.0, 0.0, 10.0};
    //Ligand ligand3 = {"Ligand3", 2, 20.0, 0.0, 0.0, 15.0};
    //Ligand ligand4 = {"Ligand4", 3, 20.0, 0.0, 0.05,20.0};
    ESettings.ligands = std::vector<Ligand> {ligand1, ligand2};

    Environment2D simEnv(ESettings);

    // Setup population 1
    LigandInteraction interaction11 = {0, 5, 0, 0};
    LigandInteraction interaction12 = {1, 20, 0, 0};
    std::vector<LigandInteraction> ligandInteractions = {interaction11, interaction12};
    BacterialParameters bactParams = {ligandInteractions, ESettings.dt, 50};
    Bacterial2DPopulation bacteria(&simEnv, bactParams, 50);

    // Setup population 2
    LigandInteraction interaction21 = {0, 10, 0, 0};
    LigandInteraction interaction22 = {1, 40, 0, 0};
    std::vector<LigandInteraction> ligandInteractions2 = {interaction21, interaction22};
    BacterialParameters bactParams2 = {ligandInteractions2, ESettings.dt, 10};
    Bacterial2DPopulation bacteria2(&simEnv, bactParams, 50);

    // Setup model
    std::vector<Bacterial2DPopulation *> populations = { &bacteria, &bacteria2 };
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
