#include "arrayfire.h"
#include "easylogging++.h"

#include "Environments/Environment2D.h"

#define GPU_DATATYPE float

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

    EnvironmentSettings ESettings;
    ESettings.resolution = 1;
    ESettings.dimensions = std::vector<double> {80, 80};
    Ligand ligand1 = {"Testligand1", 20.0, 0.0, 0.0, 0.0, 0.0, 10.0};
    Ligand ligand2 = {"Testligand2", 20.0, 0.0, 0.0, 0.0, 0.0, 10.0};
    ESettings.ligands = std::vector<Ligand> {ligand1, ligand1};
    ESettings.dt = 0.001;
    ESettings.dataType = f32;
    Window mywindow(1024, 1024,"Simple Diffusion simulation");
    mywindow.setColorMap(AF_COLORMAP_MOOD);
    ESettings.win = &mywindow;
    Environment2D simEnv(ESettings);
    //simEnv.printInternals();

    simEnv.test();

    return 0;
}
