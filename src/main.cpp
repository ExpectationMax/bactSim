#include "arrayfire.h"
#include "easylogging++.h"

#include "Environment.h"

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

    EnvironmentSettings<GPU_DATATYPE> ESettings;
    ESettings.resolution = 1;
    ESettings.dimensions = std::vector<double> {5, 5};
    Ligand<GPU_DATATYPE> ligand = {"Testligand", 0.0, 0.0, 0.0, 5.0};
    ESettings.ligands = std::vector<Ligand<GPU_DATATYPE>> {ligand};
    ESettings.dt = 0.0001;
    ESettings.dataType = f32;

    Environment<GPU_DATATYPE> simEnv(ESettings);

    simEnv.printInternals();
    simEnv.test();

    return 0;
}
