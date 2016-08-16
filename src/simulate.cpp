//
// Created by Max Horn on 15/08/16.
//

#include <H5Cpp.h>
#include "Models/Model2D.h"
int main(int argc, char** argv)
{
    Window diffusionwindow(1024, 512,"Simple Diffusion simulation");
    diffusionwindow.setColorMap(AF_COLORMAP_HEAT);
    Window populationwindow(1024,512, "Populations");
    H5::H5File input("test.h5", H5F_ACC_RDWR);
    Model2D mymodel(input);
    mymodel.setupVisualizationWindows(diffusionwindow, populationwindow);
//#ifndef NO_GRAPHICS

//
//
//
//#endif

    time_t start= time(NULL);

    for(int i =0; i < 2/mymodel.env->dt; i++) {
        std::cout << i << std::endl;
        mymodel.simulateTimestep();
        mymodel.save();
#ifndef NO_GRAPHICS
        mymodel.visualize();
#endif
    }
    mymodel.closeStorage();

}