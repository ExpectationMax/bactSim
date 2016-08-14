//
// Created by Max Horn on 18/04/16.
//

#include "Environment.h"

Environment::Environment(EnvironmentSettings settings, shared_ptr<Solver> odesolver) {
    this->init(settings, odesolver);
}

Environment::Environment(H5::Group group) {
    EnvironmentSettings envSettings;
    H5::StrType varstrtype(0, H5T_VARIABLE);

    group.openAttribute("dt").read(HDF5_GPUTYPE, &envSettings.dt);
    group.openAttribute("Resolution").read(H5::PredType::NATIVE_DOUBLE, &envSettings.resolution);
    group.openAttribute("Boundary condition").read(BoundaryCondition::getH5ReadType(), &envSettings.boundaryCondition);
    H5::Attribute dimsAttr = group.openAttribute("Dimensions");
    hsize_t ndim = 0;
    dimsAttr.getSpace().getSimpleExtentDims(&ndim);
    double *dims = new double[ndim];

    dimsAttr.read(H5::PredType::NATIVE_DOUBLE, dims);

    for(int i = 0; i < ndim; i++){
        if(dims[i] != 0)
            envSettings.dimensions.push_back(dims[i]);
    }
    // Get ligand properties
    int nLigands = group.getNumObjs();
    std::vector<Ligand> ligands(nLigands);
    for(int i = 0; i < nLigands; i++) {
        std::string name = group.getObjnameByIdx(i);
        H5::Group ligandGroup = group.openGroup(name);
        Ligand lig;
        ligandGroup.openAttribute("Name").read(varstrtype, lig.name);
        ligandGroup.openAttribute("Properties").read(Ligand::getH5ReadType(), &lig);
        ligands[i] = lig;
    }
    std::string solverName;
    group.openAttribute("Solver").read(varstrtype, solverName);
    shared_ptr<Solver> solver = SolverFactory::createInstance(solverName);
    this->init(envSettings, solver);
}


#ifndef NO_GRAPHICS
void Environment::setupVisualizationWindow(Window &win) {
    visualizationWin = &win;
    numLigands = this->ligands.size();
    if(numLigands > 1)
        rows = ceil(numLigands / 2);
        visualizationWin->grid(rows, 2);
}

void Environment::visualize(double normalizer) {
    if(numLigands > 1) {
        for(size_t i = 0; i < rows; i++) {
            for (size_t j = 0; j < 2 && 2*i + j < numLigands; j++) {
                array ligandDensity = this->getDensity(this->ligands[2*i + j].ligandId);
                visualizationWin->operator()(i, j).image(ligandDensity/normalizer, this->ligands[2*i + j].name.c_str());
            }
        }
    }
    else
        visualizationWin->image(this->densities/normalizer, this->ligands[0].name.c_str());

    visualizationWin->show();
}
#endif

void Environment::simulateTimeStep(void) {
    this->applyBoundaryCondition();
    this->calculateTimeStep();
}

array Environment::getLigandMapping(std::vector<int> ligandIds) {
    array mapping = constant(0, ligandIds.size(), u16);

    for(size_t i = 0; i < ligandIds.size(); i++) {
        mapping(i) = ligandMapping(ligandMapping(span, LIGANDID) == ligandIds[i], LIGANDINTERNAL);
    }
    return mapping;
}

void Environment::init(EnvironmentSettings settings, shared_ptr<Solver> odesolver) {
    this->odesolver = odesolver;
    this->ligands = settings.ligands;

    this->dt = settings.dt;
    this->resolution = settings.resolution;
    this->boundaryCondition = settings.boundaryCondition;

    std::vector<dim_t> internalDim(settings.dimensions.size() + 1);

    // Calculate dimensions of internal representations
    for (auto i = 0; i < settings.dimensions.size(); i++) {
        // Invert the order of axis, as providing dimensions in format zyx seems unituitive
        // closest position should be 1.5 * border size to outer border
        internalDim[settings.dimensions.size() - i - 1] = (dim_t) 2 * BORDER_SIZE + ceil(settings.dimensions[i]/settings.resolution);
    }

    // Last dimension is the number of ligands
    internalDim[internalDim.size() - 1] = settings.ligands.size();

    this->internal_dimensions = internalDim;

    // Store mapping of ligands on gpu
    this->ligandMapping = constant(0, ligands.size(), 2, u16);

    for(size_t i = 0; i < ligands.size(); i++) {
        this->hostLigandMapping[ligands[i].ligandId] = i;
        ligandMapping(i, LIGANDID) = ligands[i].ligandId;
        ligandMapping(i, LIGANDINTERNAL) = i;
    }
}


