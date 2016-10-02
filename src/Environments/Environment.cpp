//
// Created by Max Horn on 18/04/16.
//

#include "Environment.h"

Environment::Environment(EnvironmentSettings settings) {
    this->init(settings);
}

Environment::Environment(H5::Group group) {
    EnvironmentSettings envSettings;
    H5::StrType varstrtype(0, H5T_VARIABLE);

    // Environment parameters
//    group.openAttribute("dt").read(HDF5_GPUTYPE, &envSettings.dt);
    group.openAttribute("Resolution").read(H5::PredType::NATIVE_DOUBLE, &envSettings.resolution);
    group.openAttribute("Boundary condition").read(BoundaryCondition::getH5ReadType(), &envSettings.boundaryCondition);

    // Get original dimensions
    H5::Attribute dimsAttr = group.openAttribute("Dimensions");
    hsize_t ndim = 0;
    dimsAttr.getSpace().getSimpleExtentDims(&ndim);
    // Use resize to
    envSettings.dimensions.resize(ndim);
    dimsAttr.read(H5::PredType::NATIVE_DOUBLE, envSettings.dimensions.data());

    // Get ligand properties
    int nLigands = group.getNumObjs();
    envSettings.ligands.reserve(nLigands);
    for(int i = 0; i < nLigands; i++) {
        std::string name = group.getObjnameByIdx(i);
        H5::DataSet ligandData = group.openDataSet(name);
        Ligand lig;
        ligandData.openAttribute("Name").read(varstrtype, lig.name);
        ligandData.openAttribute("Properties").read(Ligand::getH5ReadType(), &lig);
        envSettings.ligands.push_back(lig);
    }
    Environment::init(envSettings);
    this->storage.reset(new H5::Group(group));
}


#ifndef NO_GRAPHICS
void Environment::setupVisualizationWindow(Window &win) {
    visualizationWin = &win;
    numLigands = this->ligands.size();
    if(numLigands > 1){
        rows = ceil(numLigands / 2);
        visualizationWin->grid(rows, 2);
    }
}

void Environment::visualize(double normalizer) {
    if(numLigands > 1) {
        for(size_t i = 0; i < rows; i++) {
            for (size_t j = 0; j < 2 && 2*i + j < numLigands; j++) {
                array ligandDensity = this->getDensity(this->ligands[2*i + j].ligandId);
                visualizationWin->operator()(i, j).image((ligandDensity/normalizer).as(af::dtype::f32), this->ligands[2*i + j].name.c_str());
            }
        }
    }
    else
        visualizationWin->image((densities/normalizer).as(af::dtype::f32), this->ligands[0].name.c_str());

    visualizationWin->show();
}
#endif

array Environment::getLigandMapping(std::vector<int> ligandIds) {
    array mapping = constant(0, ligandIds.size(), u16);

    for(size_t i = 0; i < ligandIds.size(); i++) {
        mapping(i) = ligandMapping(ligandMapping(span, LIGANDID) == ligandIds[i], LIGANDINTERNAL);
    }
    return mapping;

}

void Environment::init(EnvironmentSettings settings) {
    this->settings = settings;
    this->ligands = settings.ligands;

    this->resolution = settings.resolution;
    this->boundaryCondition = settings.boundaryCondition;

    std::vector<dim_t> internalDim(settings.dimensions.size() + 1);

    // Calculate dimensions of internal representations
    for (auto i = 0; i < settings.dimensions.size(); i++) {
        // Invert the order of axis, as providing dimensions in format zyx seems unituitive
        // closest position should be 1.5 * border size to outer border
        internalDim[settings.dimensions.size() - i - 1] = (dim_t) 4 * BORDER_SIZE + ceil(settings.dimensions[i]/settings.resolution);
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

void Environment::setupStorage(unique_ptr<H5::Group> storage) {
    // Store reference to group
    this->storage = std::move(storage);

    // Define datastypes and spaces used more often
    H5::DataSpace scalar(H5S_SCALAR);
    H5::StrType varstrtype(0, H5T_VARIABLE);

//    this->storage->createAttribute("dt", H5::PredType::IEEE_F64LE, scalar)
//            .write(HDF5_GPUTYPE, &this->dt);
    this->storage->createAttribute("Resolution", H5::PredType::IEEE_F64LE, scalar)
            .write(H5::PredType::NATIVE_DOUBLE, &this->resolution);
    this->storage->createAttribute("Boundary condition", BoundaryCondition::getH5SaveType(), scalar)
            .write(BoundaryCondition::getH5ReadType(),&this->boundaryCondition);

    hsize_t ndims = this->settings.dimensions.size();
    H5::DataSpace dimSpace(1, &ndims);
    this->storage->createAttribute("Dimensions", H5::PredType::IEEE_F64LE, dimSpace)
            .write(H5::PredType::NATIVE_DOUBLE, this->settings.dimensions.data());
}

void Environment::closeStorage() {
    this->storage.reset();
}


