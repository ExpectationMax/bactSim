//
// Created by Max Horn on 08/08/16.
//

#include "Ligand.h"

H5::CompType LigandInteraction::getH5SaveType() {
    H5::CompType interaction(sizeof(LigandInteraction));
    interaction.insertMember("ligandId", HOFFSET(LigandInteraction, ligandId), H5::PredType::STD_U32LE);
    interaction.insertMember("uptakeRate", HOFFSET(LigandInteraction, uptakeRate), H5::PredType::IEEE_F64LE);
    interaction.insertMember("Ku", HOFFSET(LigandInteraction, Ku), H5::PredType::IEEE_F64LE);
    interaction.insertMember("productionRate", HOFFSET(LigandInteraction, productionRate), H5::PredType::IEEE_F64LE);
    interaction.insertMember("Kon", HOFFSET(LigandInteraction, Kon), H5::PredType::IEEE_F64LE);
    interaction.insertMember("Koff", HOFFSET(LigandInteraction, Koff), H5::PredType::IEEE_F64LE);
    return interaction;
}

H5::CompType LigandInteraction::getH5ReadType() {
    H5::CompType interaction(sizeof(LigandInteraction));
    interaction.insertMember("ligandId", HOFFSET(LigandInteraction, ligandId), H5::PredType::NATIVE_UINT);
    interaction.insertMember("uptakeRate", HOFFSET(LigandInteraction, uptakeRate), H5::PredType::NATIVE_DOUBLE);
    interaction.insertMember("Ku", HOFFSET(LigandInteraction, Ku), H5::PredType::NATIVE_DOUBLE);
    interaction.insertMember("productionRate", HOFFSET(LigandInteraction, productionRate), H5::PredType::NATIVE_DOUBLE);
    interaction.insertMember("Kon", HOFFSET(LigandInteraction, Kon), H5::PredType::NATIVE_DOUBLE);
    interaction.insertMember("Koff", HOFFSET(LigandInteraction, Koff), H5::PredType::NATIVE_DOUBLE);
    return interaction;
}

// Sadly variable length types dont work that easily in CompType -->
// http://stackoverflow.com/questions/35477590/reading-a-hdf5-dataset-with-compound-data-type-containing-multiple-sets-with-var
// Have to set them manually

H5::CompType Ligand::getH5SaveType() {
    H5::CompType ligand(sizeof(Ligand));
    //H5::StrType varstrtype(0, H5T_VARIABLE);
    //ligand.insertMember("Name", HOFFSET(Ligand, name), varstrtype);
    ligand.insertMember("Ligand ID", HOFFSET(Ligand, ligandId), H5::PredType::STD_U32LE);
    ligand.insertMember("Initial concentration", HOFFSET(Ligand, initialConcentration), H5::PredType::IEEE_F64LE);
    ligand.insertMember("Global production rate", HOFFSET(Ligand, globalProductionRate), H5::PredType::IEEE_F64LE);
    ligand.insertMember("Global degradation rate", HOFFSET(Ligand, globalDegradationRate), H5::PredType::IEEE_F64LE);
    ligand.insertMember("Diffusion coefficient", HOFFSET(Ligand, diffusionCoefficient), H5::PredType::IEEE_F64LE);
    return ligand;
}

H5::CompType Ligand::getH5ReadType() {
    H5::CompType ligand(sizeof(Ligand));
    H5::StrType varstrtype(0, H5T_VARIABLE);
    //ligand.insertMember("Name", HOFFSET(Ligand, name), varstrtype);
    ligand.insertMember("Ligand ID", HOFFSET(Ligand, ligandId), H5::PredType::NATIVE_UINT);
    ligand.insertMember("Initial concentration", HOFFSET(Ligand, initialConcentration), H5::PredType::NATIVE_DOUBLE);
    ligand.insertMember("Global production rate", HOFFSET(Ligand, globalProductionRate), H5::PredType::NATIVE_DOUBLE);
    ligand.insertMember("Global degradation rate", HOFFSET(Ligand, globalDegradationRate), H5::PredType::NATIVE_DOUBLE);
    ligand.insertMember("Diffusion coefficient", HOFFSET(Ligand, diffusionCoefficient), H5::PredType::NATIVE_DOUBLE);
    return ligand;
}