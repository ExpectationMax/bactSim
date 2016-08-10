//
// Created by Max Horn on 08/08/16.
//

#include "Ligand.hpp"

H5::CompType LigandInteraction::getH5type() {
    H5::CompType interaction(sizeof(LigandInteraction));
    interaction.insertMember("ligandId", HOFFSET(LigandInteraction, ligandId), H5::PredType::NATIVE_UINT);
    interaction.insertMember("uptakeRate", HOFFSET(LigandInteraction, uptakeRate), H5::PredType::NATIVE_DOUBLE);
    interaction.insertMember("productionRate", HOFFSET(LigandInteraction, productionRate), H5::PredType::NATIVE_DOUBLE);
    interaction.insertMember("Kd", HOFFSET(LigandInteraction, Kd), H5::PredType::NATIVE_DOUBLE);
    return interaction;
}
