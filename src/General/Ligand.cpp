//
// Created by Max Horn on 08/08/16.
//

#include "Ligand.hpp"

H5::CompType LigandInteraction::getH5type() {
    H5::CompType interaction(sizeof(LigandInteraction));
    interaction.insertMember("ligandId", HOFFSET(LigandInteraction, ligandId), H5::PredType::STD_U32LE);
    interaction.insertMember("uptakeRate", HOFFSET(LigandInteraction, uptakeRate), H5::PredType::IEEE_F64LE);
    interaction.insertMember("productionRate", HOFFSET(LigandInteraction, productionRate), H5::PredType::IEEE_F64LE);
    interaction.insertMember("Kon", HOFFSET(LigandInteraction, Kon), H5::PredType::IEEE_F64LE);
    interaction.insertMember("Koff", HOFFSET(LigandInteraction, Koff), H5::PredType::IEEE_F64LE);
    return interaction;
}
