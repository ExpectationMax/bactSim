//
// Created by Max Horn on 10/08/16.
//

#include "BoundaryCondition.h"

H5::CompType BoundaryCondition::getH5SaveType() {
    H5::CompType bc(sizeof(BoundaryCondition));
    bc.insertMember("Boundary condition type", HOFFSET(BoundaryCondition, type), getBCEnumType());
    // x
    bc.insertMember("xpos", HOFFSET(BoundaryCondition, xpos), H5::PredType::IEEE_F64LE);
    bc.insertMember("xneg", HOFFSET(BoundaryCondition, xneg), H5::PredType::IEEE_F64LE);
    // y
    bc.insertMember("ypos", HOFFSET(BoundaryCondition, ypos), H5::PredType::IEEE_F64LE);
    bc.insertMember("yneg", HOFFSET(BoundaryCondition, yneg), H5::PredType::IEEE_F64LE);
    // z
    bc.insertMember("zpos", HOFFSET(BoundaryCondition, zpos), H5::PredType::IEEE_F64LE);
    bc.insertMember("zneg", HOFFSET(BoundaryCondition, zneg), H5::PredType::IEEE_F64LE);
    return bc;
}

H5::CompType BoundaryCondition::getH5ReadType() {
    H5::CompType bc(sizeof(BoundaryCondition));
    bc.insertMember("Boundary condition type", HOFFSET(BoundaryCondition, type), getBCEnumType());
    // x
    bc.insertMember("xpos", HOFFSET(BoundaryCondition, xpos), H5::PredType::NATIVE_DOUBLE);
    bc.insertMember("xneg", HOFFSET(BoundaryCondition, xneg), H5::PredType::NATIVE_DOUBLE);
    // y
    bc.insertMember("ypos", HOFFSET(BoundaryCondition, ypos), H5::PredType::NATIVE_DOUBLE);
    bc.insertMember("yneg", HOFFSET(BoundaryCondition, yneg), H5::PredType::NATIVE_DOUBLE);
    // z
    bc.insertMember("zpos", HOFFSET(BoundaryCondition, zpos), H5::PredType::NATIVE_DOUBLE);
    bc.insertMember("zneg", HOFFSET(BoundaryCondition, zneg), H5::PredType::NATIVE_DOUBLE);
    return bc;
}

H5::EnumType BoundaryCondition::getBCEnumType() {
    H5::EnumType boundaryCondition(sizeof(BoundaryConditionType));
    BoundaryConditionType periodic, dirichelet, neumann;
    periodic = BC_PERIODIC; dirichelet = BC_DIRICHELET; neumann = BC_NEUMANN;
    boundaryCondition.insert("PERIODIC", &periodic);
    boundaryCondition.insert("DIRICHELET", &dirichelet);
    boundaryCondition.insert("NEUMANN", &neumann);
    return boundaryCondition;
}