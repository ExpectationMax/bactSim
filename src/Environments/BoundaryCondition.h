//
// Created by Max Horn on 10/08/16.
//

#ifndef CHEMOHYBRID_GPU_BOUNDARYCONDITION_H
#define CHEMOHYBRID_GPU_BOUNDARYCONDITION_H


// Define BoundaryConditions, this does not contain the implementation of the condition but merely its definition.
enum BoundaryConditionType {
    BC_PERIODIC,
    BC_DIRICHELET,
    BC_NEUMANN
};

class BoundaryCondition {
public:
    // Default allways use Neumann boundaries
    BoundaryCondition() {
        this->type = BC_NEUMANN;
    };
    BoundaryCondition(BoundaryConditionType boundaryConditionType) {
        this->type = boundaryConditionType;
    };
    BoundaryCondition(BoundaryConditionType boundaryConditionType, double xpos, double xneg) :
            BoundaryCondition(boundaryConditionType) {
        this->xpos = xpos;
        this->xneg = xneg;
    };
    BoundaryCondition(BoundaryConditionType boundaryConditionType, double xpos, double xneg, double ypos, double yneg) :
            BoundaryCondition(boundaryConditionType, xpos, xneg) {
        this->ypos = ypos;
        this->yneg = yneg;
    };
    BoundaryCondition(BoundaryConditionType boundaryConditionType, double xpos, double xneg, double ypos, double yneg,
                      double zpos, double zneg) : BoundaryCondition(boundaryConditionType, xpos, xneg, ypos, yneg) {
        this->zpos = zpos;
        this->zneg = zneg;
    };
    BoundaryConditionType type;
    double xneg, xpos, yneg, ypos, zneg, zpos = 0;
    static H5::EnumType getBCEnumType();
    static H5::CompType getH5SaveType();
    static H5::CompType getH5ReadType();
};


#endif //CHEMOHYBRID_GPU_BOUNDARYCONDITION_H
