//
// Created by Max Horn on 09/09/16.
//

#ifndef BACTSIM_GPU_ARRAYFIREHELPER_H
#define BACTSIM_GPU_ARRAYFIREHELPER_H

#include <arrayfire.h>
using namespace af;

class ArrayFireHelper {
public:
    static array coordinateIndexing(array &A, array &x, array &y, array &z);
    static array coordinateIndexing(array &A, array &x, array &y);
    static array lastDimCoordinateIndexing(array &A, array &x, array &y, array &z);
    static array indexZAxis(array &A, array &indexes, array &z);
    static array isUnique(array A);
};


#endif //BACTSIM_GPU_ARRAYFIREHELPER_H
