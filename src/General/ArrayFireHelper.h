//
// Created by Max Horn on 09/09/16.
//

#ifndef BACTSIM_GPU_ARRAYFIREHELPER_H
#define BACTSIM_GPU_ARRAYFIREHELPER_H

#include <tuple>
#include <arrayfire.h>
using namespace af;

class ArrayFireHelper {
public:
    static array coordinateIndexing(array &A, array &x, array &y, array &z);
    static array coordinateIndexing(array &A, array &x, array &y);
    static array lastDimCoordinateIndexing(array &A, array &x, array &y, array &z);
    static array indexZAxis(array &A, array &indexes, array &z);
    static std::tuple<array, array, array> getOriginalIndexes(array &A, array &indexes);
    static array isUnique(array A);
    static array gammaSampler(unsigned int n, int shape, double scale, double location);

};


#endif //BACTSIM_GPU_ARRAYFIREHELPER_H
