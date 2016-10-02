//
// Created by Max Horn on 24/09/16.
//

#ifndef BACTSIM_GPU_COORDINATE_INDEXER_H
#define BACTSIM_GPU_COORDINATE_INDEXER_H

#include <arrayfire.h>

using namespace af;

class CoordinateIndexer {
    dim4 dims;
public:
    CoordinateIndexer();
    CoordinateIndexer(array &A);
    array operator()(array &x, array &y);
    array operator()(array &x, array &y, array &z);
};


#endif //BACTSIM_GPU_COORDINATE_INDEXER_H
