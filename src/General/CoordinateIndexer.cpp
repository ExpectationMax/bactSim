//
// Created by Max Horn on 24/09/16.
//

#include "CoordinateIndexer.h"

CoordinateIndexer::CoordinateIndexer(array &A) {
    dims = A.dims();
}

array CoordinateIndexer::operator()(array &x, array &y) {
    return dims[0]*y + x;
}

array CoordinateIndexer::operator()(array &x, array &y, array &z) {
    return dims[0]*dims[1]*z + CoordinateIndexer::operator()(x,y);
}

CoordinateIndexer::CoordinateIndexer() {
    // Create a dummy indexer
    dims = array(1,1,1,1).dims();
}



