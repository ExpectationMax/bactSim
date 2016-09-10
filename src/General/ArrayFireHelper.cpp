//
// Created by Max Horn on 09/09/16.
//

#include "ArrayFireHelper.h"

array ArrayFireHelper::coordinateIndexing(array &A, array &x, array &y, array &z) {
    dim4 dim = A.dims();
    return dim[0]*dim[1]*z + dim[0]*y + x;
}

array ArrayFireHelper::lastDimCoordinateIndexing(array &A, array &x, array &y, array &z) {
    dim4 dim = A.dims();
    // reorder tile and flat required to ensure repitition of same value
    array scaledZIndexes = flat(tile(reorder(dim[0]*dim[1]*z, 1, 0), x.dims(0)));
    array scaledXYIndexes = tile(dim[0]*y + x, z.dims(0));
    return scaledZIndexes + scaledXYIndexes;
}

array ArrayFireHelper::coordinateIndexing(array &A, array &x, array &y) {
    dim4 dim = A.dims();
    return dim[0]*y + x;
}

array ArrayFireHelper::isUnique(array A) {
    array uniqueValues = setUnique(A);
    auto nUniques = uniqueValues.elements();
    auto nElements = A.elements();
    // TODO: Alternative? Requires a lot of memory and probably is not efficient
    array comparison = tile(A, 1, nUniques) == tile(reorder(uniqueValues, 1, 0), nElements);
    array counts = reorder(sum(comparison), 1, 0);
    array index = counts == 1;
    array unique = sum(comparison(span, index), 1);

    return unique;

}

array ArrayFireHelper::indexZAxis(array &A, array &indexes, array &z) {
    dim4 dim = A.dims();
    // reorder tile and flat required to ensure repitition of same value
    array scaledZIndexes = flat(tile(reorder(dim[0]*dim[1]*z, 1, 0), indexes.dims(0)));
    array scaledXYIndexes = tile(indexes, z.dims(0));
    return scaledZIndexes + scaledXYIndexes;
}

