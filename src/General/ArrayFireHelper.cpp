//
// Created by Max Horn on 09/09/16.
//

#include "ArrayFireHelper.h"
#include "Types.h"
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

std::tuple<array, array, array> ArrayFireHelper::getOriginalIndexes(array &A, array &indexes) {
    dim4 dim = A.dims();
    array x, y, z, zremain;
    z = af::floor(indexes/(dim[0]*dim[1]));
    zremain = mod(indexes, (dim[0]*dim[1]));
    y = af::floor(zremain/dim[0]);
    x = mod(zremain, dim[0]);
    return std::make_tuple(x, y, z);
}

array ArrayFireHelper::gammaSampler(unsigned int n, int shape, double scale, double location) {
    array x = constant(1.0, n, AF_GPUTYPE);
    for(auto i = 0; i < shape; i++)
        x *= randu(n, AF_GPUTYPE);
    x = scale*-log(x); + location;
    eval(x);
    return x;
}

