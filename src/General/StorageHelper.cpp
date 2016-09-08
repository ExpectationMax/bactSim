//
// Created by Max Horn on 16/08/16.
//

#include "StorageHelper.h"

H5::DataSpace StorageHelper::H5Scalar = {H5S_SCALAR};

H5::StrType StorageHelper::H5VariableString = {0, H5T_VARIABLE};

