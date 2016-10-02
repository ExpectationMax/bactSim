//
// Created by Max Horn on 08/09/16.
//

#ifndef BACTSIM_GPU_TYPES_H
#define BACTSIM_GPU_TYPES_H
#define FLOATTYPE 1
#define DOUBLETYPE 2

#define REALTYPE FLOATTYPE


#if REALTYPE == FLOATTYPE
#define GPU_REALTYPE float
#define HDF5_GPUTYPE H5::PredType::NATIVE_FLOAT
#define AF_GPUTYPE af::dtype::f32
#else
#define GPU_REALTYPE double
#define HDF5_GPUTYPE H5::PredType::NATIVE_DOUBLE
#define AF_GPUTYPE af::dtype::f64
#endif
#endif //BACTSIM_GPU_TYPES_H
