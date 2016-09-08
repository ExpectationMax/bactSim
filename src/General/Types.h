//
// Created by Max Horn on 08/09/16.
//

#ifndef BACTSIM_GPU_TYPES_H
#define BACTSIM_GPU_TYPES_H
#define GPU_REALTYPE float

#if GPU_REALTYPE == float
#define HDF5_GPUTYPE H5::PredType::NATIVE_FLOAT
#define AF_GPUTYPE f32
#else
#define HDF5_GPUTYPE H5::PredType::NATIVE_DOUBLE
#define AF_GPUTYPE f64
#endif
#endif //BACTSIM_GPU_TYPES_H
