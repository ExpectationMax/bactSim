//
// Created by Max Horn on 16/08/16.
//

#ifndef BACTSIM_GPU_STORAGEHELPER_H
#define BACTSIM_GPU_STORAGEHELPER_H

#include <H5Cpp.h>
#include <arrayfire.h>
using namespace af;
using namespace H5;

class StorageHelper {
public:
    template <class T>  static array loadLastDataToGpu(DataSet data, DataType H5MemoryType, dtype arrayfireType)
    {
        // Determine dimensions of data
        DataSpace inputSpace = data.getSpace();
        hsize_t ndims = inputSpace.getSimpleExtentNdims();
        hsize_t outputdims = ndims - 1;
        hsize_t *dims = new hsize_t[ndims];
        hsize_t *start = new hsize_t[ndims];
        hsize_t *count = new hsize_t[ndims];
        inputSpace.getSimpleExtentDims(dims);

        // Only select the last timepoint
        start[0] = dims[0]-1;
        count[0] = 1;
        for(int i = 1; i < ndims; i++) {
            start[i] = 0;
            count[i] = dims[i];
        }
        inputSpace.selectHyperslab(H5S_SELECT_SET, count, start);

        // create temporary 1D space
        hsize_t dataLength = 1;
        for(int i = 1; i < ndims; i++){
            dataLength *= dims[i];
        }
        T *temp = new T[dataLength];
        DataSpace tempSpace(1, &dataLength);

        // Read data to memory
        data.read(temp, H5MemoryType, tempSpace, inputSpace);

        // create GpuArray of appropriate size
        array output;
        switch(outputdims) {
            case 1:
                output = array(dims[ndims-1], temp);
                break;
            case 2:
                // Invert order due to col first/row first differences
                output = array(dims[ndims-1], dims[ndims-2], temp).T();
                break;
            case 3:
                // Don't know if transpose trick works that well here...
                output = array(dims[ndims-1], dims[ndims-2], dims[ndims-3], temp).T();
                break;
            default:
                throw Exception("Unexpected number of dimensions (loadLastDataToGpu)");
        }

        delete[] temp;
        delete[] count;
        delete[] start;
        delete[] dims;

        // Account for col first/row first differences between arrayfire and HDF5
        if(outputdims < 3)
            output = output.as(arrayfireType);
        else
            output = output.T().as(arrayfireType);
        output.eval();
        return output;
    }

    template <class T> static void appendDataToDataSet(array data, DataSet target, DataType H5MemoryType){
        // read dims from HDF5
        DataSpace targetSpace = target.getSpace();
        hsize_t ndims = targetSpace.getSimpleExtentNdims();
        hsize_t *dims = new hsize_t[ndims];
        targetSpace.getSimpleExtentDims(dims);

        // extend dims to contain another timepoint
        hsize_t *newdims = new hsize_t[ndims];
        std::copy(dims, dims+ndims, newdims);
        newdims[0] += 1;
        target.extend(newdims);
        targetSpace = target.getSpace();

        // select region for storage
        hsize_t *start = new hsize_t[ndims];
        hsize_t *count = new hsize_t[ndims];
        start[0] = dims[0];
        for(int i = 1; i < ndims; i++)
            start[i] = 0;
        std::copy(dims, dims+ndims, count);
        count[0] = 1;
        targetSpace.selectHyperslab(H5S_SELECT_SET, count, start);

        // calculate source space
        hsize_t nelements = 1;
        for(int  i = 1; i < ndims; i++)
            nelements *= dims[i];
        H5::DataSpace sourceSpace(1, &nelements);

        // Copy data from gpu
        T *hostMem;
//        std::cout << ndims << std::endl;
        if((ndims - 1) == 2) {
            // 2D data must be transposed due to col-major storage
            hostMem = data.T().host<T>();
        } else
            // 1D data, must not be transposed
            hostMem = data.host<T>();

        // Write data to Dataset
        target.write(hostMem, H5MemoryType, sourceSpace, targetSpace);
        af::freeHost(hostMem);
    }

    static H5::DataSpace H5Scalar;
    static H5::StrType H5VariableString;
private:

};


#endif //CHEMOHYBRID_GPU_GPUHELPER_H
