find_package(PkgConfig)

find_path(LIBHDF5_INCLUDE_DIR H5Cpp.h
        HINTS /usr/local/include)

find_library(LIBHDF5Cpp_LIBRARY NAMES hdf5_cpp)

find_library(LIBHDF5_LIBRARY NAMES hdf5)


include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set LIBHDF5_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(LibHDF5  DEFAULT_MSG
        LIBHDF5Cpp_LIBRARY LIBHDF5_LIBRARY LIBHDF5_INCLUDE_DIR)

mark_as_advanced(LIBHDF5_INCLUDE_DIR LIBHDF5Cpp_LIBRARY LIBHDF5_LIBRARY )

set(LIBHDF5_LIBRARIES ${LIBHDF5_LIBRARY} ${LIBHDF5Cpp_LIBRARY} )
set(LIBHDF5_INCLUDE_DIRS ${LIBHDF5_INCLUDE_DIR} )