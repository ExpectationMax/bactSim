CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
project(BactSim-gpu)

# Compile with C++11 support.
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)

if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
    message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

CHECK_CXX_COMPILER_FLAG("-pg" COMPILER_SUPPORTS_PROFILER)
if(COMPILER_SUPPORTS_PROFILER)
    set(CMAKE_CXX_FLAGS_DEBUG "-pg ${CMAKE_CXX_FLAGS_DEBUG}")
else()
    message(WARNING "Profiling is not supported with this compiler.")
endif()

if(APPLE)
    SET(CMAKE_SKIP_BUILD_RPATH  FALSE) # orig value FALSE
    SET(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
    SET(CMAKE_INSTALL_RPATH "/usr/local/lib/")
endif()

# indicate we wish for the binaries to be written to the bin directory
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
# we have custom CMake configuration files, tell CMake to use them:
LIST(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMakeModules)

# build items in subdirectories
add_subdirectory(src)
