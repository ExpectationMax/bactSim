bactSim
=======
A library for the simulation of bacterial agents based on the ArrayFire GPU acceleration library

Requirements
------------

* CMake
* The ArrayFire acceleration library available at https://github.com/arrayfire/arrayfire
* The HDF5 library with C++ bindings


A few example videos can also be found [here](https://drive.google.com/drive/folders/0B5C_8h9Ht_tOUFh0TWdnNklPLWs)


Project structure
-----------------

```
src
├── BacterialPopulations
│   ├── BacterialPopulation.cpp       <-- parent class for bacterial populations
│   ├── BacterialPopulation.h
│   ├── Matthaeus2009Population.cpp   <-- bacterial population following the publication of
│   ├── Matthaeus2009Population.h         Matthäus et al. E. coli superdiffusion and chemotaxis-search strategy,
                                          precision, and motility. Biophysical Journal, 2009
│   ├── SimplePopulation.cpp          <-- Simple example population, constantly produces molecules and constantly
│   └── SimplePopulation.h                degrades them, keeps initial direction forever
│
├── CMakeLists.txt
├── Environments
│   ├── BoundaryCondition.cpp         <-- boundary conditions that implement PDE boundary
│   ├── BoundaryCondition.h               constraints and bacterial behavior at boundary
│   ├── ConstantEnvironment.h         <-- environment with constant never chainging concentraition
│   ├── Environment.cpp               <-- environment with PDEs governing ligand diffucion
│   ├── Environment.h
│   ├── EnvironmentBase.cpp           <-- parent class for environments
│   └── EnvironmentBase.h
│
├── Examples
│   ├── AnalysisandVisualisation.ipynb  <-- example of how to read data from simulation run
│   ├── CMakeLists.txt
│   ├── Example1.cpp                  <-- example for initializing an environment with two diffusing molecules
|   |                                     and a single molecule degrading and producing population
│   ├── Example2.cpp                  <-- example of how to implement a new bacterial population
│   ├── Example3.cpp                  <-- example of Matthäus2009 population in constant environment with
|   |                                     gaussians
│   ├── Example4.cpp                  <-- example of Matthäus2009 population producing and degrading same molecule
│   └── Example5.cpp                  <-- example of Matthäus2009Population and SimplePopulation in same
|                                         environment, producing and degrading molecules
│
├── General                           <-- convenience functions
│   ├── ArrayFireHelper.cpp
│   ├── ArrayFireHelper.h
│   ├── CoordinateIndexer.cpp
│   ├── CoordinateIndexer.h
│   ├── Ligand.cpp
│   ├── Ligand.h
│   ├── ScalableCoordinateIndexer.cpp
│   ├── ScalableCoordinateIndexer.h
│   ├── StorageHelper.cpp
│   ├── StorageHelper.h
│   └── Types.h
│
├── Models                            <-- implementation model that combines bacterial populations with environments
│   ├── Model2D.cpp
│   └── Model2D.h
│
├── Solvers                           <-- solvers implemented for ODEs and PDEs
│   ├── ForwardEulerSolver.cpp
│   ├── ForwardEulerSolver.h
│   ├── RungeKuttaSolver.cpp
│   ├── RungeKuttaSolver.h
│   ├── Solver.cpp
│   └── Solver.h
│
└── simulate.cpp                      <-- loads a model from a stored hdf5 file and runs the simulation
```
