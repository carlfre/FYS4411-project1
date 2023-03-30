# FYS4411 - Project 1

In this project we use the variational Monte Carlo method to estimate the ground state energy of a gas of weakly interacting bosons in spherical and elliptic potentials. The main code is written in C++, and uses the Armadillo library for linear algebra. The code is parallelized using OpenMP.

## Code structure
The repository is structured as follows:

``` bash
├── main.cpp
├── makefile
├── configs
│   ├── analytical.txt
│   ├── density.txt
│   ├── gradient.txt
│   ├── interactions_gradient.txt
│   ├── interactions.txt
│   ├── numerical.txt
│   ├── parallel.txt
│   ├── statistics.txt
│   └── timing.txt
├── include
│   ├── interaction.hpp
│   ├── no_interaction.hpp
│   ├── parallelization.hpp
│   └── vmc_walker.hpp
└── src
    ├── initialization.cpp
    ├── interaction.cpp
    ├── no_interaction.cpp
    ├── parallelization.cpp
    ├── sampling.cpp
    ├── statistical_analysis.py
    └── vmc_walker.cpp
```
The main.cpp file contains the main function, and calls the functions in the src folder. The configs folder contains the settings for each simulation problem. The include folder contains the header files for the classes used in the code. The src folder contains the implementation of the classes.

## Usage

The code is written in C++ and can be compiled using the makefile. The code can be compiled using

```bash
make
```

The code can be run using

```bash
./main.exe [name of simulation problem]
```
Simulation problems include 
* analytical (run vmc-calculation for alpha-values in interval [0.1, 1])
* numerical (compute energy using numerical estimate for kinetic energy)
* density (writes particle positions to file, used for plotting one-body density)
* statistics (writes energies to file, used for estimating ground state energy + computing error)
* gradient (gradient descent for finding optimal alpha-value)
* interactions (run with interactions)
* interactions_gradient (gradient descent run with interactions)
* parallel (evaluate speedup for parallelization)
* timing (for timing of numerical energy estimate)

The settings for each simulation can be changed in configs/[problem name].txt
