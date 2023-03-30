# FYS4411-project1

In this project we use the variational Monte Carlo method to simulate Bose-Einstein condensate in spherical and elliptic potentials and estimate the ground state energy. The main code is written in C++, and uses the Armadillo library for linear algebra. The code is parallelized using OpenMP.

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