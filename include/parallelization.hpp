#pragma once 
#include <iostream>
#include <armadillo>
#include <cmath>
#include <chrono>
#include <string>
#include <fstream>
#include <vector>

#include "../include/vmc_walker.hpp"

/* Merge filenames of form [filename]_1, [filename]_2, [filename]_3, ..., [filename]_n 
to a file called [filename]*/
void merge_files(std::string filename, int n_walkers);


/* Creates multiple instances of the VMCWalker class, to run parallelized monte carlo computation. 
Can optionally send in filenames for density problem and energy estimation problem. */
arma::vec parallelized_mcmc(
    double alpha, 
    double beta,
    double gamma, 
    double step, // for no importance sampling
    double time_step, // for importance sampling
    double hard_core_radius, // for interactions
    double ndd_h, // h for finite difference double derivative
    int N_particles, 
    int N_dimensions, 
    bool importance_sampling, 
    bool numerical_double_derivative, 
    bool interactions,
    int MC_cycles,
    int n_walkers,
    std::string density_filename,
    std::string energy_filename
);


arma::vec parallelized_mcmc(
    double alpha, 
    double beta,
    double gamma, 
    double step, // for no importance sampling
    double time_step, // for importance sampling
    double hard_core_radius, // for interactions
    double ndd_h, // h for finite difference double derivative
    int N_particles, 
    int N_dimensions, 
    bool importance_sampling, 
    bool numerical_double_derivative, 
    bool interactions,
    int MC_cycles,
    int n_walkers
);
