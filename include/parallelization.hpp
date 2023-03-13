#pragma once 
#include <iostream>
#include <armadillo>
#include <cmath>
#include <chrono>
#include <string>
#include <fstream>
#include <vector>

#include "../include/vmc_walker.hpp"

void merge_files(std::string filename, int n_walkers);


// TODO: change type to vector
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
