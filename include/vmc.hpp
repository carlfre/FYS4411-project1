#pragma once 
#include <iostream>
#include <armadillo>
#include <cmath>
#include <chrono>
#include <string>
#include <fstream>
#include <vector>


arma::vec monte_carlo(double alpha, int MC_cycles, double step, int N_particles, int N_dimensions, bool importance_sampling, double time_step, bool numerical_double_derivative, double ndd_h, bool interactions, double gamma, double beta, double hard_core_radius);

arma::vec monte_carlo(double alpha, int MC_cycles, double step, int N_particles, int N_dimensions, bool importance_sampling, double time_step, bool numerical_double_derivative, double ndd_h);

arma::vec monte_carlo(double alpha, int MC_cycles, double step, int N_particles, int N_dimensions, bool importance_sampling, double time_step);

arma::vec monte_carlo(double alpha, int MC_cycles, double step, int N_particles, int N_dimensions, bool importance_sampling);

arma::vec monte_carlo(double alpha, int MC_cycles, double step, int N_particles, int N_dimensions, bool importance_sampling, double time_step, bool interactions, double gamma, double beta, double hard_core_radius);

arma::vec monte_carlo_parallelized(int N_walkers, double alpha, int MC_cycles, double step, int N_particles, int N_dimensions, bool importance_sampling, double time_step, bool numerical_double_derivative, double ndd_h, bool interactions, double gamma, double beta, double hard_core_radius);

void minimize_parameters(int MC_cycles, double step, int N_particles, int N_dimensions, bool importance_sampling, double time_step, double learning_rate, int max_iter, bool interactions, double gamma, double beta, double hard_core_radius);

void minimize_parameters(int MC_cycles, double step, int N_particles, int N_dimensions, bool importance_sampling, double time_step, double learning_rate, int max_iter);
