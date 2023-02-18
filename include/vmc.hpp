#pragma once 
#include <iostream>
#include <armadillo>
#include <cmath>
#include <chrono>
#include <string>
#include <fstream>
#include <vector>

//function to calculate the local energy 
double local_energy(arma::mat& position, double alpha);

//function to calculate the probability ratio used in the Metropolis algorithm (without the interaction)
double probability_ratio(arma::mat& old_position, arma::mat& new_position, double alpha, int particle_index);

arma::vec monte_carlo(double alpha, int MC_cycles, double step, int N_particles, int N_dimensions, bool importance_sampling, double time_step, bool numerical_double_derivative, double ndd_h, bool interactions, double gamma, double beta, double hard_core_radius);

arma::vec monte_carlo(double alpha, int MC_cycles, double step, int N_particles, int N_dimensions, bool importance_sampling, double time_step, bool numerical_double_derivative, double ndd_h);

arma::vec monte_carlo(double alpha, int MC_cycles, double step, int N_particles, int N_dimensions, bool importance_sampling, double time_step);

arma::vec monte_carlo(double alpha, int MC_cycles, double step, int N_particles, int N_dimensions, bool importance_sampling);

arma::vec monte_carlo(double alpha, int MC_cycles, double step, int N_particles, int N_dimensions, bool importance_sampling, double time_step, bool interactions, double gamma, double beta, double hard_core_radius);

void sampling(arma::mat& position, arma::mat& new_position, double alpha, int k, double step, double time_step, double D, bool importance_sampling, int& accepted_moves, arma::vec& stepping_vector, double random_number, bool interactions, arma::mat& relative_position, arma::mat& relative_position_new, double gamma, double beta, double hard_core_radius);

double local_energy_numerical(arma::mat& position, double alpha, double h);

double trial_wavefunction(arma::mat& position, double alpha);

void minimize_parameters(int MC_cycles, double step, int N_particles, int N_dimensions, bool importance_sampling, double time_step, double learning_rate, int max_iter);

double greens_ratio(arma::vec& qf_old, arma::vec& qf_new, arma::mat& position_old, arma::mat& position_new, double time_step, double D, int index);

arma::vec quantum_force(arma::mat& position, int index, double alpha);

double dwf_dalpha(arma::mat& position);