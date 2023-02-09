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
double probability_ratio(arma::mat& old_position, arma::mat& new_position, double alpha, int index);

arma::vec monte_carlo(double alpha, int MC_cycles, double step, int N_particles, int N_dimensions, bool importance_sampling, double time_step, bool numerical_double_derivative, double ndd_h);

arma::vec monte_carlo(double alpha, int MC_cycles, double step, int N_particles, int N_dimensions, bool importance_sampling, double time_step);

double local_energy_numerical(arma::mat& position, double alpha, double h);

double trial_wavefunction(arma::mat& position, double alpha);

void minimize_parameters(int MC_cycles, double step, int N_particles, int N_dimensions, bool importance_sampling, double time_step, double learning_rate, int max_iter);

double greens_ratio(arma::vec& qf_old, arma::vec& qf_new, arma::mat& position_old, arma::mat& position_new, double time_step, double D, int index);

arma::vec quantum_force(arma::mat& position, int index, double alpha);

double dwf_dalpha(arma::mat& position);