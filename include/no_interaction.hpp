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

double trial_wavefunction(arma::mat& position, double alpha);

double dwf_dalpha(arma::mat& position, double beta,  bool interactions);

double local_energy_numerical(arma::mat& position, double alpha, double h);

double numerical_double_derivative(arma::mat& position, double alpha, double h, int dimension_index, int particle_index);

arma::vec quantum_force(arma::mat& position, int index, double alpha);

double greens_ratio(arma::vec& qf_old, arma::vec& qf_new, arma::mat& position_old, arma::mat& position_new, double time_step, double D, int index);
