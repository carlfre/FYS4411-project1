#pragma once 
#include <iostream>
#include <armadillo>
#include <cmath>
#include <chrono>
#include <string>
#include <fstream>
#include <vector>

/* Computes the local energy for non-interacting particles in spherical potential*/ 
double local_energy(arma::mat& position, double alpha);

/* function to calculate the probability ratio used in the Metropolis algorithm (without the interaction) -
it is used to decide if a move is accepted */
double probability_ratio(arma::mat& old_position, arma::mat& new_position, double alpha, int particle_index);

/* Computes the trial wave function for non-interacting particles in spherical potential */
double trial_wavefunction(arma::mat& position, double alpha);

/* Computes the derivative of the trial wave function with respect to alpha */
double dwf_dalpha(arma::mat& position, double beta,  bool interactions);

/* Computes the local energy numerically */
double local_energy_numerical(arma::mat& position, double alpha, double h);

/* Computes the double derivative of the trial wave function numerically */
double numerical_double_derivative(arma::mat& position, double alpha, double h, int dimension_index, int particle_index);

/* Computes the quantum force for importance sampling, for non-interacting particles */
arma::vec quantum_force(arma::mat& position, int index, double alpha);

/* Computes ratio of greens functions, used for deciding whether to accept a move with importance sampling*/
double greens_ratio(arma::vec& qf_old, arma::vec& qf_new, arma::mat& position_old, arma::mat& position_new, double time_step, double D, int index);
