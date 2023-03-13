#pragma once 
#include <iostream>
#include <armadillo>
#include <cmath>
#include <cassert>

double local_energy_interactions(arma::mat& position, arma::mat& relative_position, double alpha, double beta, double a, double gamma);

double local_energy_naive(arma::mat& position, arma::mat& relative_position, double alpha, double beta, double a, double gamma);

double probability_ratio_naive(arma::mat& old_position, arma::mat& new_position,
arma::mat& relative_position, arma::mat& relative_position_new, double alpha, double beta, double hard_core_radius, int particle_index);

arma::vec quantum_force_naive(arma::mat& position, arma::mat& relative_position, double alpha, double beta, double hard_core_radius, int particle_index);
