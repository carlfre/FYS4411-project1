#pragma once 
#include <iostream>
#include <armadillo>
#include <cmath>
#include <cassert>

/* Computes local energy for interacting particles in elliptic potential */
double local_energy_interaction(arma::mat& position, arma::mat& relative_position, double alpha, double beta, double a, double gamma);

/* Computes probability ratio for testing if move is accepted, for interacting particles*/
double probability_ratio_interaction(
    arma::mat& old_position, 
    arma::mat& new_position,
    arma::mat& relative_position, 
    arma::mat& relative_position_new, 
    double alpha, 
    double beta, 
    double hard_core_radius, 
    int particle_index);

/* Computes quantum force for importance sampling, for interacting particles */
arma::vec quantum_force_interaction(arma::mat& position, arma::mat& relative_position, double alpha, double beta, double hard_core_radius, int particle_index);
