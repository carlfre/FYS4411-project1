#pragma once 
#include <iostream>
#include <armadillo>
#include <cmath>
#include <chrono>
#include <string>
#include <fstream>
#include <vector>

//function to calculate the local energy 
double local_energy(arma::mat position, double alpha);

//function to calculate the probability ratio used in the Metropolis algorithm (without the interaction)
double probability_ratio(arma::mat old_position, arma::mat new_position, double alpha);

void monte_carlo(int MC_cycles, double step, int N_particles, int N_dimensions);