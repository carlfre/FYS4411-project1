#pragma once 
#include <iostream>
#include <armadillo>
#include <cmath>
#include <chrono>

//function to calculate the local energy 
double local_energy(double position, double alpha);

//function to calculate the probability ratio used in the Metropolis algorithm (without the interaction)
double probability_ratio(double old_position, double new_position, double alpha);

void monte_carlo(int MC_cycles, double step);