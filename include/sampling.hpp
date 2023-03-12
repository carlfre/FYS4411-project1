#pragma once 
#include <iostream>
#include <armadillo>
#include <cmath>
#include <chrono>
#include <string>
#include <fstream>
#include <vector>

void importance_sampling_with_interactions(
    arma::mat& position, 
    arma::mat& new_position, 
    arma::mat& relative_position, 
    arma::mat& relative_position_new,
    double alpha,
    double beta,
    double hard_core_radius,
    double time_step,
    double D,
    double random_number,
    int k,
    int& accepted_moves
    );

void importance_sampling_without_interactions(
    arma::mat& position, 
    arma::mat& new_position, 
    double alpha,
    double beta,
    double time_step,
    double D,
    double random_number,
    int k,
    int& accepted_moves);

void brute_force_sampling_with_interactions(    
    arma::mat& position, 
    arma::mat& new_position, 
    arma::mat& relative_position, 
    arma::mat& relative_position_new,
    arma::vec& stepping_vector,
    double alpha,
    double beta,
    double hard_core_radius,
    double step,
    double random_number,
    int k,
    int& accepted_moves);

void brute_force_sampling_without_interactions(
    arma::mat& position, 
    arma::mat& new_position, 
    arma::vec& stepping_vector,
    double alpha,
    double step,
    double random_number,
    int k,
    int& accepted_moves);

void sampling(arma::mat& position, arma::mat& new_position, double alpha, int k, double step, double time_step, double D, bool importance_sampling, int& accepted_moves, arma::vec& stepping_vector, double random_number, bool interactions, arma::mat& relative_position, arma::mat& relative_position_new, double gamma, double beta, double hard_core_radius);
