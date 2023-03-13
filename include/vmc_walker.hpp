#pragma once 
#include <iostream>
#include <armadillo>
#include <cmath>
#include <chrono>
#include <string>
#include <fstream>
#include <vector>

class VMCWalker{
    public:
        arma::mat position;
        arma::mat new_position;
        arma::mat relative_position;
        arma::mat relative_position_new;
        arma::vec qf;
        arma::vec qf_new;
        // arma::vec gaussian_random;
        arma::vec stepping_vector;
        int N_particles;
        int N_dimensions;
        double alpha;
        double beta;
        double gamma;
        double hard_core_radius;
        double time_step;
        double step;
        double D;
        double random_number;
        double ndd_h;
        int accepted_moves;
        bool interactions;
        bool importance_sampling;
        bool numerical_double_derivative;
        // random number generation
        // mersenne twister
        std::mt19937_64 generator;
        // uniform distribution
        std::uniform_real_distribution<double> rng_double;
        // normal distribution
        std::normal_distribution<double> rng_normal;


        // vmc_walker.cpp
        VMCWalker(
            double alpha, 
            double beta,
            double gamma, 
            double step, // for no importance sampling
            double time_step, // for importance sampling
            double hard_core_radius, // for interactions
            double ndd_h, // h for finite difference double derivative
            int N_particles, 
            int N_dimensions, 
            bool importance_sampling, 
            bool numerical_double_derivative, 
            bool interactions
        );
        arma::vec walk(int MC_cycles, std::string density_filename);
        arma::vec walk(int MC_cycles);
        void minimize_parameters(int MC_cycles, double learning_rate, int max_iter);

        // initialization.cpp
        void set_initial_state_interaction();
        void set_initial_state_no_interaction();
        void set_initial_state();
        void burn_in_interaction();
        void initialize();

        // sampling.cpp
        void importance_sampling_with_interactions(int k);
        void importance_sampling_without_interactions(int k);
        void brute_force_sampling_with_interactions(int k);
        void brute_force_sampling_without_interactions(int k);
        void sampling(int k);
};