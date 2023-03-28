#include "../include/vmc_walker.hpp"
#include "../include/interaction.hpp"
#include "../include/no_interaction.hpp"

using namespace arma;
using namespace std;


void VMCWalker::importance_sampling_with_interactions(int k){
    vec qf = quantum_force_naive(position, relative_position, alpha, beta, hard_core_radius, k);
    // vec qf = quantum_force(position, k, alpha);
    vec gaussian_random = vec(N_dimensions).randn();
    new_position.col(k) = position.col(k) + D * qf * time_step + gaussian_random*sqrt(time_step);
    double rel_pos = 0.0;
    for (int i = 0; i < N_particles; i++ ){
        if (i != k){
            rel_pos = norm(new_position.col(k) - new_position.col(i));
            if (rel_pos < hard_core_radius && hard_core_radius > 0){ //if particles are too close, reject move
                new_position.col(k) = position.col(k); //move particle k back to original position
                relative_position_new = relative_position;
                return; //break out of function
            }
            else{
                relative_position_new(max(i, k), min(i ,k)) = rel_pos;
            }
        }
    }
    vec qf_new = quantum_force_naive(new_position, relative_position_new, alpha, beta, hard_core_radius, k);
    // vec qf_new = quantum_force(new_position, k, alpha);
    double greens = greens_ratio(qf, qf_new, position, new_position, time_step, D, k);
    double ratio = greens*probability_ratio_naive(position, new_position, relative_position, relative_position_new, alpha, beta, hard_core_radius, k);
    //cout << ratio << endl;
    double random_number = rng_double(generator);
    if (random_number <= ratio){ //check whether or not to move particle k
        accepted_moves += 1;
        position.col(k) = new_position.col(k);
        relative_position = relative_position_new;
    }
    else{ //move particle k back to original position
        new_position.col(k) = position.col(k); 
        relative_position_new = relative_position;
    }
}

void VMCWalker::importance_sampling_without_interactions(int k)
{
    vec qf = quantum_force(position, k, alpha);
    vec gaussian_random = vec(N_dimensions).randn();
    new_position.col(k) = position.col(k) + D * qf * time_step + gaussian_random*sqrt(time_step);
    vec qf_new = quantum_force(new_position, k, alpha);
    double greens = greens_ratio(qf, qf_new, position, new_position, time_step, D, k);

    double ratio = greens * probability_ratio(position, new_position, alpha, k);
    // cout << ratio << endl;
    double random_number = rng_double(generator);
    if (random_number <= ratio){ //check whether or not to move particle k
        accepted_moves += 1;
        position.col(k) = new_position.col(k);
    }
    else{ //move particle k back to original position
        new_position.col(k) = position.col(k); 
    }
}

void VMCWalker::brute_force_sampling_with_interactions(int k)
{
    vec random_walk = vec(N_dimensions).randu() - stepping_vector;
    new_position.col(k) = position.col(k) + step*random_walk; //move particle k with a random walk
    double rel_pos = 0.0;
    for (int i = 0; i < N_particles; i++ ){
        if (i != k){
            rel_pos = norm(new_position.col(k) - new_position.col(i));
            if (rel_pos < hard_core_radius && hard_core_radius > 0){ //if particles are too close, reject move
                new_position.col(k) = position.col(k); //move particle k back to original position
                relative_position_new = relative_position;
                return; //break out of function
            }
            else{
                relative_position_new(max(i, k), min(i, k)) = rel_pos;
            }
        }
    }
    double ratio = probability_ratio_naive(position, new_position, relative_position, relative_position_new, alpha, beta, hard_core_radius, k);
    // cout << ratio << endl;
    double random_number = rng_double(generator);
    if (random_number <= ratio){
        accepted_moves += 1;
        //check whether or not to move particle k
        position.col(k) = new_position.col(k);
        relative_position = relative_position_new;
    }
    else{
        new_position.col(k) = position.col(k); //move particle k back to original position
        relative_position_new = relative_position;
    }
}

void VMCWalker::brute_force_sampling_without_interactions(int k)
{
    int N_dimensions = position.n_rows;
    vec random_walk = vec(N_dimensions).randu() - stepping_vector;
    new_position.col(k) = position.col(k) + step*random_walk; //move particle k with a random walk
    double ratio = probability_ratio(position, new_position, alpha, k);
    double random_number = rng_double(generator);
    if (random_number <= ratio){
        accepted_moves += 1;
        //check whether or not to move particle k
        position.col(k) = new_position.col(k);
    }
    else{
        new_position.col(k) = position.col(k); //move particle k back to original position
    }
}

void VMCWalker::sampling(int k)
{
    if (importance_sampling && interactions){
        importance_sampling_with_interactions(k);
    }
    else if (importance_sampling && !interactions){
        importance_sampling_without_interactions(k);
    }
    else if (!importance_sampling && interactions){
        brute_force_sampling_with_interactions(k);
    }
    else if (!importance_sampling && !interactions){
        brute_force_sampling_without_interactions(k);
    }
    else{
        cout << "Error: something went wrong." << endl;
    }
}
