#include "../include/sampling.hpp"
#include "../include/no_interaction.hpp"
#include "../include/interaction.hpp"

using namespace arma;
using namespace std;


void importance_sampling_with_interactions(
    mat& position, 
    mat& new_position, 
    mat& relative_position, 
    mat& relative_position_new,
    double alpha,
    double beta,
    double hard_core_radius,
    double time_step,
    double D,
    double random_number,
    int k,
    int& accepted_moves)
{
    int N_dimensions = position.n_rows;
    int N_particles = position.n_cols;

    vec qf = quantum_force_naive(position, relative_position, alpha, beta, hard_core_radius, k);
    vec gaussian_random = vec(N_dimensions).randn();
    new_position.col(k) = position.col(k) + D * qf * time_step + gaussian_random*sqrt(time_step);
    double rel_pos = 0.0;
    for (int i = 0; i < N_particles; i++ ){
        if (i != k){
            rel_pos = norm(new_position.col(k) - new_position.col(i));
            if (rel_pos < hard_core_radius){ //if particles are too close, reject move
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
    double greens = greens_ratio(qf, qf_new, position, new_position, time_step, D, k);
    double ratio = greens * probability_ratio_naive(position, new_position, relative_position, relative_position_new, alpha, beta, hard_core_radius, k);
        
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

void importance_sampling_without_interactions(
    mat& position, 
    mat& new_position, 
    double alpha,
    double beta,
    double time_step,
    double D,
    double random_number,
    int k,
    int& accepted_moves)
{
    int N_dimensions = position.n_rows;
    //importance sampling for non-interactive case
    vec qf = quantum_force(position, k, alpha);
    vec gaussian_random = vec(N_dimensions).randn();
    new_position.col(k) = position.col(k) + D * qf * time_step + gaussian_random*sqrt(time_step);

    vec qf_new = quantum_force(new_position, k, alpha);
    double greens = greens_ratio(qf, qf_new, position, new_position, time_step, D, k);
    double ratio = greens * probability_ratio(position, new_position, alpha, k);

    if (random_number <= ratio){ //check whether or not to move particle k
        accepted_moves += 1;
        position.col(k) = new_position.col(k);
    }
    else{ //move particle k back to original position
        new_position.col(k) = position.col(k); 
    }
}

void brute_force_sampling_with_interactions(    
    mat& position, 
    mat& new_position, 
    mat& relative_position, 
    mat& relative_position_new,
    vec& stepping_vector,
    double alpha,
    double beta,
    double hard_core_radius,
    double step,
    double random_number,
    int k,
    int& accepted_moves)
{
    int N_dimensions = position.n_rows;
    int N_particles = position.n_cols;

    vec random_walk = vec(N_dimensions).randu() - stepping_vector;
    new_position.col(k) = position.col(k) + step*random_walk; //move particle k with a random walk
    double rel_pos = 0.0;
    for (int i = 0; i < N_particles; i++ ){
        if (i != k){
            rel_pos = norm(new_position.col(k) - new_position.col(i));
            if (rel_pos < hard_core_radius){ //if particles are too close, reject move
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

void brute_force_sampling_without_interactions(
    mat& position, 
    mat& new_position, 
    vec& stepping_vector,
    double alpha,
    double step,
    double random_number,
    int k,
    int& accepted_moves)
{
    int N_dimensions = position.n_rows;
    vec random_walk = vec(N_dimensions).randu() - stepping_vector;
    new_position.col(k) = position.col(k) + step*random_walk; //move particle k with a random walk
    double ratio = probability_ratio(position, new_position, alpha, k);
    if (random_number <= ratio){
        accepted_moves += 1;
        //check whether or not to move particle k
        position.col(k) = new_position.col(k);
    }
    else{
        new_position.col(k) = position.col(k); //move particle k back to original position
    }
}

void sampling(mat& position, mat& new_position, double alpha, int k, double step, double time_step, double D, bool importance_sampling, int& accepted_moves, vec& stepping_vector, double random_number, bool interactions, mat& relative_position, mat& relative_position_new, double gamma, double beta, double hard_core_radius){
    if (importance_sampling && interactions){
        importance_sampling_with_interactions(
            position, 
            new_position, 
            relative_position, 
            relative_position_new,
            alpha,
            beta,
            hard_core_radius,
            time_step,
            D,
            random_number,
            k,
            accepted_moves);  
    }
    else if (importance_sampling && !interactions){
        importance_sampling_without_interactions(
            position, 
            new_position, 
            alpha,
            beta,
            time_step,
            D,
            random_number,
            k,
            accepted_moves);
    }
    else if (!importance_sampling && interactions){
        brute_force_sampling_with_interactions(
            position, 
            new_position, 
            relative_position, 
            relative_position_new,
            stepping_vector,
            alpha,
            beta,
            hard_core_radius,
            step,
            random_number,
            k,
            accepted_moves);
    }
    else if (!importance_sampling && !interactions){
        brute_force_sampling_without_interactions(
            position, 
            new_position, 
            stepping_vector,
            alpha,
            step,
            random_number,
            k,
            accepted_moves);
    }
}