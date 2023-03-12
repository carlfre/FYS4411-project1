#include "../include/no_interaction.hpp"

using namespace arma;
using namespace std;

double local_energy(mat& position, double alpha){
    //k is the index of the particle we are calculating the energy for
    double energy = (0.5 - 2*alpha*alpha)
                        *accu(position % position) 
                        + position.n_cols*position.n_rows*alpha;
    //accu(A%A) is elementwise multiplication and sum over all elements
    return energy;
}

double probability_ratio(mat& old_position, mat& new_position, double alpha, int particle_index){
    double exponent = dot(old_position.col(particle_index), old_position.col(particle_index)) 
                    - dot(new_position.col(particle_index), new_position.col(particle_index));

    return exp(2*alpha*exponent);
}

double trial_wavefunction(mat& position, double alpha){
    double exponent = accu(position%position);
    //accu(A%A) is elementwise multiplication and sum over all elements
    return exp(-alpha*exponent);
}

double dwf_dalpha(arma::mat& position){
    return -accu(position % position);
}

double local_energy_numerical(mat& position, double alpha, double h){
    double energy = 0;
    //loop over all particles
    for (int i = 0; i < position.n_cols; i++){
        //loop over all dimensions
        for (int j = 0; j < position.n_rows; j++){
            energy += -0.5*numerical_double_derivative(position, alpha, h, j, i);
        }
    }
    energy /= trial_wavefunction(position, alpha);
    //potential energy
    energy += 0.5*accu(position % position);
    return energy;
}

double numerical_double_derivative(mat& position, double alpha, double h, int dimension_index, int particle_index){
    double wf_mid = trial_wavefunction(position, alpha);
    //move particle in dimension dimension_index h forwards
    position(dimension_index, particle_index) += h;
    double wf_plus = trial_wavefunction(position, alpha);
    //move particle in dimension dimension_index h backwards
    position(dimension_index, particle_index) -= 2*h;
    double wf_minus = trial_wavefunction(position, alpha);
    //move particle back to original position
    position(dimension_index, particle_index) += h;
    double double_derivative = (wf_plus - 2*wf_mid + wf_minus)/(h*h);
    return double_derivative;
}

vec quantum_force(mat& position, int index, double alpha){
    // expression: 2* grad psi_T / psi_T
    // for non-interactive case, simplifies to -4*alpha*position
    // set force for particle index
    vec force = -4*alpha*position.col(index);
    return force;
}

double greens_ratio(vec& qf_old, vec& qf_new, mat& position_old, mat& position_new, double time_step, double D, int index){
    // print dimensions of matrix
    vec tmp = 0.5*(0.5*D*time_step*(qf_old - qf_new) + position_old.col(index) - position_new.col(index));
    double exponent = accu(tmp%(qf_old + qf_new)); //elementwise multiplication
    
    return exp(exponent);
}