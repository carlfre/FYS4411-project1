import numpy as np
import subprocess



N_values = [10, 50, 100]

for N in N_values:
    print("run for N =", N)
    msg=f"""# alpha_0 (initial guess for alpha)
0.5
# beta
2.82843
#gamma
2.82843
# step
0.3
# time_step
0.3
# hard_core_radius
0.0043
# N_particles
{N}
# N_dimensions
3
# importance_sampling
1
# log_10(MC_cycles)
4
# log_10(Learning rate for gradient descent)
0.02
# Max iterations for gradient descent
50"""


    # write to config file
    with open("configs/interactions_gradient.txt", "w") as f:
        f.write(msg)
    
    # run program
    subprocess.run(["./main.exe", "interactions_gradient"])

    print("Finished run for N =", N)