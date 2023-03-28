import numpy as np
import subprocess

N_vals = [10, 50, 100]
alpha_vals = [0.49745384375, 0.48920212500000004, 0.48308106250000005]
MC_cycles = [24, 21, 19]

for N, alpha, MC in zip(N_vals, alpha_vals, MC_cycles):
    msg = f"""# alpha
{alpha}
# beta 
2.82843
# gamma
2.82843
# step
1.0
# time_step
0.1
# hard_core_radius
0.0043
# N_particles
{N}
# N_dimensions
3
# importance_sampling
1
# interactions
1
# log_2 (MC_cycles) - 
{MC}
# N_walkers
16"""
    print(N, alpha, MC)


    # write to file
    with open("configs/density.txt", "w") as f:
        f.write(msg)
    
    # run program
    subprocess.run(["./main.exe", "density"])