import yaml
import subprocess
import numpy as np

def write_CV_file(US_data, dof, writefile=True):
    """
    Write a suitable array of starting values to CV0_vals.list
    """
    dof_names = {
        'thetaA' : 'theta_A_0',
        'thetaB' : 'theta_B_0',
        'phiA' : 'phi_A_0',
        'phiB' : 'phi_B_0',
        'phiC' : 'phi_C_0'
    }

    eq_value = np.round(US_data['Boresch equilibrium values'][dof_names[dof]], 1)
    starting_value = eq_value - 0.4

    CV_string = str(np.round(starting_value,2))

    for i in range(1, 10): # Write 9 windows per dof
        CV_string += f", {np.round(starting_value + i*0.1, 2)}"

    if writefile == True:
        with open('CV0_vals.list', 'w') as file:
            file.write(CV_string)

    return CV_string

with open('../US_config.yaml', 'r') as file:
    US_data = yaml.safe_load(file)

for dof in ['thetaA', 'thetaB', 'phiA', 'phiB', 'phiC']:

    # Generate CV file
    write_CV_file(US_data, dof)

    # Read in file
    with open("params.in", "r") as file:
        lines = [line.strip() for line in file]

    # Submit triplicate run
    for n_run in [1,2,3]:

        # Write paramfiles and submit triplicate run 
        with open("params.in", "w") as file:
            for n in range(len(lines)):
                if n == 1:
                    file.write(f"run_number = {n_run}\n")
                elif n == 3:
                    file.write(f"dof = {dof}\n")
                else:
                    file.write(lines[n] + "\n")

        try:
            subprocess.run("./submitrun", check=True, shell=True)
        except subprocess.CalledProcessError as e:
            print(f"An error occurred: {e}")