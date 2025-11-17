import subprocess

systems = ['receptor_glue', 'ligand', 'receptor_glue_only', 'ligand_only', 'receptor_gluewithligand', 'ligandwithreceptor_glue'] # List of all systems

for species in systems:

    for n_run in [1,2,3]:

        # Read in file
        with open("params.in", "r") as file:
            lines = [line.strip() for line in file]

        # Write paramfiles and submit triplicate run 
        with open("params.in", "w") as file:
            for n in range(len(lines)):
                if n == 1:
                    file.write(f"run_number = {n_run}\n")
                elif n == 3:
                    file.write(f"species = {species}\n")
                else:
                    file.write(lines[n] + "\n")

        try:
            subprocess.run("./submitrun", check=True, shell=True)
        except subprocess.CalledProcessError as e:
            print(f"An error occurred: {e}")






    