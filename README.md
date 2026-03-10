# Umbrella Sampling Ternary Complex Free Energy Calculations

This repository contains a skeleton template for calculating the binding free energy of a ternary complex. The methodology is based on the geometric route of [Woo and Roux ](https://www.pnas.org/doi/10.1073/pnas.0409005102), with [Boresch-style restraints](https://pubs.acs.org/doi/10.1021/jp0217839) being used to restrain the relative orientation of the two proteins, as per the framework of [Notari *et al*](https://pubs.acs.org/doi/10.1021/acs.jctc.4c01695). Umbrella Sampling (US) simulations are used to calculate the Potential of Mean Force (PMF) curves corresponding to each stage of the binding/unbinding process. The MD engine used for all simulations is [OpenMM](https://openmm.org/).

The protocol follows three main stages: 

1. **Unrestrained MD** for obtaining equilibrated starting structures and equilibrium values for the various Boresch DOFs
2. **Steered MD** for generating input configurations for the separation US simulations
3. **Umbrella Sampling** simulations for calculating the various free energy contributions required for calculating the overall binding free energy

## Requirements for running the protocol

- PyMol for editing and generating input structures
- [Ambertools](https://anaconda.org/conda-forge/ambertools) for generating input topologies and coordinates
- Python environment with openmm<8.3, pandas, pyyaml and MDAnalysis

## Running unrestrained MD

1. Prepare your input files as **complex.prmtop**/**complex.inpcrd** within the ```UnrestrainedMD/structures``` directory,  e.g. using ```tleap -f tleap.in```. **Important** - assign the residue name "MOL" to the small molecule when generating input parameters in antechamber. The input files should be equilibrated under NPT conditions. The script ```scripts/NPT_equil.py``` can be used to perform this equilibration if desired.
2. Specify the receptor and ligand residue indices in config.yaml using **0-indexing**. This is needed to specify the receptor and ligand selections in the analysis scripts.
3. Submit an MD simulation to slurm using ```sbatch submitMD.sh -r $RUN_NUMBER```. It is recommended to perform three runs.
4. Once the MD simulations have finished, calculate the RMSD and RMSF for the trajectories by running ```analyse_RMSF.py```. By default, this script assumes there are three runs of unrestrained MD, you can change this if necessary.
5. Visualise the RMSF in the notebook ```select_anchors.ipynb``` and pick stable residues which will become the anchor points for the Boresch angles and dihedrals. Ensure that the anchor point CA indices have been written to ```config.yaml```
6. Run the script ```generate_yaml.py``` to extract the interface residues, interface equilbration time and Boresch DOF equilibrium values. These values are written to ```../umbrella_sampling/US_config.yaml```. Plots of the Boresch DOF histograms and the interface RMSD are written to the ```plots``` folder. A pdb snapshot of the equilibrated reference frame is written to the ```equilibrated_structures``` directory in both the ```Unrestrained_MD``` and ```umbrella_sampling``` locations. By default, the snapshot is saved from run 1 of the unrestrained MD simulations. Within this repository, the ```equilibrated_structures``` has been populated with example inputs for the CRBN-lenalidomide-CK1a ternary complex (PDB [5FQD](https://www.rcsb.org/structure/5FQD)).

## Running Steered MD

The SMD simulation is performed in the ```umbrella_sampling/separation``` directory.
During the analysis of the unrestained MD simulation, a snapshot of the equilibrated system is saved to ```umbrella_sampling/equilibrated_structures```. You have to generate Amber input files from this snapshot under the name ```complex_eq.prmtop/.inpcrd```. See [Amber LEaP tutorial](https://ambermd.org/tutorials/pengfei/index.php) and [GAFF tutorial](https://ambermd.org/tutorials/basic/tutorial4b/index.php) for help. The script ```scripts/NPT_equil.py``` can be used to perform for NPT equilibration if desired.

**Important** - assign the residue name "MOL" to the small molecule when generating input parameters in antechamber. This allows ```SMD.py``` and ```run_window.py``` to add the ligand into the receptor interface.

The steered MD simulation is submitted using the ```SMD.py``` script. You can play around with this script, adding an equilibration period or modifying the force constants if necessary.

## Submitting jobs to slurm

The ```scripts``` directory contains slurm/grid engine submission scripts for the US simulations. Modify ```submit_window_workstation.sh``` to match the slurm configuration on your local workstation. 

## Generating the separation PMF

An Umbrella Sampling simulation can be submitted using the ```submitrun``` executable. This executable calls the script ```run_window.py``` to submit each the simulation for each umbrella. Note that if you want to perform the resolvation step for each window, which can increase simulation speed by $\sim 20$ %, you will need to modify ```windows/modify_pdb.pdb``` to match your specific system. To then use your resolvated ```system.prmtop```/```system.inpcrd``` input files, use the ```run_window_resolvated.py``` script.

Before submitting the ```submitrun``` executable, modify the ```params.in``` file to specify all parameters (force constants, sampling time, etc.) You can specify the window spacing by providing comma separated values in the ```r0_vals.list``` file. 

## RMSD and Boresch US

- Amber .prmtop/.inpcrd inputs should be prepared for the isolated ligand and glue/PROTAC-bound receptor within the ```equilibrated_structures``` directory. These inputs should be named ```ligand.prmtop```/```receptor.prmtop``` respectively.
- US simulations are submitted in the same way as for the separation PMF (i.e. via the ```submitrun``` executable that reads in the simulation parameters from ```params.in```/```r0_vals.list```). For convenience, a python script called ```submit_triplicate.py``` can be used to submit triplicate US simulations for every system (Boresch and RMSD US).

## Analysis

```umbrella_sampling/analysis.py``` contains useful functions for generating PMFs using WHAM, analysing the individual US windows and calculating the corresponding free energy contributions for the various stages of the thermodynamic cycle.

The notebook ```umbrella_sampling/analysis.ipynb``` contains some useful analysis code, including visualisation of average PMFs and automated calculation of free energy contributions.


