import openmm as mm
import openmm.app as app
import openmm.unit as unit
import sys
import os
from sys import stdout
import numpy as np
import pandas as pd
import yaml

"""r0 value for window bias"""

r0 = np.round(float(sys.argv[1]), 4) # Take r0 value as command line argument
jobid = int(sys.argv[2])

"""Read in US configuration from yaml"""

with open('../US_config.yaml', 'r') as file:
    US_data = yaml.safe_load(file)

"""Read global params from params.in"""

def read_param(param_str, jobid):
    """
    Read in a specific parameter and assign the parameter value to a variable
    """
    with open(f'jobs/{jobid}.in', 'r') as file:
        for line in file:
            line = line.strip()

            if line.startswith(param_str):
                parts = line.split(' = ')
                part = parts[1].strip()

    values = part.split(' ')
    value = values[0].strip()

    # Attempt int conversion
    try:
        value = int(value)
    except:
        value = str(value)

    return value

# MD parameters
timestep = read_param('timestep', jobid)
save_traj = read_param('save_traj', jobid)
equil_steps = int(read_param('equilibration_time', jobid)//(timestep*1e-6))
sampling_steps = int(read_param('sampling_time', jobid)//(timestep*1e-6))
record_steps = read_param('n_steps_between_sampling', jobid)

# Force constants
k_Boresch = read_param('k_Boresch', jobid)
k_RMSD = read_param('k_RMSD', jobid)
k_sep = read_param('k_sep', jobid)

# Directory to save all results
run_number = read_param('run_number', jobid)
restraint_type = read_param('RMSD_restraint', jobid)
inputdir = f"windows/{r0}"
savedir = f"results/{restraint_type}_RMSD/run{run_number}"

if not os.path.exists(inputdir): # Check if directory exists
    raise FileNotFoundError(f"Input directory does not exist for r0 = {r0}")

if os.path.exists(f'{savedir}/{r0}.txt'): # Check if a CV sample file already exists
    raise FileExistsError(f"A file of CV samples already exists for r0 = {r0}")

if not os.path.exists(savedir): # Make save directory if it doesn't yet exist
    os.makedirs(savedir)

if not os.path.exists('jobs'): # Make jobs directory if it doesn't yet exist
    os.makedirs('jobs')

if not os.path.exists('outfiles'): # Make outfiles directory if it doesn't yet exist
    os.makedirs('outfiles')

"""Selection tuple for restraint type"""

if restraint_type == 'backbone':
    restraint_selection = ['C', 'N', 'CA']

elif restraint_type == 'CA':
    restraint_selection = ['CA']

elif restraint_type != 'heavy_atom':
    raise ValueError('Select one of the following restraint type options:  heavy_atom, backbone, CA')

"""System setup"""

dt = timestep*unit.femtoseconds 

# Load param and coord files
prmtop = app.AmberPrmtopFile(f'../equilibrated_structures/complex_eq.prmtop')

system = prmtop.createSystem(nonbondedMethod=app.PME, hydrogenMass=1.5*unit.amu, nonbondedCutoff=1.0*unit.nanometer, constraints=app.HBonds)  
integrator = mm.LangevinMiddleIntegrator(0.0000*unit.kelvin, 1.0000/unit.picosecond, dt)

simulation = app.Simulation(prmtop.topology, system, integrator)

# Set positions to the correct window
pdb = app.PDBFile(f"windows/{r0}/{r0}.pdb")
simulation.context.setPositions(pdb.positions)

# Add reporters to output data
simulation.reporters.append(app.StateDataReporter(f'{savedir}/{r0}.csv', 1000, step=True, time=True, potentialEnergy=True, kineticEnergy=True, totalEnergy=True, temperature=True, volume=True, density=True, speed=True))
simulation.reporters.append(app.StateDataReporter(stdout, 2000, step=True, time=True, potentialEnergy=True, temperature=True, speed=True))

if save_traj=='True':
    simulation.reporters.append(app.DCDReporter(f'{savedir}/{r0}.dcd', 2000))

# Minimise energy 
simulation.minimizeEnergy()

"""System heating"""

for i in range(50):
    integrator.setTemperature(6*(i+1)*unit.kelvin)
    simulation.step(1000)

simulation.step(1000)
simulation.context.setVelocitiesToTemperature(300.0000*unit.kelvin)

"""Find indices of all small molecule heavy atoms"""

mol_indices = []
for atom in simulation.topology.atoms():
    if atom.residue.name == 'MOL':
        if not atom.name.startswith('H'):
            mol_indices.append(atom.index)

"""RMSD Restraints"""

reference_positions = simulation.context.getState(getPositions=True).getPositions()
# 0-indexed residues
receptor_residues = US_data['Receptor residues']
ligand_residues = US_data['Ligand residues']

receptor_atoms = [
    atom.index for atom in simulation.topology.atoms()
    if atom.residue.index in receptor_residues and atom.name in restraint_selection
]
ligand_atoms = [
    atom.index for atom in simulation.topology.atoms()
    if atom.residue.index in ligand_residues and atom.name in restraint_selection
]

# Add small molecule heavy atom indices to restrained ligand atoms
receptor_atoms = receptor_atoms + mol_indices

# Add restraining forces for receptor and ligand rmsd
receptor_rmsd_force = mm.CustomCVForce('0.5*k_rec*rmsd^2')
receptor_rmsd_force.addGlobalParameter('k_rec', k_RMSD * unit.kilocalories_per_mole / unit.angstrom**2)
receptor_rmsd_force.addCollectiveVariable('rmsd', mm.RMSDForce(reference_positions, receptor_atoms))
system.addForce(receptor_rmsd_force)

ligand_rmsd_force = mm.CustomCVForce('0.5*k_lig*rmsd^2')
ligand_rmsd_force.addGlobalParameter('k_lig', k_RMSD * unit.kilocalories_per_mole / unit.angstrom**2)
ligand_rmsd_force.addCollectiveVariable('rmsd', mm.RMSDForce(reference_positions, ligand_atoms))
system.addForce(ligand_rmsd_force)

simulation.context.reinitialize(preserveState=True)

for atom in simulation.topology.atoms():
    if atom.index==receptor_atoms[-1] and atom.residue.name!='MOL':
        raise ValueError(f'Small molecule is missing!')
    
"""Radial separation CV"""

# 1-indexing from MDAnalysis
rec_interface_res = US_data['Receptor interface residues']
lig_interface_res =  US_data['Ligand interface residues']

# Account for OpenMM residue 0-indexing
rec_interface_res = -1 + np.array(rec_interface_res) 
lig_interface_res = -1 + np.array(lig_interface_res)

rec_group = [
    atom.index for atom in simulation.topology.atoms()
    if atom.residue.index in rec_interface_res and atom.name=='CA'
]

lig_group = [
    atom.index for atom in simulation.topology.atoms()
    if atom.residue.index in lig_interface_res and atom.name=='CA'
]

# Add small molecule heavy atom indices to restrained ligand atoms
rec_group = rec_group + mol_indices

# Define radial distance as collective variable which we will vary
cv = mm.CustomCentroidBondForce(2, "distance(g1,g2)")
cv.addGroup(np.array(rec_group))
cv.addGroup(np.array(lig_group))

# Specify bond groups
bondGroups = [0, 1]
cv.addBond(bondGroups)

# Define biasing potential
bias_pot = mm.CustomCVForce('0.5 * k_r * (cv-r0)^2')
bias_pot.addGlobalParameter('k_r', k_sep * unit.kilocalories_per_mole / unit.angstrom**2)
bias_pot.addGlobalParameter('r0', r0* unit.nanometers)

bias_pot.addCollectiveVariable('cv', cv)
system.addForce(bias_pot)

simulation.context.reinitialize(preserveState=True)

"""Boresch restraints"""

def obtain_CA_idx(res_idx):

    """Function to obtain the index of the alpha carbon for a given residue index"""
    
    atom_idx = None

    for atom in simulation.topology.atoms():
        if atom.residue.index == res_idx and atom.name=='CA':
            atom_idx = atom.index
    
    return atom_idx
    
# Boresch_residues = [13, 8, 142, 370, 332, 244]

# Define anchor points (1-indexing)
res_b = US_data['Boresch anchor points']['res_b']
res_c = US_data['Boresch anchor points']['res_c']
res_B = US_data['Boresch anchor points']['res_B']
res_C = US_data['Boresch anchor points']['res_C']

# Account for OpenMM 0-indexing 
res_b -=1 
res_c -=1
res_B -=1
res_C -=1

# Find atomic indices
idx_b = obtain_CA_idx(res_b)
idx_c = obtain_CA_idx(res_c)
idx_B = obtain_CA_idx(res_B)
idx_C = obtain_CA_idx(res_C)

print('Anchor points:')
for atom in simulation.topology.atoms():
    if atom.index in [idx_b, idx_c, idx_B, idx_C]:
        print(atom)

# Check that we have only selected CA anchor points
all_atoms = [idx_b] + [idx_c] + lig_group + [idx_B] + [idx_C]
for atom in simulation.topology.atoms():
    if atom.index in all_atoms and atom.name != 'CA':
        raise ValueError('Select only CA atoms as anchorpoints')
    
# Equilibrium values of Boresch dof
theta_A_0 = US_data['Boresch equilibrium values']['theta_A_0']
theta_B_0 = US_data['Boresch equilibrium values']['theta_B_0']
phi_A_0 = US_data['Boresch equilibrium values']['phi_A_0']
phi_B_0 = US_data['Boresch equilibrium values']['phi_B_0']
phi_C_0 = US_data['Boresch equilibrium values']['phi_C_0']

k_Boresch = k_Boresch * unit.kilocalories_per_mole / unit.radians**2 #Set global force constant

theta_A_pot = mm.CustomCentroidBondForce(3, '0.5 * k_Boresch * (angle(g1,g2,g3)-theta_A_0)^2')
theta_A_pot.addGlobalParameter('theta_A_0', theta_A_0)
theta_A_pot.addGlobalParameter('k_Boresch', k_Boresch)

# Add the particle groups
theta_A_pot.addGroup([idx_b])
theta_A_pot.addGroup(np.array(rec_group))
theta_A_pot.addGroup(np.array(lig_group))

# Add the centroid angle bond
theta_A_pot.addBond([0, 1, 2])

system.addForce(theta_A_pot)

theta_B_pot = mm.CustomCentroidBondForce(3, '0.5 * k_Boresch * (angle(g1,g2,g3)-theta_B_0)^2')
theta_B_pot.addGlobalParameter('theta_B_0', theta_B_0)
theta_B_pot.addGlobalParameter('k_Boresch', k_Boresch)

# Add the particle groups
theta_B_pot.addGroup(np.array(rec_group))
theta_B_pot.addGroup(np.array(lig_group))
theta_B_pot.addGroup([idx_B])

# Add the centroid angle bond
theta_B_pot.addBond([0, 1, 2])

system.addForce(theta_B_pot)

phi_A_pot = mm.CustomCentroidBondForce(4, "0.5*k_Boresch*min(dtheta, 2*pi-dtheta)^2; dtheta = abs(dihedral(g1,g2,g3,g4)-phi_A_0); pi = 3.1415926535")
phi_A_pot.addGlobalParameter('phi_A_0', phi_A_0)
phi_A_pot.addGlobalParameter('k_Boresch', k_Boresch)

# Add the particle groups
phi_A_pot.addGroup([idx_c])
phi_A_pot.addGroup([idx_b])
phi_A_pot.addGroup(np.array(rec_group))
phi_A_pot.addGroup(np.array(lig_group))

# Add the centroid angle bond
phi_A_pot.addBond([0, 1, 2, 3])

system.addForce(phi_A_pot)

phi_B_pot = mm.CustomCentroidBondForce(4, "0.5*k_Boresch*min(dtheta, 2*pi-dtheta)^2; dtheta = abs(dihedral(g1,g2,g3,g4)-phi_B_0); pi = 3.1415926535")
phi_B_pot.addGlobalParameter('phi_B_0', phi_B_0)
phi_B_pot.addGlobalParameter('k_Boresch', k_Boresch)

# Add the particle groups
phi_B_pot.addGroup([idx_b])
phi_B_pot.addGroup(np.array(rec_group))
phi_B_pot.addGroup(np.array(lig_group))
phi_B_pot.addGroup([idx_B])

# Add the centroid angle bond
phi_B_pot.addBond([0, 1, 2, 3])

system.addForce(phi_B_pot)

phi_C_pot = mm.CustomCentroidBondForce(4, "0.5*k_Boresch*min(dtheta, 2*pi-dtheta)^2; dtheta = abs(dihedral(g1,g2,g3,g4)-phi_C_0); pi = 3.1415926535")
phi_C_pot.addGlobalParameter('phi_C_0', phi_C_0)
phi_C_pot.addGlobalParameter('k_Boresch', k_Boresch)

# Add the particle groups
phi_C_pot.addGroup(np.array(rec_group))
phi_C_pot.addGroup(np.array(lig_group))
phi_C_pot.addGroup([idx_B])
phi_C_pot.addGroup([idx_C])

# Add the centroid angle bond
phi_C_pot.addBond([0, 1, 2, 3])

system.addForce(phi_C_pot)

simulation.context.reinitialize(preserveState=True)

"""Collecting CV samples"""

def obtain_RMSDs(simulation):
    """
    Return the instantaneous values of the RMSDs restrained
    by the previously applied bias potentials
    """
    # Convert RMSD units to Angstrom
    RMSD_rec = 10*receptor_rmsd_force.getCollectiveVariableValues(simulation.context)[0]
    RMSD_lig = 10*ligand_rmsd_force.getCollectiveVariableValues(simulation.context)[0]

    return [RMSD_rec, RMSD_lig]

print('running window for r0 = ', r0)

# Equilibration if specified
if equil_steps>0:
    simulation.step(equil_steps) 

# Prepare sampling
n_samples = int(sampling_steps//record_steps) # Total number of samples
cv_values = np.zeros((n_samples, 2)) # Empty array to store samples
dof_data = np.zeros((n_samples, 4))

# Run the simulation and record the value of the CV.
for i in range(n_samples):

    simulation.step(record_steps)

    # get the current value of the cv
    current_cv_value = bias_pot.getCollectiveVariableValues(simulation.context)
    sample = current_cv_value[0]
    cv_values[i] = [i, sample]

    # Save the other dofs
    dofs = obtain_RMSDs(simulation)
    dof_data[i] = [i, i*timestep*1e-6] + dofs

# Final save
np.savetxt(f'{savedir}/{r0}.txt', cv_values)

# Save the RMSDs to csv
df = pd.DataFrame(data=dof_data, columns=[
    'Steps',
    'Time (ns)',
    'RMSD_rec',
    'RMSD_lig'
])

df.to_csv(f'{savedir}/{r0}_RMSD.csv', index=False)

print('Completed window for r0 = ', r0)




