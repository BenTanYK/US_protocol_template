"""
Calculate the interface equilibration time, Boresch equilibrium values 
and interface residue indices. These parameters are written to the 
US_config.yaml file in the ../umbrella_sampling directory
"""

import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import rms, align
import matplotlib.pyplot as plt
from MDAnalysis.analysis.distances import dist
from MDAnalysis.transformations.nojump import NoJump
from tqdm import tqdm
import pandas as pd
import pickle
import red
import os
import yaml

def obtain_RMSD(run_numbers, species='complex'):

    if species == 'complex':
        filename = 'RMSD.csv'
    elif species == 'ligand':
        filename = 'RMSD_lig.csv'
    else:
        filename = 'RMSD_rec.csv'

    if isinstance(run_numbers, int):
        run_numbers = [run_numbers]

    RMSDs = []

    for run_number in run_numbers:
        df = pd.read_csv(f"results/run{run_number}/{filename}")
        time = df['Time (ns)'].to_numpy()
        RMSDs.append(df['RMSD (Angstrom)'].to_numpy())

    RMSDs = np.array(RMSDs)

    return time, np.average(RMSDs, axis=0)

def obtain_RMSF(run_number, species='complex'):

    if species == 'complex':
        filename = 'RMSF.csv'
    elif species == 'ligand':
        filename = 'RMSF_lig.csv'
    else:
        filename = 'RMSF_rec.csv'

    df = pd.read_csv(f"results/run{run_number}/{filename}")

    residx = df['Residue index'].to_numpy()
    RMSF = df['RMSF (Angstrom)'].to_numpy()

    return residx, RMSF

def obtain_av_RMSF(run_numbers, species='complex'):

    RMSFs = []

    for run_number in run_numbers:

        res, rmsf = obtain_RMSF(run_number, species)
        RMSFs.append(rmsf)

    RMSFs = np.array(RMSFs)

    return res, np.average(RMSFs, axis=0), np.std(RMSFs, axis=0)

def obtain_CA_idx(u, res_idx):
    """Function to obtain the index of the alpha carbon for a given residue index"""
    
    selection_str = f"protein and resid {res_idx} and name CA"
    
    selected_CA = u.select_atoms(selection_str)

    if len(selected_CA.indices) == 0:
        print('CA not found for the specified residue...')
    
    elif len(selected_CA.indices) > 1:
        print('Multiple CAs found, uh oh...')

    else:  
        return selected_CA.indices[0]

def closest_residue_to_point(atoms, point):
    """Find the closest residue in a selection of atoms to a given point"""
    residues = atoms.residues
    distances = np.array([np.linalg.norm(res.atoms.center_of_mass() - point) for res in residues])

    # Find the index of the smallest distance
    closest_residue_index = np.argmin(distances)

    # Return the closest residue
    return residues[closest_residue_index], distances[closest_residue_index]

def getDistance(idx1, idx2, u):
    """
    Get the distance between two atoms in a universe.

    Parameters
    ----------
    idx1 : int
        Index of the first atom
    idx2 : int
        Index of the second atom
    u : MDAnalysis.Universe
        The MDA universe containing the atoms and
        trajectory.

    Returns
    -------
    distance : float
        The distance between the two atoms in Angstroms.
    """
    distance = dist(
        mda.AtomGroup([u.atoms[idx1]]),
        mda.AtomGroup([u.atoms[idx2]]),
        box=u.dimensions,
    )[2][0]
    return distance

def obtain_angle(pos1, pos2, pos3):

    return mda.lib.distances.calc_angles(pos1, pos2, pos3)

def obtain_dihedral(pos1, pos2, pos3, pos4):
    
    return mda.lib.distances.calc_dihedrals(pos1, pos2, pos3, pos4)

def generate_selection_str(indices):
    """
    Generate the MDAnalysis string for a given array of indices
    """
    selection_str = "resid "

    for index in indices:
        selection_str+=f"{index} "

    selection_str+='and name N CA C'

    return selection_str

def obtain_Boresch_dof(dof, u, rec_group, lig_group, res_b, res_c, res_B, res_C):

    """
    Calculate a Boresch DOF (thetaA, thetaB, phiA, phiB, phiC) 
    Specify the recepter and ligand interface as rec_group and lig_group
    The other anchor points are given as res_b, res_c, res_B and res_C
    """

    group_a = u.atoms[rec_group]
    group_b = u.atoms[[obtain_CA_idx(u, res_b)]]
    group_c = u.atoms[[obtain_CA_idx(u, res_c)]]
    group_A = u.atoms[lig_group]
    group_B = u.atoms[[obtain_CA_idx(u, res_B)]]
    group_C = u.atoms[[obtain_CA_idx(u, res_C)]]

    pos_a = group_a.center_of_mass()
    pos_b = group_b.center_of_mass()
    pos_c = group_c.center_of_mass()
    pos_A = group_A.center_of_mass()
    pos_B = group_B.center_of_mass()
    pos_C = group_C.center_of_mass()

    dof_indices = {
        'thetaA' : [pos_b, pos_a, pos_A],
        'thetaB' : [pos_a, pos_A, pos_B],
        'phiA' : [pos_c, pos_b, pos_a, pos_A],
        'phiB': [pos_b, pos_a, pos_A, pos_B],
        'phiC': [pos_a, pos_A, pos_B, pos_C]
    }

    indices = dof_indices[dof]

    if len(indices) == 3:
        return obtain_angle(indices[0], indices[1], indices[2])

    else:
        return obtain_dihedral(indices[0], indices[1], indices[2], indices[3])

def obtain_interface_RMSD(run_number, rec_interface_res, lig_interface_res):
    """
    Obtain the interface RMSD 
    """

    # Add lenalidomide to interface
    len_selection = 'resname MOL and not name H*'
    
    interface_res = np.append(rec_interface_res, lig_interface_res)
    interface_selection_str = generate_selection_str(interface_res)

    prmtop = "structures/complex.prmtop"
    dcd = f"results/run{run_number}/traj.dcd"

    u = mda.Universe(prmtop, dcd)

    dry_system = u.select_atoms("not resname HOH or resname WAT")

    ref = dry_system
    R_u =rms.RMSD(dry_system, ref, select=f"({interface_selection_str}) or ({len_selection})")
    R_u.run()

    rmsd_u = R_u.rmsd.T #take transpose
    time = rmsd_u[1]/1000 # Units of ns
    rmsd = rmsd_u[2]

    return time, rmsd

def obtain_av_interface_RMSD(run_numbers, rec_interface_res, lig_interface_res):
    """
    Obtain average interface RMSD
    """ 

    RMSDs = []

    for run_number in run_numbers:

        time, rmsd = obtain_interface_RMSD(run_number, rec_interface_res, lig_interface_res)
        RMSDs.append(rmsd)

    RMSDs = np.array(RMSDs)

    return time, np.average(RMSDs, axis=0)

"""Read in anchor points and selection bounds"""
print("Reading in residue selections and anchor points\n")

with open('config.yaml', 'r') as infile:
    initial_config = yaml.safe_load(infile)

# Boresch anchor points
res_b = initial_config['Boresch anchor points']['res_b']
res_c = initial_config['Boresch anchor points']['res_c']
res_B = initial_config['Boresch anchor points']['res_B']
res_C = initial_config['Boresch anchor points']['res_C']

boresch_anchor_points = {
    'res_b' : res_b,
    'res_c' : res_c,
    'res_B' : res_B,
    'res_C' : res_C
}

# 0-indexing
rec_start = initial_config['Receptor residue range'][0]
rec_end = initial_config['Receptor residue range'][-1]
lig_start = initial_config['Ligand residue range'][0]
lig_end = initial_config['Ligand residue range'][-1]

receptor_res = [int(res) for res in range(rec_start, rec_end+1)] # no lenalidomide
ligand_res = [int(res) for res in range(lig_start, lig_end+1)]

"""Initialise dictionary to store parameters"""
print("Initialising dictionary to store config params...\n")

config_data = {
    'Boresch anchor points' : boresch_anchor_points,
    'Receptor residues' : receptor_res,
    'Ligand residues' : ligand_res
}

"""Calculating initial interface"""
print('Calculating the initial interface for RMSD calculation...\n')

# Identify interface from starting structure

u_ref = mda.Universe('structures/complex.prmtop', f'structures/complex.inpcrd')

receptor_selection_str=f"resid {rec_start+1}-{rec_end+1} and name CA" # receptor
ligand_selection_str = f"resid {lig_start+1}-{lig_end+1} and name CA" # ligand

lenalidomide = u_ref.select_atoms('resname MOL') 
ligand = u_ref.select_atoms(ligand_selection_str)
receptor = u_ref.select_atoms(receptor_selection_str)

# Find indices of lenalidomide heavy atoms
lenalidomide_indices = []
for atom in lenalidomide.atoms:
    if atom.name[0] != ('H'):
        lenalidomide_indices.append(atom.index)

ligand_interface = []
receptor_interface = []

# Save all ligand CA indices that are within 12 AA of a receptor CA
for lig_CA in ligand.atoms:
    
    for rec_CA in receptor.atoms:
        distance = getDistance(lig_CA.index, rec_CA.index, u_ref)
        if distance <=12.000: 
            ligand_interface.append(int(lig_CA.index))
            receptor_interface.append(int(rec_CA.index))

# Returns interface residues
ligand_interface=set(ligand_interface)
receptor_interface=set(receptor_interface)   

rec_interface_res_init = [
    int(u_ref.select_atoms(f"index {index}").residues.resids[0]) for index in receptor_interface 
]

lig_interface_res_init = [
    int(u_ref.select_atoms(f"index {index}").residues.resids[0]) for index in ligand_interface 
]

"""Calculating the interface RMSD"""
print('Calculating the interface RMSD...\n')

if not os.path.exists('results/av_iRMSD.csv'):
    time, av_RMSD = obtain_av_interface_RMSD([1,2,3], rec_interface_res_init, lig_interface_res_init)
    df = pd.DataFrame()
    df['Time (ns)'] = time
    df['Interface RMSD (Angstrom)'] = av_RMSD
    df.to_csv('results/av_iRMSD.csv')

"""Identifying interface equilibration"""
print('Calculating interface equilibration time')

df = pd.read_csv('results/av_iRMSD.csv')
time_interface, iRMSD = df['Time (ns)'].to_numpy(), df['Interface RMSD (Angstrom)'].to_numpy()

equil_idx, g, ess = red.detect_equilibration_window(iRMSD,
                                              method="min_sse",
                                              plot=False)

# If equilibration not detected in first 60 ns, apply an equilibration time of 60 ns
u1 = mda.Universe("structures/complex.prmtop", "results/run1/traj.dcd")
n_frames = len(u1.trajectory)

if not isinstance(equil_idx, np.int64) or equil_idx > int(0.6 * n_frames): # RED unable to assign idx
    equil_idx = int(0.6 * n_frames)

interface_equilibration_time = time_interface[equil_idx]
config_data['Interface equilibration time'] = float(np.round(interface_equilibration_time, 4))
print(f"Equilibration time detected by RED is {interface_equilibration_time} ns")

# Generate RMSD plot
run_numbers = [1,2,3]

# Define systems for plotting
systems = ['receptor', 'ligand']

start_idx = 0 # cutoff index for RMSD plotting
stop_idx = -1

# Create subplots for each RMSD plot
fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(6, 8), dpi=500, sharex=True)

# Plot Interface RMSD
axs[0].plot(time_interface[start_idx:stop_idx], iRMSD[start_idx:stop_idx], c='grey', label='Interface')
axs[0].vlines(x=time_interface[equil_idx], ymin=0, ymax=3, colors='k', linestyle='dotted')
axs[0].set_ylabel('RMSD (Å)')
axs[0].set_title('Interface RMSD')

colors = {'ligand':'g', 'receptor':'b'}

# Plot each system's RMSD
for i, species in enumerate(systems, start=1):
    time, RMSD = obtain_RMSD(run_numbers, species)
    axs[i].plot(time[start_idx:stop_idx], RMSD[start_idx:stop_idx], c=colors[species], label=species)
    axs[i].set_ylabel('RMSD (Å)')
    axs[i].set_title(f'{species} RMSD')

plt.xlabel('Time (ns)')
plt.tight_layout()
plt.savefig('plots/RMSD.png')

"""Calculating equilibrated interface residues"""
print('Calculating the interface residues after equilibration...\n')

u = mda.Universe('structures/complex.prmtop', f'results/run1/traj.dcd')

u.trajectory[equil_idx]

receptor_selection_str=f"resid {rec_start+1}-{rec_end+1} and name CA" # receptor
ligand_selection_str = f"resid {lig_start+1}-{lig_end+1} and name CA" # ligand

lenalidomide = u.select_atoms('resname MOL') 
ligand = u.select_atoms(ligand_selection_str)
receptor = u.select_atoms(receptor_selection_str)

# Find indices of lenalidomide heavy atoms
lenalidomide_indices = []
for atom in lenalidomide.atoms:
    if atom.name[0] != ('H'):
        lenalidomide_indices.append(atom.index)

ligand_interface = []
receptor_interface = []

# Save all ligand CA indices that are within 12 AA of a receptor CA
for lig_CA in ligand.atoms:
    
    for rec_CA in receptor.atoms:
        distance = getDistance(lig_CA.index, rec_CA.index, u)
        if distance <=12.000: 
            ligand_interface.append(int(lig_CA.index))
            receptor_interface.append(int(rec_CA.index))

# Returns interface residues
ligand_interface=set(ligand_interface)
lig_group = sorted(list(ligand_interface)) 

receptor_interface=set(receptor_interface)   
rec_group = sorted(list(receptor_interface)) + lenalidomide_indices 

receptor_interface_res = [
    int(u.select_atoms(f"index {index}").residues.resids[0]) for index in receptor_interface 
]

ligand_interface_res = [
    int(u.select_atoms(f"index {index}").residues.resids[0]) for index in ligand_interface 
]

# Add lenalidomide resid to the list of indices
# receptor_interface_res.append(lenalidomide.resids[0])

print(f"\nReceptor interface residues = {sorted(receptor_interface_res)}")
print(f"Ligand interface residues = {sorted(ligand_interface_res)}\n")

config_data['Receptor interface residues'] = sorted(receptor_interface_res)
config_data['Ligand interface residues'] = sorted(ligand_interface_res)

"""Check that there is no anchor point collinearity and 
save the equilibrated structure"""

u = mda.Universe('structures/complex.prmtop', f'results/run1/traj.dcd')

# Save structure after equil_time + 5 ns, RED tends to identify short equilibration times
save_idx = equil_idx + int(5/100.2 * len(u.trajectory))
u.trajectory[save_idx]

# Check for collinearity
pos_a = u.atoms[rec_group].center_of_mass()
pos_b = u.atoms[[idx_b]].center_of_mass()
pos_c = u.atoms[[idx_c]].center_of_mass()
pos_A = u.atoms[lig_group].center_of_mass()
pos_B = u.atoms[[idx_B]].center_of_mass()
pos_C = u.atoms[[idx_C]].center_of_mass()

tolerance = np.radians(10) # Raise error for angle between 170 and 190 deg 

angle = obtain_angle(pos_c, pos_b, pos_a)
if np.pi - tolerance < angle < np.pi + tolerance:
    raise ValueError(f"Possible collinearity detected between anchor points c, b and a...")

angle = obtain_angle(pos_b, pos_a, pos_A)
if np.pi - tolerance < angle < np.pi + tolerance:
    raise ValueError(f"Possible collinearity detected between anchor points b, a and A...")

angle = obtain_angle(pos_a, pos_A, pos_B)
if np.pi - tolerance < angle < np.pi + tolerance:
    raise ValueError(f"Possible collinearity detected between anchor points a, A and B...")

angle = obtain_angle(pos_A, pos_B, pos_C)
if np.pi - tolerance < angle < np.pi + tolerance:
    raise ValueError(f"Possible collinearity detected between anchor points A, B and C...")

# Save the equilibrated frame to Unrestrained_MD and umbrella_sampling directories
u.atoms.write('equilibrated_structures/eq_frame.pdb')
u.atoms.write('../umbrella_sampling/equilibrated_structures/eq_frame.pdb')

print(f"\nSave equilibrated snapshot from run 1 after {equil_time+5.0} ns\n")

"""Calculating Boresch distributions"""
print('Calculating Boresch histograms...\n')

# Calculate and save Boresch DOF samples
for run_number in [1,2,3]:

    for dof in ['thetaA', 'thetaB', 'phiA', 'phiB', 'phiC']:

        if os.path.exists(f'results/run{run_number}/{dof}.pkl'):
            continue
        else:
            print(f"Performing Boresch analysis for {dof} run {run_number}")

            u = mda.Universe('structures/complex.prmtop', f'results/run{run_number}/traj.dcd')
            transformation = NoJump()
            u.trajectory.add_transformations(transformation)

            vals = []

            for ts in tqdm(u.trajectory, total=(u.trajectory.n_frames), desc='Frames analysed'):
                vals.append(obtain_Boresch_dof(dof, u, rec_group, lig_group, res_b, res_c, res_B, res_C))

            frames = np.arange(1, len(vals) + 1)

            dof_data = {
                'Frames': frames,
                'Time (ns)': np.round(0.01 * frames, 6),
                'DOF values': vals
            }

            # Save interface data to pickle
            file = f'results/run{run_number}/{dof}.pkl'
            with open(file, 'wb') as f:
                pickle.dump(dof_data, f)

# Functions to read in data for plotting
def obtain_Boresch_dof(run_number, dof):

    boreschfile = f'results/run{run_number}/{dof}.pkl'

    with open(boreschfile, 'rb') as f:
        loaded_data = pickle.load(f)

    frames = loaded_data['Frames']
    time = loaded_data['Time (ns)']
    vals = loaded_data['DOF values']

    return time, vals

def correct_torsions(vals):
    """
    Correct an array/list of dihedral angles so that the domain 
    -pi -> pi is adhered to
    """

    if (max(vals)-min(vals))>6.0:
        if np.average(vals) > 0:
            for n in range(len(vals)):
                if vals[n] < 0:
                    vals[n] = vals[n] + 2*np.pi
        else:
            for n in range(len(vals)):
                if vals[n] > 0:
                    vals[n] = vals[n] - 2*np.pi    

    return vals 

# Titles for plotting
dof_titles = {
    'thetaA': r'$\Theta _A$', 
    'thetaB': r'$\Theta _B$', 
    'phiA': r'$\phi _A$', 
    'phiB': r'$\phi _B$', 
    'phiC': r'$\phi _C$'
}

# Plot time evolution of each DOF
fig, axes = plt.subplots(nrows=5, ncols=1, figsize=(3, 12), dpi=200)

for idx, dof in enumerate(dof_titles.keys()):

    for run_number in [1,2,3]:

        time, vals = obtain_Boresch_dof(run_number, dof)
        vals = correct_torsions(vals)

        axes[idx].plot(time, vals, label=f"run {run_number}")
        axes[idx].set_title(dof_titles[dof])
        axes[idx].set_ylabel('radians')
        axes[idx].set_xlabel('time (ns)')
        axes[idx].legend()

plt.tight_layout()
plt.savefig("plots/Boresch_sampling.png")

# Plot the histograms
eq_values = {}

equil_time = interface_equilibration_time # ns
cutoff_time = 100.0  # ns

dof_titles = {
    'thetaA': r'$\Theta _A$', 
    'thetaB': r'$\Theta _B$', 
    'phiA': r'$\phi _A$', 
    'phiB': r'$\phi _B$', 
    'phiC': r'$\phi _C$'
}

fig, axes = plt.subplots(nrows=5, ncols=1, figsize=(3, 12), dpi=500)

for idx, dof in enumerate(dof_titles.keys()):

    vals_all = []

    for run_number in [1, 2, 3]:
        time, vals = obtain_Boresch_dof(run_number, dof)

        start_idx = np.searchsorted(time, equil_time, side='left')
        end_idx = np.searchsorted(time, cutoff_time, side='right')

        time_slice = time[start_idx:end_idx]
        vals_slice = vals[start_idx:end_idx]

        vals_all.append(vals_slice)

    if len(vals_all) == 0:
        raise RuntimeError("No data in the requested time window.")

    vals_all = np.concatenate(vals_all)
    vals_all = correct_torsions(vals_all)

    counts, bins = np.histogram(vals_all, bins=50)
    bin_centers = 0.5 * (bins[:-1] + bins[1:])
    peak_center = bin_centers[counts.argmax()]

    axes[idx].hist(vals_all, bins=50)
    axes[idx].vlines(peak_center, ymin=0, ymax=counts.max() * 1.1, colors='k', linestyles='dotted', label='Equilibrium value')
    axes[idx].set_title(dof_titles[dof])
    axes[idx].set_yticks([])
    axes[idx].set_xlabel('radians')

    eq_values[f"{dof[:-1]}_{dof[-1]}_0"] = float(peak_center)

plt.tight_layout()
plt.savefig("plots/Boresch_hist.png")

config_data['Boresch equilibrium values'] = eq_values

# Check all the relevant info has been saved
keys = [
    'Interface equilibration time',
    'Boresch anchor points',
    'Boresch equilibrium values',
    'Receptor residues',
    'Ligand residues',
    'Receptor interface residues',
    'Ligand interface residues'
]

for key in keys:
    if config_data.get(key) is None:
        raise ValueError(f"Key '{key}' missing or has no value in yaml...")

print('\nFinal dictionary to save:\n')

print(config_data)

# Save final dictionary to yaml
with open('../umbrella_sampling/US_config.yaml', 'w') as file:
    yaml.safe_dump(config_data, file)


