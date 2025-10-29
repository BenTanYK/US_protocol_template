"""This script calculates the RMSD and RMSF of the unrestrained MD trajectories.
We assume there are three runs - change this if necessary"""

import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis import rms, align
import yaml

# Read in the initial residue selection
with open('config.yaml', 'r') as infile:
    initial_config = yaml.safe_load(infile)

# 0-indexing
rec_start = int(initial_config['Receptor residue range'][0])
rec_end = int(initial_config['Receptor residue range'][-1])
lig_start = int(initial_config['Ligand residue range'][0])
lig_end = int(initial_config['Ligand residue range'][-1])

# Convert to 1-indexing for MDAnalysis selections
rec_start += 1
rec_end += 1
lig_start += 1
lig_end += 1

def obtain_RMSD(run_number, res_range=[rec_start, lig_end]):
    u = mda.Universe('structures/complex.prmtop', f'results/run{run_number}/traj.dcd')
    protein = u.select_atoms("protein")
    ref = protein

    R_u =rms.RMSD(protein, ref, select=f'backbone and resid {res_range[0]}-{res_range[1]}')
    R_u.run()

    rmsd_u = R_u.rmsd.T #take transpose
    time = rmsd_u[1]/1000
    rmsd= rmsd_u[2]

    return time, rmsd

def save_RMSD(run_number, res_range=[rec_start, lig_end], filename='RMSD.csv'):
    """
    Save the RMSD of a given run in a .csv file
    """
    time, RMSD = obtain_RMSD(run_number, res_range)

    df = pd.DataFrame()
    df['Time (ns)'] = time
    df['RMSD (Angstrom)'] = RMSD

    df.to_csv(f"results/run{run_number}/{filename}")

    return df

def obtain_RMSF(run_number, res_range=[rec_start, lig_end]):
    u = mda.Universe('structures/complex.prmtop', f'results/run{run_number}/traj.dcd')
    
    start, end = int(res_range[0]), int(res_range[1])

    alignment_selection = f'protein and name CA and resid {start}-{end}'
    c_alphas = u.select_atoms(alignment_selection)
    if len(c_alphas) == 0:
        raise ValueError(f"No atoms selected with selection: '{alignment_selection}'")

    # build average structure 
    avg = align.AverageStructure(u, select=alignment_selection, ref_frame=0)
    avg.run()
    ref = avg.results.universe

    # align trajectory in memory 
    align.AlignTraj(u, ref, select=alignment_selection, in_memory=True).run()

    # compute RMSF
    R = rms.RMSF(c_alphas)
    R.run()

    return c_alphas.resids, R.results.rmsf

def save_RMSF(run_number, res_range=[rec_start, lig_end], filename='RMSF.csv'):
    """
    Save the RMSF of a given run in a .csv file
    """
    residx, RMSF = obtain_RMSF(run_number, res_range)

    df = pd.DataFrame()
    df['Residue index'] = residx
    df['RMSF (Angstrom)'] = RMSF
    df.to_csv(f"results/run{run_number}/{filename}")

    return df

for n_run in [1,2,3]:
    # complex
    print(f"\nCalculating complex RMSD for run {n_run}")
    save_RMSD(n_run)

    # receptor
    print(f"\nCalculating complex RMSD for run {n_run}")
    save_RMSD(n_run, [rec_start, rec_end], 'RMSD_rec.csv')

    # ligand
    print(f"\nCalculating complex RMSD for run {n_run}")
    save_RMSD(n_run, [lig_start, lig_end], 'RMSD_lig.csv')

    # complex
    print(f"\nGenerating RMSF for  run {n_run}")
    save_RMSF(n_run)   

    # receptor
    print(f"\nCalculating complex RMSF for run {n_run}")
    save_RMSF(n_run, [rec_start, rec_end], 'RMSF_rec.csv')

    # ligand
    print(f"\nCalculating complex RMSF for run {n_run}")
    save_RMSF(n_run, [lig_start, lig_end], 'RMSF_lig.csv')

