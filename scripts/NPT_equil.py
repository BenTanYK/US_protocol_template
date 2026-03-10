"""This is a short MD simulation script used for performing a short NPT 
simulation for a set of Amber input files $FILENAME.prmtop/$FILENAME.inpcrd.

The script takes the name of the input files as the first argument,
e.g. 'python NPT_equil.py input' -> this yields the output input_NPT_equil.inpcrd"""

import openmm as mm
import openmm.app as app
import openmm.unit as unit
import sys

system = str(sys.argv[1])

"""System setup"""

dt = 2*unit.femtoseconds 

# Load param and coord files
prmtop = app.AmberPrmtopFile(f'{system}.prmtop')
inpcrd = app.AmberInpcrdFile(f'{system}.inpcrd')

system = prmtop.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=1.0*unit.nanometer, constraints=app.HBonds)  
integrator = mm.LangevinMiddleIntegrator(6.0000*unit.kelvin, 1.0000/unit.picosecond, dt)

simulation = app.Simulation(prmtop.topology, system, integrator)
simulation.context.setPositions(inpcrd.positions)

# Add reporter
simulation.reporters.append(app.StateDataReporter(sys.stdout, 1000, step=True, time=True, potentialEnergy=True, temperature=True, speed=True))

# Minimise energy 
simulation.minimizeEnergy()
simulation.context.setVelocitiesToTemperature(6.0000*unit.kelvin)

"""System heating"""

for i in range(50):
    integrator.setTemperature(6*(i+1)*unit.kelvin)
    simulation.step(1000)

"""NPT simulation"""

system.addForce(mm.MonteCarloBarostat(1.0*unit.atmospheres, 300*unit.kelvin, 25))
simulation.context.reinitialize(preserveState=True)

# 200 ps simulation
simulation.step(100000)

"""Save NPT-equilibrated coordinates"""

state = simulation.context.getState(getPositions=True, enforcePeriodicBox=True)
positions = state.getPositions(asNumpy=True)   
a, b, c = state.getPeriodicBoxVectors()        

structure = pmd.openmm.load_topology(
    simulation.topology,   # or prmtop.topology
    system=system,
    xyz=positions
)

structure.box_vectors = (a, b, c)

# Save NPT-equilibrated coordinates
structure.save(f"{system}_NPT_eq.inpcrd", format="rst7", overwrite=True)