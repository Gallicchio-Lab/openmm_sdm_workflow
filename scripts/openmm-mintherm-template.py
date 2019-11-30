from __future__ import print_function

from simtk import openmm as mm
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import os, re,time, shutil, math
from simtk.openmm.app.desmonddmsfile import *
from datetime import datetime
from SDMplugin import *

print("Started at: " + str(time.asctime()))
start=datetime.now()

shutil.copyfile('{jobname}_rcpt.dms','{jobname}_rcpt_0.dms') 
shutil.copyfile('{jobname}_lig.dms','{jobname}_lig_0.dms') 

testDes = DesmondDMSFile(['{jobname}_lig_0.dms','{jobname}_rcpt_0.dms']) 
#not using a cutoff because AGBNP does not have a smooth switching function
#stability of binding energy is poor
#system = testDes.createSystem(nonbondedMethod=CutoffNonPeriodic,nonbondedCutoff=12.0*nanometer, OPLS = True, implicitSolvent= '{implicitsolvent}')
system = testDes.createSystem(nonbondedMethod=NoCutoff, OPLS = True, implicitSolvent= '{implicitsolvent}')


natoms_ligand = {nlig_atoms}
lig_atoms = range(natoms_ligand)
# atom indexes here refer to indexes in either lig or rcpt dms file, rather than in the complex 
lig_atom_restr = {rest_ligand_cmlig}   #indexes of ligand atoms for CM-CM Vsite restraint
rcpt_atom_restr = {rest_ligand_cmrec}   #indexes of rcpt atoms for CM-CM Vsite restraint
kf = {rest_ligand_cmkf} * kilocalorie_per_mole/angstrom**2 #force constant for Vsite CM-CM restraint 
r0 = {rest_ligand_cmtol} * angstrom #radius of Vsite sphere

#these can be 'None" if not using orientational restraints
lig_ref_atoms = {ligand_ref_atoms} # the 3 atoms of the ligand that define the coordinate system of the ligand
rcpt_ref_atoms = {receptor_ref_atoms} # the 3 atoms of the receptor that define the coordinate system of the receptor
angle_center = {angle_center} * degrees
kfangle = {kfangle} * kilocalorie_per_mole/degrees**2
angletol = {angletol} * degrees
dihedral1center = {dihedral1center} * degrees
kfdihedral1 = {kfdihedral1} * kilocalorie_per_mole/degrees**2
dihedral1tol = {dihedral1tol} * degrees
dihedral2center = {dihedral2center} * degrees
kfdihedral2 = {kfdihedral2} * kilocalorie_per_mole/degrees**2
dihedral2tol = {dihedral2tol} * degrees

#transform indexes of receptor atoms
for i in range(len(rcpt_atom_restr)):
    rcpt_atom_restr[i] += natoms_ligand
if rcpt_ref_atoms:
    for i in range(len(rcpt_ref_atoms)):
        rcpt_ref_atoms[i] += natoms_ligand

sdm_utils = SDMUtils(system, lig_atoms)
sdm_utils.addRestraintForce(lig_cm_particles = lig_atom_restr,
                            rcpt_cm_particles = rcpt_atom_restr,
                            kfcm = kf,
                            tolcm = r0,
                            lig_ref_particles = lig_ref_atoms,
                            rcpt_ref_particles = rcpt_ref_atoms,
                            angle_center = angle_center,
                            kfangle = kfangle,
                            angletol = angletol,
                            dihedral1center = dihedral1center,
                            kfdihedral1 = kfdihedral1,
                            dihedral1tol = dihedral1tol,
                            dihedral2center = dihedral2center,
                            kfdihedral2 = kfdihedral2,
                            dihedral2tol = dihedral2tol)

initial_temperature = 50 * kelvin
final_temperature = {temperature} * kelvin

temperature = initial_temperature
frictionCoeff = 0.5 / picosecond
MDstepsize = 0.001 * picosecond
    
integrator = LangevinIntegrator(temperature/kelvin, frictionCoeff/(1/picosecond), MDstepsize/ picosecond)

platform_name = '{platform}'
platform = Platform.getPlatformByName(platform_name)

properties = {{}}

simulation = Simulation(testDes.topology, system, integrator,platform, properties)
print ("Using platform %s" % simulation.context.getPlatform().getName())
simulation.context.setPositions(testDes.positions)
simulation.context.setVelocities(testDes.velocities)
context=simulation.context

state = simulation.context.getState(getEnergy = True)
print(state.getPotentialEnergy())
print("Energy minimizing the system ...")
simulation.minimizeEnergy()
state = simulation.context.getState(getEnergy = True)
print(state.getPotentialEnergy())

stepId = 5000
totalSteps = 50000
loopStep = totalSteps/stepId
delta_temperature = (final_temperature - initial_temperature)/loopStep
simulation.reporters.append(StateDataReporter(stdout, stepId, step=True, potentialEnergy = True, temperature=True))

#MD with temperature ramp
for i in range(loopStep):
    simulation.step(stepId)
    #prepare system for new temperature
    temperature = temperature + delta_temperature
    integrator.setTemperature(temperature)

positions = simulation.context.getState(getPositions=True).getPositions()
velocities = simulation.context.getState(getVelocities=True).getVelocities()
print( "Updating positions and velocities")
testDes.setPositions(positions)
testDes.setVelocities(velocities)
testDes.close()

#save a pdb file that can be used as a topology to load .dcd files in vmd
with open('{jobname}.pdb', 'w') as output:
  PDBFile.writeFile(simulation.topology, positions, output)

end=datetime.now()
elapsed=end - start
print("elapsed time="+str(elapsed.seconds+elapsed.microseconds*1e-6)+"s")
