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

binding_file = '{jobname}_@n@.out'
f = open(binding_file, 'w')

temperature = @temperature@ * kelvin
lmbd = @lambda@
lambda1 = @lambda1@
lambda2 = @lambda2@
alpha = @alpha@ / kilocalorie_per_mole
u0 = @u0@ * kilocalorie_per_mole
w0coeff = @w0coeff@ * kilocalorie_per_mole

umsc = {soft_core_umax} * kilocalorie_per_mole
acore = {soft_core_acore}

print("temperature = ", temperature)
print("lambda = ", lmbd)
print("lambda1 = ", lambda1)
print("lambda2 = ", lambda2)
print("alpha = ", alpha)
print("u0 = ", u0)
print("w0coeff = ", w0coeff)
print("soft core method = ", '{soft_core_method}')
print("umax = ", umsc)
print("acore = ", acore)

rcptfile_input  = '{jobname}_rcpt_@nm1@.dms'
ligfile_input   = '{jobname}_lig_@nm1@.dms'
rcptfile_output = '{jobname}_rcpt_tmp.dms'
ligfile_output  = '{jobname}_lig_tmp.dms'
rcptfile_result = '{jobname}_rcpt_@n@.dms'
ligfile_result  = '{jobname}_lig_@n@.dms'

shutil.copyfile(rcptfile_input, rcptfile_output)
shutil.copyfile(ligfile_input, ligfile_output)

testDes = DesmondDMSFile([ligfile_output, rcptfile_output]) 
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

platform_name = '@platform@'
platform = Platform.getPlatformByName(platform_name)

properties = {{}}

if platform_name =='OpenCL':
    #expected "platformid:deviceid" or empty
    device = "@pn@"
    m = re.match("(\d+):(\d+)", device)
    if m:
        platformid = m.group(1)
        deviceid = m.group(2)
        properties["OpenCLPlatformIndex"] = platformid
        properties["DeviceIndex"] = deviceid
        print("Using platform id: %s, device id: %s" % ( platformid ,  deviceid) )

frictionCoeff = 0.5 / picosecond
MDstepsize = 0.001 * picosecond

integrator = LangevinIntegratorSDM(temperature/kelvin, frictionCoeff/(1/picosecond), MDstepsize/ picosecond, lig_atoms)
integrator.setBiasMethod(sdm_utils.ILogisticMethod)
integrator.setLambda(lmbd)
integrator.setLambda1(lambda1)
integrator.setLambda2(lambda2)
integrator.setAlpha(alpha*kilojoule_per_mole)
integrator.setU0(u0/ kilojoule_per_mole)
integrator.setW0coeff(w0coeff / kilojoule_per_mole)    

soft_core_method = sdm_utils.{soft_core_method}
integrator.setSoftCoreMethod(soft_core_method)
integrator.setUmax(umsc / kilojoule_per_mole)
integrator.setAcore(acore)

simulation = Simulation(testDes.topology, system, integrator,platform, properties)
simulation.context.setPositions(testDes.positions)
simulation.context.setVelocities(testDes.velocities)

totalSteps = {nsteps}
nprnt = {nprnt}
ntrj = {ntrj}
simulation.reporters.append(StateDataReporter(stdout, nprnt, step=True, temperature=True))
simulation.reporters.append(DCDReporter("{jobname}_@n@.dcd", ntrj))
simulation.reporters.append(PDBReporter("{jobname}_@n@.pdb", totalSteps))

loops = totalSteps/nprnt
start=datetime.now()
step = 0
for i in range(loops):
    simulation.step(nprnt)
    pot_energy = (integrator.getPotEnergy()*kilojoule_per_mole).value_in_unit(kilocalorie_per_mole)
    bind_energy = (integrator.getBindE()*kilojoule_per_mole).value_in_unit(kilocalorie_per_mole)
    print("%f %f %f %f %f %f %f %f %f" % (temperature/kelvin,lmbd, lambda1, lambda2, alpha*kilocalorie_per_mole, u0/kilocalorie_per_mole, w0coeff/kilocalorie_per_mole, pot_energy, bind_energy), file=f )
    f.flush()
    step += nprnt
end=datetime.now()

positions = simulation.context.getState(getPositions=True).getPositions()
velocities = simulation.context.getState(getVelocities=True).getVelocities()
testDes.setPositions(positions)
testDes.setVelocities(velocities)
testDes.close()

f.close()

shutil.copyfile(rcptfile_output, rcptfile_result)
shutil.copyfile(ligfile_output, ligfile_result)

elapsed=end - start
print("MD time="+str(elapsed.seconds+elapsed.microseconds*1e-6)+"s")
