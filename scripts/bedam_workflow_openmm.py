from __future__ import print_function

# bedam python class
__doc__="""
$Revision: 0.1 $

A class to prepare OpenMM AsyncRE jobs

"""
# Contributors: Emilio Gallicchio

import os, sys, time, re, glob
from schrodinger.utils import cmdline
import schrodinger.utils.log
import shutil
import signal
import glob
import time
import sqlite3 as sqlite
import StringIO

from bedam_asyncre import bedam_job_asyncre

class bedam_job_openmm_asyncre(bedam_job_asyncre):

    #asyncre override control file for openmm
    def writeCntlFile(self):
        input = ""

        job_transport = self.keywords.get('JOB_TRANSPORT')
        if job_transport is None:
            msg = "writeCntlFile: JOB_TRANSPORT is not specified"
            self.exit(msg)
        if not (job_transport == "SSH"):
            msg = "writeCntlFile: invalid JOB_TRANSPORT: %s Must be 'SSH'." % job_transport
            self.exit(msg)
        input += "JOB_TRANSPORT = '%s'\n" % job_transport
        
        re_type = self.keywords.get('RE_TYPE')
        if re_type is None:
            msg = "writeCntlFile: RE_TYPE is not specified"
            self.exit(msg)
        if not (re_type == 'BEDAMTEMPT'):
            msg = "writeCntlFile: invalid RE_TYPE. Must be 'BEDAMTEMPT'."
            self.exit(msg)
        input += "RE_TYPE = '%s'\n" % re_type

        engine = self.keywords.get('ENGINE')
        if engine is None:
            engine = "OPENMM"
        input += "ENGINE = '%s'\n" % engine

        input += "ENGINE_INPUT_BASENAME = '%s'\n" % self.jobname
		
        input += "RE_SETUP = 'YES'\n"
        
        extfiles = self.keywords.get('ENGINE_INPUT_EXTFILES')
        required_files = "runopenmm"
        if extfiles is None:
            extfiles = required_files
        else:
            extfiles += ",%s" % required_files
        input_file = "%s.py" % self.jobname

        if re_type == 'BEDAMTEMPT':
            rcptfile =  self.jobname + '_rcpt_0' + '.dms'
            ligfile =  self.jobname + '_lig_0' + '.dms'
            if job_transport == 'SSH':
		extfiles += ",%s,%s" % (rcptfile,ligfile)
            input += "ENGINE_INPUT_EXTFILES = '%s'\n" % extfiles

        temperatures = self.keywords.get('TEMPERATURES')
        if temperatures is not None:
            input += "TEMPERATURES = '%s'\n" % temperatures
        
        lambdas = self.keywords.get('LAMBDAS')
        if lambdas is not None:
            input += "LAMBDAS = '%s'\n" % lambdas
        else:
            msg = "writeCntlFile: 'LAMBDAS' is required."
            self.exit(msg)
            
        nlambdas = len(lambdas.split(","))
        zerosdefault = "0.000"
        for i in range(1,nlambdas):
            zerosdefault += ",0.000"

        lambda1 = self.keywords.get('LAMBDA1')
        if lambda1 is not None:
            input += "LAMBDA1 = '%s'\n" % lambda1
        else:
            input += "LAMBDA1 = '%s'\n" % lambdas     
        
        lambda2 = self.keywords.get('LAMBDA2')
        if lambda2 is not None:
            input += "LAMBDA2 = '%s'\n" % lambda2
        else:
            input += "LAMBDA2 = '%s'\n" % lambdas

        alpha = self.keywords.get('ALPHA')
        if alpha is not None:
            input += "ALPHA = '%s'\n" % alpha
        else:
            input += "ALPHA = '%s'\n" % zerosdefault
        
        u0 = self.keywords.get('U0')
        if u0 is not None:
            input += "U0 = '%s'\n" % u0
        else:
            input += "U0 = '%s'\n" % zerosdefault
            
        w0coeff = self.keywords.get('W0COEFF')
        if w0coeff is not None:
            input += "W0COEFF = '%s'\n" % w0coeff
        else:
            input += "W0COEFF = '%s'\n" % zerosdefault
            
        wall_time = self.keywords.get('WALL_TIME')
        if wall_time is not None:
            input += "WALL_TIME = %d\n" % int(wall_time)

        replica_run_time = self.keywords.get('REPLICA_RUN_TIME')
        if replica_run_time is not None:
            input += "REPLICA_RUN_TIME = %d\n" % int(replica_run_time)

        cycle_time = self.keywords.get('CYCLE_TIME')
        if cycle_time is not None:
            input += "CYCLE_TIME = %d\n" % int(cycle_time)
            
        if job_transport == 'SSH':
            input += "NODEFILE = 'nodefile'\n"
            total_cores = self.keywords.get('TOTAL_CORES')
            if total_cores is None:
                msg = "writeCntlFile: TOTAL_CORES is required"
                self.exit(msg)
            input += "TOTAL_CORES = %d\n" % int(total_cores) 
            subjob_cores = self.keywords.get('SUBJOB_CORES')
            if subjob_cores is not None:
                input += "SUBJOB_CORES = %d\n" % int(subjob_cores)

        subjobs_buffer_size = self.keywords.get('SUBJOBS_BUFFER_SIZE')
        if subjobs_buffer_size is not None:
            input += "SUBJOBS_BUFFER_SIZE = '%f'\n" % float(subjobs_buffer_size)

        nsteps = self.keywords.get('PRODUCTION_STEPS')
        nprnt = self.keywords.get('PRNT_FREQUENCY')
        if (not nsteps) or (not nprnt):
            msg = "writeCntlFile: PRODUCTION_STEPS and PRNT_FREQUENCY are required"
            self.exit(msg)
        input += "PRODUCTION_STEPS = '%d'\n" % int(nsteps)
        input += "PRNT_FREQUENCY = '%d'\n" % int(nprnt)

        input += "IMPLICITSOLVENT = 'AGBNP'\n"
            
        verbose = self.keywords.get('VERBOSE')
        if verbose is not None:
            input += "VERBOSE = '%s'\n" % verbose
        
        cntlfile = "%s_asyncre.cntl" % self.jobname
        f = open(cntlfile, "w")
        f.write(input)
        f.close

#asyncre override control file for openmm
    def writeOpenMMDriverFile(self):
        input_openmm = """
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

lmbd = @lambda@
lambda1 = @lambda1@
lambda2 = @lambda2@
alpha = @alpha@ / kilocalorie_per_mole
u0 = @u0@ * kilocalorie_per_mole
w0coeff = @w0coeff@ * kilocalorie_per_mole

umsc = {soft_core_umax} * kilocalorie_per_mole
acore = {soft_core_acore}

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

temperature = @temperature@ * kelvin
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
    f.write("%f %f %f %f %f %f %f %f\\n" % (lmbd, lambda1, lambda2, alpha*kilocalorie_per_mole, u0/kilocalorie_per_mole, w0coeff/kilocalorie_per_mole, pot_energy, bind_energy))
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
"""

        #figure out number of atoms and indexes of Vsite restrained atoms
        rcptfile = self.jobname + '_rcpt' + '.dms'
        conn = sqlite.connect(rcptfile)
        q = """SELECT id FROM particle"""
        n_rcpt = 0;
        for atom_id in conn.execute(q):
            n_rcpt += 1
        receptor_sql =  self.keywords.get('REST_LIGAND_CMRECSQL')
        q = "SELECT id FROM particle WHERE %s" % receptor_sql
        rcpt_ids = []
        for atom_id in conn.execute(q):
            rcpt_ids.append( int(atom_id[0]) )
        rcpt_atom_id = str(rcpt_ids)    
        conn.close()
        
        
        ligfile =  self.jobname + '_lig' + '.dms'
        conn = sqlite.connect(ligfile)
        q = """SELECT id FROM particle"""
        n_lig = 0;
        for atom_id in conn.execute(q):
            n_lig += 1
        ligand_sql =  self.keywords.get('REST_LIGAND_CMLIGSQL')
        q = "SELECT id FROM particle WHERE %s" % ligand_sql
        lig_ids = []
        for atom_id in conn.execute(q):
            lig_ids.append( int(atom_id[0]) )
        lig_atom_id = str(lig_ids)
        conn.close()
        
        #force constant etc.
        kf = float(self.keywords.get('REST_LIGAND_CMKF'))
        tol = float(self.keywords.get('REST_LIGAND_CMTOL'))

        #number of MD steps etc
        nsteps = int(self.keywords.get('PRODUCTION_STEPS'))
        nprnt =  int(self.keywords.get('PRNT_FREQUENCY'))
        ntrj =  int(self.keywords.get('TRJ_FREQUENCY'))

        #implicit solvent
        if self.keywords.get('IMPLICIT_SOLVENT') is None:
            implicitsolvent = 'AGBNP'
        else:
            implicitsolvent = self.keywords.get('IMPLICIT_SOLVENT')
        

        #soft core settings
        soft_core_method = self.keywords.get('SOFT_CORE_METHOD')
        if soft_core_method is None:
            soft_core_method = 'NoSoftCoreMethod'
        soft_core_umax = self.keywords.get('SOFT_CORE_UMAX')
        if soft_core_umax is None:
            soft_core_umax = 50.0 #kcal/mol
        else:
            soft_core_umax = float(soft_core_umax)
        soft_core_acore = self.keywords.get('SOFT_CORE_ACORE')
        if soft_core_acore is None:
            soft_core_acore = 1.0
        else:
            soft_core_acore = float(soft_core_acore)

        inputr = input_openmm.format(
            jobname = self.jobname,
            nlig_atoms = n_lig,
            rest_ligand_cmrec = rcpt_atom_id,
            rest_ligand_cmlig = lig_atom_id,
            rest_ligand_cmkf = kf,
            rest_ligand_cmtol = tol,
            ligand_ref_atoms = None,
            receptor_ref_atoms = None,
            angle_center = None,
            angletol = None,
            kfangle = None,
            dihedral1center = None,
            kfdihedral1 = None,
            dihedral1tol = None,
            dihedral2center = None,
            kfdihedral2 = None,
            dihedral2tol = None,
            implicitsolvent = implicitsolvent,
            soft_core_method = soft_core_method,
            soft_core_umax = soft_core_umax,
            soft_core_acore = soft_core_acore,
            nsteps = nsteps,
            nprnt = nprnt,
            ntrj = ntrj
        )
        
        driverfile = "%s.py" % self.jobname
        f = open(driverfile, "w")
        f.write(inputr)
        f.close()


#asyncre override control file for openmm
    def writeThermInputFile(self):
        mintherm_openmm = """
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

end=datetime.now()
elapsed=end - start
print("elapsed time="+str(elapsed.seconds+elapsed.microseconds*1e-6)+"s")
"""

        temperature = self.keywords.get('TEMPERATURE')

        #figure out number of atoms and indexes of Vsite restrained atoms
        rcptfile = self.jobname + '_rcpt' + '.dms'
        conn = sqlite.connect(rcptfile)
        q = """SELECT id FROM particle"""
        n_rcpt = 0;
        for atom_id in conn.execute(q):
            n_rcpt += 1
        receptor_sql =  self.keywords.get('REST_LIGAND_CMRECSQL')
        q = "SELECT id FROM particle WHERE %s" % receptor_sql
        rcpt_ids = []
        for atom_id in conn.execute(q):
            rcpt_ids.append( int(atom_id[0]) )
        rcpt_atom_id = str(rcpt_ids)    
        conn.close()
                
        ligfile =  self.jobname + '_lig' + '.dms'
        conn = sqlite.connect(ligfile)
        q = """SELECT id FROM particle"""
        n_lig = 0;
        for atom_id in conn.execute(q):
            n_lig += 1
        ligand_sql =  self.keywords.get('REST_LIGAND_CMLIGSQL')
        q = "SELECT id FROM particle WHERE %s" % ligand_sql
        lig_ids = []
        for atom_id in conn.execute(q):
            lig_ids.append( int(atom_id[0]) )
        lig_atom_id = str(lig_ids)
        conn.close()
        
        #implicit solvent
        if self.keywords.get('IMPLICIT_SOLVENT') is None:
            implicitsolvent = 'AGBNP'
        else:
            implicitsolvent = self.keywords.get('IMPLICIT_SOLVENT')

        #platform etc.
        if self.keywords.get('OPENMM_PLATFORM') is None:
            platform_name = 'Reference'
        else:
            platform_name = self.keywords.get('OPENMM_PLATFORM')

        #force constant etc.
        kf = float(self.keywords.get('REST_LIGAND_CMKF'))
        tol = float(self.keywords.get('REST_LIGAND_CMTOL'))
            
        inputr = mintherm_openmm.format(
            jobname = self.jobname,
            temperature = temperature,
            implicitsolvent = implicitsolvent,
            platform = platform_name,
            nlig_atoms = n_lig,
            rest_ligand_cmrec = rcpt_atom_id,
            rest_ligand_cmlig = lig_atom_id,
            rest_ligand_cmkf = kf,
            rest_ligand_cmtol = tol,
            ligand_ref_atoms = None,
            receptor_ref_atoms = None,
            angle_center = None,
            angletol = None,
            kfangle = None,
            dihedral1center = None,
            kfdihedral1 = None,
            dihedral1tol = None,
            dihedral2center = None,
            kfdihedral2 = None,
            dihedral2tol = None
        )
        
        driverfile = "%s_mintherm.py" % self.jobname
        f = open(driverfile, "w")
        f.write(inputr)
        f.close()

    def writeLigStructureFile(self):
        if self.ligidxfile is None :
            msg = "writeLigStructureFile: Internal error: Structure file not found"
            self.exit(msg)

        if not os.path.exists(self.ligidxfile):
            msg = 'File does not exist: %s' % self.ligidxfile
            self.exit(msg)

        con = sqlite.connect(self.ligidxfile)

        #remove global_cell table: not using periodic boundary conditions
        with con:
            cur = con.cursor()
            cur.execute("DROP TABLE IF EXISTS global_cell")

            
    #override for writing flat-bottom positional restraints into receptor file
    #end removing global cell
    def writeRecStructureFile(self):
        if self.recidxfile is None :
            msg = "writeRecStructureFile: Internal error: Structure file not found"
            self.exit(msg)
            
        if not os.path.exists(self.recidxfile):
            msg = 'File does not exist: %s' % self.recidxfile
            self.exit(msg)

        rest_sql =  self.keywords.get('REST_RECEPTOR_SQL')

        if (rest_sql is None) or (rest_sql is ""):
            return
        
        if self.keywords.get('REST_RECEPTOR_KF') is not None:
            rest_kf = float(self.keywords.get('REST_RECEPTOR_KF'))
        else:
            rest_kf = float('0.6')

        if self.keywords.get('REST_RECEPTOR_TOL') is not None:
            rest_tol = float(self.keywords.get('REST_RECEPTOR_TOL'))
        else:
            rest_tol = float('1.0')

        table_name = "posre_harmflatbottom"
        
        con = sqlite.connect(self.recidxfile)
        with con:
            cur = con.cursor()
            cur.execute("DROP TABLE IF EXISTS %s_term" % table_name)
            cur.execute("DROP TABLE IF EXISTS %s_param" % table_name)
            cur.execute("CREATE TABLE IF NOT EXISTS %s_term (p0 INTEGER PRIMARY KEY, x0 REAL, y0 REAL, z0 REAL, param INTEGER )" % table_name)
            cur.execute("CREATE TABLE IF NOT EXISTS %s_param (id INTEGER PRIMARY KEY, fc REAL, tol REAL)" % table_name)
            cur.execute("INSERT INTO %s_param (fc, tol, id) VALUES (%f, %f, 0)" % (table_name, rest_kf, rest_tol)) 
            cur.execute("SELECT id, x, y, z FROM particle WHERE " + rest_sql)
            rows = cur.fetchall()
            for row in rows:
                atom = row[0]
                x0 = row[1]
                y0 = row[2]
                z0 = row[3]
                cur.execute("INSERT INTO %s_term (p0, x0, y0, z0, param) VALUES (%d, %f, %f, %f, 0)" % (table_name, atom, x0, y0, z0))

            #remove global_cell table: not using periodic boundary conditions
            cur.execute("DROP TABLE IF EXISTS global_cell")            
        
        
##################### MAIN CODE ##########################
if __name__ == '__main__':

    # Setup the logger
    logger = schrodinger.utils.log.get_output_logger("bedam_prep")

    # Parse arguments:
    usage = "%prog [options] <inputfile>"
    parser = cmdline.SingleDashOptionParser(usage)
    (options, args) = parser.parse_args(sys.argv[1:])
    
    if len(args) != 2:
        parser.error("usage= python bedam_workflow_openmm.py <commandfile>")
    
    commandFile = args[0]

    print( "")
    print( "====================================")
    print( "       SDM Job Preparation for OpenMM        ")
    print( "====================================")
    print( "")
    print( "SCHRODINGER: " + os.environ['SCHRODINGER'])
    print( "Started at: " + str(time.asctime()))
    print( "Input file:", commandFile)
    print( "")
    sys.stdout.flush()
    
    print ("Reading options")
    sys.stdout.flush()
    bedam = bedam_job_openmm_asyncre(commandFile, options)
    
    print ("Set put templates for input files ...")
    sys.stdout.flush()
    bedam.setupTemplatesASyncRE()
    
    print ("Analyzing structure files ...")
    sys.stdout.flush()
    bedam.getDesmondDMSFiles()

    print ("Adding atomic restraints etc to structures files ...")
    sys.stdout.flush()
    bedam.writeRecStructureFile()
    bedam.writeLigStructureFile()
    
    print ("Writing mintherm driver file ...")
    sys.stdout.flush()
    bedam.writeThermInputFile()

    print ("Writing job input files ...")
    sys.stdout.flush()
    bedam.writeCntlFile()
    bedam.writeOpenMMDriverFile()
    
    
    print("")
    print( "Job preparation complete")
    print( "")
    print( "When completed run the production calculation with:")
    print( "<path_to_asyncre>/bedamtempt_async_re.py %s_asyncre.cntl" % bedam.jobname)
    print( "")
