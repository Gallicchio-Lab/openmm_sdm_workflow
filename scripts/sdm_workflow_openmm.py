from __future__ import print_function

# bedam python class
__doc__="""
$Revision: 0.1 $

A class to prepare OpenMM AsyncRE jobs

"""
# Contributors: Emilio Gallicchio

import os, sys, time
import sqlite3 as sqlite
from configobj import ConfigObj

class sdm_job_openmm_asyncre(object):

    def __init__(self, command_file, options):
        self.command_file = command_file
        self.jobname = os.path.splitext(os.path.basename(command_file))[0]
        self.keywords = ConfigObj(self.command_file)

    def exit(self, text):
        """ Print an error and exit """
        print(text)
        sys.exit(1)
        
    #ASyncRE control file
    def writeCntlFile(self):
        input = ""

        job_transport = self.keywords.get('JOB_TRANSPORT')
        if job_transport is None:
            msg = "writeCntlFile: JOB_TRANSPORT is not specified"
            self.exit(msg)
        if not (job_transport == "SSH" or job_transport == "LOCAL_OPENMM"):
            msg = "writeCntlFile: invalid JOB_TRANSPORT: %s Should be SSH or LOCAL_OPENMM." % job_transport
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

        checkpoint_time = self.keywords.get('CHECKPOINT_TIME')
        if checkpoint_time is not None:
            input += "CHECKPOINT_TIME = %d\n" % int(checkpoint_time)
            
        input += "NODEFILE = 'nodefile'\n"

        subjobs_buffer_size = self.keywords.get('SUBJOBS_BUFFER_SIZE')
        if subjobs_buffer_size is not None:
            input += "SUBJOBS_BUFFER_SIZE = '%f'\n" % float(subjobs_buffer_size)

        nsteps = self.keywords.get('PRODUCTION_STEPS')
        nprnt = self.keywords.get('PRNT_FREQUENCY')
        if (not nsteps) or (not nprnt):
            msg = "writeCntlFile: PRODUCTION_STEPS PRNT_FREQUENCY are required"
            self.exit(msg)
        input += "PRODUCTION_STEPS = '%d'\n" % int(nsteps)
        input += "PRNT_FREQUENCY = '%d'\n" % int(nprnt)
        ntrj =  int(self.keywords.get('TRJ_FREQUENCY'))
        if not ntrj:
            ntrj = nprnt
        input += "TRJ_FREQUENCY = '%d'\n" % int(ntrj)

        #implicit solvent
        implicitsolvent = self.keywords.get('IMPLICIT_SOLVENT');        
        if implicitsolvent is None:
            msg = "writeCntlFile: Specify implicit solvent to continue"
            self.exit(msg)
        input += "IMPLICITSOLVENT = '%s'\n" % str(implicitsolvent)

        #if the job_transport is 'LOCAL_OPENMM', add additional information to the asyncre cntl file
        if job_transport == 'LOCAL_OPENMM':
            rcpt_ids = self.getCMRcptAtoms()
            if len(rcpt_ids) > 0:
                rcpt_atom_id = ','.join(map(str,rcpt_ids))
            else:
                rcpt_atom_id = ""
            input += "REST_LIGAND_CMREC_ATOMS = %s\n" % rcpt_atom_id


            lig_ids = self.getCMLigAtoms()
            if len(lig_ids) > 0:
                lig_atom_id = ','.join(map(str,lig_ids))
            else:
                lig_atom_id = ""
            input += "REST_LIGAND_CMLIG_ATOMS = %s\n" % lig_atom_id
            input += "NATOMS_LIGAND = '%d'\n" % int(self.n_lig)
            
            #force constant etc.
            (kf, tol) = self.getVsiteParams()
            if kf is not None:
                input += "CM_KF = %.2f\n" % float(kf)
            if tol is not None:
                input += "CM_TOL = %.2f\n" % float(tol)

            #soft core settings
            (soft_core_method,  soft_core_umax, soft_core_acore) = self.getSoftCoreParams()
            input += "SOFT_CORE_METHOD = '%s'\n" % str(soft_core_method)
            input += "UMAX = %.2f\n" % soft_core_umax
            input += "ACORE = %f\n" % soft_core_acore

            #MD stepsize and relaxation parameters
            (friction_coeff, stepsize) = self.getMDParams()
            input += "FRICTION_COEFF = %f\n" % friction_coeff
            input += "TIME_STEP = %.3f\n" % stepsize

        verbose = self.keywords.get('VERBOSE')
        if verbose is not None:
            input += "VERBOSE = '%s'\n" % verbose

        cntlfile = "%s_asyncre.cntl" % self.jobname
        f = open(cntlfile, "w")
        f.write(input)
        f.close

    def getCMRcptAtoms(self):
        rcptfile = self.jobname + '_rcpt' + '.dms'
        conn = sqlite.connect(rcptfile)
        receptor_sql =  self.keywords.get('REST_LIGAND_CMRECSQL')
        if not receptor_sql:
            self.exit("REST_LIGAND_CMRECSQL is required")
        q = "SELECT id FROM particle WHERE %s" % receptor_sql
        rcpt_ids = []
        rcpt_atom_id = ""
        for atom_id in conn.execute(q):
            rcpt_ids.append( int(atom_id[0]) )
        conn.close()
        return rcpt_ids

    def getCMLigAtoms(self):
        ligfile =  self.jobname + '_lig' + '.dms'
        conn = sqlite.connect(ligfile)
        q = """SELECT id FROM particle"""
        self.n_lig = 0;
        for atom_id in conn.execute(q):
            self.n_lig += 1
        ligand_sql =  self.keywords.get('REST_LIGAND_CMLIGSQL')
        if not ligand_sql:
            self.exit("REST_LIGAND_CMLIGSQL is required")
        q = "SELECT id FROM particle WHERE %s" % ligand_sql
        lig_ids = []
        for atom_id in conn.execute(q):
            lig_ids.append( int(atom_id[0]) )
        conn.close()
        return lig_ids

    def getVsiteParams(self):
        kf = float(self.keywords.get('REST_LIGAND_CMKF'))
        tol = float(self.keywords.get('REST_LIGAND_CMTOL'))
        return (kf, tol)

    def getSoftCoreParams(self):
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
        return (soft_core_method, soft_core_umax, soft_core_acore)

    def getMDParams(self):
        friction_coeff = self.keywords.get('FRICTION_COEFF')
        if friction_coeff is None:
            friction_coeff = 0.5 #in picosecond
        else:
            friction_coeff = float(friction_coeff)
        stepsize = self.keywords.get('MD_TIME_STEP')
        if stepsize is None:
            stepsize = 0.001 #in picosecond
        else:
            stepsize = float(stepsize)
        return (friction_coeff, stepsize)
    
#the template driver file for each MD cycle
    def writeOpenMMDriverFile(self):
        job_transport = self.keywords.get('JOB_TRANSPORT')
        if job_transport != 'SSH':
            return

        with open('openmm-driver-template.py', 'r') as f:
            input_openmm = f.read()

        #Vsite settings
        rcpt_ids = self.getCMRcptAtoms()
        if len(rcpt_ids) > 0:
            rcpt_atom_id = str(rcpt_ids)
        else:
            rcpt_atom_id = ""

        lig_ids = self.getCMLigAtoms()
        if len(lig_ids) > 0:
            lig_atom_id = str(lig_ids)
        else:
            lig_atom_id = ""

        (kf, tol) = self.getVsiteParams()

        #MD parameters
        (friction_coeff, stepsize) = self.getMDParams()
        
        #steps etc.
        nsteps = int(self.keywords.get('PRODUCTION_STEPS'))
        nprnt =  int(self.keywords.get('PRNT_FREQUENCY'))
        ntrj =  int(self.keywords.get('TRJ_FREQUENCY'))

        #implicit solvent
        implicitsolvent = self.keywords.get('IMPLICIT_SOLVENT');        
        if implicitsolvent is None:
            implicitsolvent = None
        
        #soft core settings
        (soft_core_method,  soft_core_umax, soft_core_acore) = self.getSoftCoreParams()

        inputr = input_openmm.format(
            jobname = self.jobname,
            nlig_atoms = self.n_lig,
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
        #only write Python driver file if job trasnport is 'SSH'
        f = open(driverfile, "w")
        f.write(inputr)
        f.close()

#minimization/thermalization .py driver file
    def writeThermInputFile(self):
        with open('openmm-mintherm-template.py', 'r') as f:
            mintherm_openmm = f.read()

        temperature = self.keywords.get('TEMPERATURES').split(",")[0]

        #Vsite settings
        rcpt_ids = self.getCMRcptAtoms()
        if len(rcpt_ids) > 0:
            rcpt_atom_id = str(rcpt_ids)
        else:
            rcpt_atom_id = ""
        lig_ids = self.getCMLigAtoms()
        if len(lig_ids) > 0:
            lig_atom_id = str(lig_ids)
        else:
            lig_atom_id = ""
        (kf, tol) = self.getVsiteParams()
        
        #implicit solvent
        implicitsolvent = self.keywords.get('IMPLICIT_SOLVENT')
        
        #platform etc.
        if self.keywords.get('OPENMM_PLATFORM') is None:
            platform_name = 'Reference'
        else:
            platform_name = self.keywords.get('OPENMM_PLATFORM')

        inputr = mintherm_openmm.format(
            jobname = self.jobname,
            temperature = temperature,
            implicitsolvent = implicitsolvent,
            platform = platform_name,
            nlig_atoms = self.n_lig,
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
        ligfile =  self.jobname + '_lig' + '.dms'
        con = sqlite.connect(ligfile)
        #remove global_cell table: not using periodic boundary conditions
        with con:
            cur = con.cursor()
            cur.execute("DROP TABLE IF EXISTS global_cell")


    def writeUWHAMFile(self):
        with open('uwham_analysis-template.R', 'r') as f:
            input_uwham = f.read()

        temperatures = self.keywords.get('TEMPERATURES')
        if temperatures is None:
            msg = "writeUWAHMFile: 'TEMPERATURES' is required."
            self.exit(msg)
        
        lambdas = self.keywords.get('LAMBDAS')
        if lambdas is None:
            msg = "writeUWAHMFile: 'LAMBDAS' is required."
            self.exit(msg)
            
        nlambdas = len(lambdas.split(","))
        zerosdefault = "0.000"
        for i in range(1,nlambdas):
            zerosdefault += ",0.000"

        lambda1 = self.keywords.get('LAMBDA1')
        if lambda1 is None:
            lambda1 = lambdas
        
        lambda2 = self.keywords.get('LAMBDA2')
        if lambda2 is None:
            lambda2 = lambdas

        alpha = self.keywords.get('ALPHA')
        if alpha is None:
            alpha = zerosdefault
        
        u0 = self.keywords.get('U0')
        if u0 is None:
            u0 = zerosdefault
            
        w0coeff = self.keywords.get('W0COEFF')
        if w0coeff is None:
            w0coeff = zerosdefault

        vsiterad = self.keywords.get('REST_LIGAND_CMTOL')
        if vsiterad is None:
            msg = "writeUWAHMFile: 'REST_LIGAND_CMTOL' is required."
            self.exit(msg)
        vsiterad = float(vsiterad)
            
        #can't use .format() because of curly brackets in R, order is important
        inputr = input_uwham % (
            lambdas,
            lambda1,
            lambda2,
            alpha,
            u0,
            w0coeff,
            temperatures,
            vsiterad
        )
            
        uwham_file = "uwham_analysis.R"
        f = open(uwham_file, "w")
        f.write(inputr)
        f.close()
            
##################### MAIN CODE ##########################
if __name__ == '__main__':
    # Parse arguments:
    usage = "%prog <ConfigFile>"
    
    if len(sys.argv) != 2:
        print("Please specify ONE input file")
        sys.exit(1)

    commandFile = sys.argv[1]

    sys.stdout.flush()
    
    print ("Reading options")
    sys.stdout.flush()
    sdm = sdm_job_openmm_asyncre(commandFile, options=None)
    
    print ("Writing mintherm driver file ...")
    sys.stdout.flush()
    sdm.writeThermInputFile()

    print ("Writing job input files ...")
    sys.stdout.flush()
    sdm.writeCntlFile()
    sdm.writeOpenMMDriverFile()

    print("Writing analysis UWHAM R script")
    sys.stdout.flush()
    sdm.writeUWHAMFile()
    
