# bedam python class
__doc__="""
$Revision: 0.1 $

A class to prepare Temperature AsyncRE jobs

"""
# Contributors: Emilio Gallicchio, Junchao Xia

import os, sys, time, re, glob
from schrodinger.utils import cmdline
import schrodinger.utils.log
import shutil
import signal
import glob
import time

import sqlite3 as sqlite

from math import *

from tempt_asyncre import tempt_job_asyncre

class tempt_job_openmm_asyncre(tempt_job_asyncre):

    #override for writing flat-bottom positional restraints into receptor file
    def writeStructureDMSFile(self):
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



##################### MAIN CODE ##########################
if __name__ == '__main__':

    # Setup the logger
    logger = schrodinger.utils.log.get_output_logger("bedam_prep")

    # Parse arguments:
    usage = "%prog [options] <inputfile>"
    parser = cmdline.SingleDashOptionParser(usage)
    (options, args) = parser.parse_args(sys.argv[1:])
    
    if len(args) != 1:
        parser.error("Please specify ONE input file")
    
    commandFile = args[0]

    print ""
    print "===================================="
    print "       BEDAM Job Preparation        "
    print "===================================="
    print ""
    print "SCHRODINGER: " + os.environ['SCHRODINGER']
    print "Started at: " + str(time.asctime())
    print "Input file:", commandFile
    print ""
    sys.stdout.flush()
    
    print "Reading options"
    tempt = tempt_job_asyncre(commandFile, options)

    print "Set up templates for input files ..."
    tempt.setupTemplatesASyncRE()
    print "Analyzing structure files ..."
    tempt.getDesmondDMSFiles()
    print "Adding atomic restraints to receptor file ..."
    tempt.writeStructureDMSFile()
    print "Writing job input files ..."
    tempt.writeCntlFile()
    tempt.writeThermInputFile()    
    tempt.writeRemdInputFile()
    tempt.writeRunimpactFile()
    print "Write the submission scripts for SuperMIC and Stampede"
    tempt.writeQueueFiles()

    print
    print "Job preparation complete"
    print ""
    print "To run the minimization/thermalization calculation do:"
    exec_directory =  tempt.keywords.get('EXEC_DIRECTORY')
    if exec_directory is not None:
        print "export IMPACTHOME=" + exec_directory
    else:
        print "export IMPACTHOME=" + "<path_to_academic_impact_directory>"
    print "export OMP_NUM_THREADS=" + "<number_of_CPU_cores>"
    print "./runimpact_i %s_mintherm.inp" % tempt.jobname
    print ""
    print "When completed run the production calculation with:"
    print "<path_to_asyncre>/tempt_async_re.py %s_asyncre.cntl" % tempt.jobname
    print ""
