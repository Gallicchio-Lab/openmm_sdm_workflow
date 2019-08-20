from __future__ import print_function

# add positional flat bottom harmonic restraints to a .dms file
usage = """
python add_fbh_restraints -r <sqlite_selection> -t <tolerance> -k <force constant> -f <dms file>

<sqlite selection>: an sqlite selection string such as "name GLOB 'CA'" for the target atoms    
<tolerange>: tolerance in angstroms
<force constant>: force constant in kcal/mol/angstrom^2
<dms file>: a .dms file
"""
# Contributors: Emilio Gallicchio <emilio.gallicchio@gmail.com>

import os, sys, getopt
import sqlite3 as sqlite

try:
    opts, args = getopt.getopt(sys.argv[1:],"hr:t:k:f:")
except getopt.GetoptError:
    print(usage)
    sys.exit(1)
    
if len(opts) != 4:
    print(usage)
    sys.exit(1)
    
for opt, arg in opts:
    if opt == '-r':
        sqlite_selection = arg
    elif opt == '-t':
        rest_tol = float(arg)
    elif opt == '-k':
        rest_kf = float(arg)
    elif opt == '-f':
        dms_file = arg
    else:
        print(usage)
        sys.exit(1)
        
if not os.path.exists(dms_file):
    print('%s does not exist' % dms_file)
    sys.exit(1)

if (sqlite_selection is None) or (sqlite_selection is ""):
    sys.exit(0)
        
table_name = "posre_harmflatbottom"
        
con = sqlite.connect(dms_file)
with con:
    cur = con.cursor()
    cur.execute("DROP TABLE IF EXISTS %s_term" % table_name)
    cur.execute("DROP TABLE IF EXISTS %s_param" % table_name)
    cur.execute("CREATE TABLE IF NOT EXISTS %s_term (p0 INTEGER PRIMARY KEY, x0 REAL, y0 REAL, z0 REAL, param INTEGER )" % table_name)
    cur.execute("CREATE TABLE IF NOT EXISTS %s_param (id INTEGER PRIMARY KEY, fc REAL, tol REAL)" % table_name)
    cur.execute("INSERT INTO %s_param (fc, tol, id) VALUES (%f, %f, 0)" % (table_name, rest_kf, rest_tol)) 
    cur.execute("SELECT id, x, y, z FROM particle WHERE " + sqlite_selection)
    rows = cur.fetchall()
    for row in rows:
          atom = row[0]
          x0 = row[1]
          y0 = row[2]
          z0 = row[3]
          cur.execute("INSERT INTO %s_term (p0, x0, y0, z0, param) VALUES (%d, %f, %f, %f, 0)" % (table_name, atom, x0, y0, z0))
