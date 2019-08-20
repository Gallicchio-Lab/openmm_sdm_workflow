"""
Add a table with agbnp2 parameters to a .dms file.

Copyright Schrodinger, LLC. All rights reserved.
"""
# Author: ivan.tubert-brohman@schrodinger.com

"""
usage: $SCHRODINGER/run add_agbnp2.py [options] <dms_file>

Adds a table with agbnp2 parameters to the .dms file specified in the command
line.

positional arguments:
  dms_file              .dms file to parameterize

optional arguments:
  -v, -version          Show the program's version and exit.
  -h, -help             Show this help message and exit.
  -param_file PARAM_FILE
                        path to parameter file (default: 'agbnp2.param')
"""

import sys, os, math, sqlite3
from schrodinger import structure
from schrodinger.infra import mm
from schrodinger.utils import cmdline
from schrodinger.structutils.analyze import evaluate_smarts

TABLE_NAME = 'agbnp2';

# these are the parameters *in the order in which they appear* in the
# agbnp2.param file. They must also match the names of the columns in the
# SQL table.
PARAM_NAMES = \
    'radius igamma ialpha idelta sgamma salpha sdelta hbtype hbw'.split()

FORMAL_CHARGE_PROP = 'i_m_formal_charge'
ID_PROP = 'i_dms_id'
SIGMA_PROP = 'r_dms_sigma'
EPSILON_PROP = 'r_dms_epsilon'

# Constants from agbnp.c
SIGMAW = 3.15365 # LJ sigma of TIP4P water oxygen
EPSILONW = 0.155 # LJ epsilon of TIP4P water oxygen
RHO = 0.033428   # water number density

class Agbnp2Error(Exception):
    """A class for exceptions raised by this module."""
    pass

def read_structure(conn):
    """
    Read a structure from the particle and bond tables of a database
    connection. Returns a schrodinger.structure.Structure object.
    The atoms get properties added with their database id and their
    sigma and epsilon parameters.
    """
    st = structure.create_new_structure()
    c = conn.cursor()
    atoms = {}
    for id, anum, x, y, z, sigma, epsilon, formal_charge in c.execute(
            "select particle.id, anum, x, y, z, sigma, epsilon, formal_charge "
            "from particle, nonbonded_param "
            "where particle.nbtype==nonbonded_param.id"):
        el = mm.mmat_get_element_by_atomic_number(anum)
        atom = atoms[id] = st.addAtom(el, x, y, z)
        atom.property[FORMAL_CHARGE_PROP] = formal_charge
        atom.property[ID_PROP] = id
        atom.property[SIGMA_PROP] = sigma
        atom.property[EPSILON_PROP] = epsilon
    for p0, p1, order in c.execute('select p0, p1, "order" from bond'):
        st.addBond(atoms[p0], atoms[p1], int(order))
    return st

def read_param_file(fname):
    """
    Read an agbnp2.param file. Returns a list of atom types, where each atom
    type is a tuple (smarts, params_dict).
    """
    # Header and sample line:
    #type  atom  radius  igamma ialpha   idelta    sgamma   salpha       sdelta  HB H
    #  28 [#1][C;X4]               1.250  0.000    0.800    0.000  0.00000  0.00000    0.00000      0  0.0
    atom_types = []
    for line in open(fname):
        if line.startswith('#') or line.startswith('dielectric'): continue
        cols = line.split()
        smarts = cols[1]
        params_list = [float(s) for s in cols[2:]]
        if len(params_list) != len(PARAM_NAMES):
            raise Agbnp2Error("invalid number of parameters in line:\n %s"
                % line)
        params_dict = dict(zip(PARAM_NAMES, params_list))
        atom_types.append((smarts, params_dict))
    return atom_types
        
def add_agbnp2_table(conn, atom_types, table_name):
    """
    Create or replace a table with the agbnp2 parameters for the atoms in the
    structure described by the particle and bond tables of the database
    conection.
    """
    st = read_structure(conn)
    c = conn.cursor()
    c.execute("drop table if exists %s" % table_name)
    c.execute("CREATE TABLE %s (id INTEGER NOT NULL REFERENCES particle, "
        "radius FLOAT,igamma FLOAT,ialpha FLOAT,idelta FLOAT, "
        "sgamma FLOAT,salpha FLOAT,sdelta FLOAT,hbtype INTEGER,hbw FLOAT)"
        % table_name)
    atoms = {}
    for smarts, params in atom_types:
        for match in evaluate_smarts(st, smarts):
            iatom = match[0]
            if iatom not in atoms:
                atom = atoms[iatom] = params.copy()
                atom['id'] = st.atom[iatom].property[ID_PROP]
                sigma = st.atom[iatom].property[SIGMA_PROP]
                epsilon = st.atom[iatom].property[EPSILON_PROP]
                # NOTE: we are hard-coding the OPLS combining rule here.
                # I hope that's a reasonable assumption.
                epsiloniw = (epsilon*EPSILONW)**0.5
                s6 = (sigma*SIGMAW)**3
                # alpha scaling formula adapted from agbnp.c
                ws = -16.0 * math.pi * RHO * epsiloniw * s6 / 3.0
                atom['ialpha'] *= ws
                atom['salpha'] *= ws
    for iatom in xrange(1, len(st.atom)+1): 
        if iatom not in atoms:
            raise Agbnp2Error("parameters not found for atom %d" % iatom)
        cols = ['id'] + PARAM_NAMES
        cmd = "INSERT into %s (%s) values(%s)" % (table_name, 
            ','.join(cols),
            ','.join([':%s' % name for name in cols]))
        c.execute(cmd, atoms[iatom])
    conn.commit()

def add_agbnp2_to_dms_file(dms_file, param_file):
    """
    Create or replace a table with the agbnp2 parameters for the atoms in the
    structure contained in dms_file, using the parameters specified in
    param_file.
    """
    if not os.path.isfile(dms_file):
        raise Agbnp2Error("file does not exist: %s" % dms_file)
    conn = sqlite3.connect(dms_file)
    atom_types = read_param_file(param_file)
    add_agbnp2_table(conn, atom_types, TABLE_NAME)
    conn.close()

def parse_args():
    parser = cmdline.create_argument_parser(
        usage="$SCHRODINGER/run add_agbnp2.py [options] <dms_file>",
        description="Adds a table with agbnp2 parameters to the .dms file "
            "specified in the command line.")
    parser.add_argument("dms_file", help=".dms file to parameterize")
    parser.add_argument("-param_file",
        help="path to parameter file (default: 'agbnp2.param')",
        default="agbnp2.param")
    return parser.parse_args()

def main():
    args = parse_args()
    try:
        add_agbnp2_to_dms_file(args.dms_file, args.param_file)
    except (Agbnp2Error, IOError) as e:
        print >> sys.stderr, "Error: %s" % e
        sys.exit(1)

if __name__ == '__main__':
    main()

