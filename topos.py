#!/usr/bin/env python
import os, sys, re, math
import numpy as np
import time
from datetime import date
from reference import quotes as q
from modify_structure import man

# Writes topology files in different formats, including:
# lammps data file format
# pdb file
# xyz file

# Converts atom ID number to symbol
def ids_symbols(ids):
   
    atom_ids = {
    	1 : 'H',
    	2 : 'O',
    	3 : 'C',
    	4 : 'S',
    	5 : 'Na',
	}
    symbols = [] 
    for i in ids:    
        symbols.append(atom_ids[i])
    symbols =  np.array(symbols).reshape(len(ids))
        
    return symbols
    
def symbols_ids(symbols,new_symbols={}):
   
    atom_ids = {
        'H' : 1 ,
        'O' : 2 ,
	}
    for key, val in new_symbols.items():
        atom_ids[key] = val
    ids = [] 
    for i in symbols:    
        ids.append(atom_ids[i])
    ids =  np.array(ids).reshape(len(ids))
        
    return ids

def read_lammps(
               fn_in
               ):

    data = {

        # Cell parameters
        "xhi" : 0.0,
        "yhi" : 0.0,
        "zhi" : 0.0,  
        "xlo" : 0.0,
        "ylo" : 0.0,
        "zlo" : 0.0,  

        # Atom parameters   
        "atoms" : 0,
        "bonds" : 0,
        "angles" : 0,
        "dihedrals" : 0,
        "impropers" : 0,
        "atom_types" : 2,
        "bond_types" : 0, 
        "angle_types" : 0,
        "dihedral_types" : 0,
        "improper_types" : 0,

        # Masses
        "masses" : [],

        # Symbols (Order should correspond to order of Masses)
        "symbols" : ['H','O'],
        "ids" : [],
    }

    f = open(fn_in,'r')
    header1 = f.readline()
    header2 = f.readline()
   
    natoms = 0 
    ndata = len(data)

    index = []
    molecule_tag = []
    atom_type = []
    charge = []
    coordinates = []
    velocities = []
    nxyz = []
    labels = []

    for line in f:
        line = line.strip()
        columns = line.split()

        for ind,col in enumerate(columns):
            if '#' in col:
                columns = np.delete(columns,obj=ind)

        if len(columns) == 0:
            continue
        elif len(columns) == 1:
            labels.append(columns[0])

        elif len(columns) == 2:
            if len(labels) >= 1:
                if labels[-1] == "Masses":
                    data["ids"] = np.hstack((int(columns[0]),data["ids"]))
                    data["masses"] = np.hstack((np.array(columns[1],dtype=float),data["masses"]))
            else:
                data[columns[1]] = int(columns[0])
            if natoms == 0:
                natoms = int(data["atoms"])

        elif len(columns) == 3 and columns[2] == "types":
            key = "%s_types" % (columns[1])
            data[key] = int(columns[0])

        elif len(columns) == 4:
            if 'lo' in columns[2] and 'hi' in columns[3]:
                data[columns[2]] = float(columns[0])
                data[columns[3]] = float(columns[1])
            elif len(labels) >= 1:
                if labels[-1] == "Velocities":
                    velocities.append(columns[1:4])

        elif len(columns) >= 7 and labels[-1] == "Atoms":
            index.append(columns[0])
            molecule_tag.append(columns[1])
            atom_type.append(columns[2])
            charge.append(columns[3])
            coordinates.append(columns[4:7])
            if len(columns) == 10:
                nxyz.append(columns[7:-1])
        else:
            labels.append(columns[0])
            print "WARNING: line ignored. Data not registered in lammps-file dictionary:" 
            print columns

        if len(data) != ndata:
            print "WARNING: Length of lammps dictionary changed!"
            print columns

        ndata = len(data)

    if len(coordinates) == 0:
        print "WARNING: No coordinate data!"
    if len(coordinates) != natoms:
        print "WARNING: Coordinate data does not equal number of atoms"

    index = np.array(index,dtype=int).reshape(natoms)
    molecule_tag = np.array(molecule_tag,dtype=int).reshape(natoms)
    atom_type = np.array(atom_type,dtype=int).reshape(natoms)
    charge = np.array(charge,dtype=float).reshape(natoms)
    coordinates = np.array(coordinates,dtype=float).reshape(natoms,3)
    if len(velocities) == len(coordinates):
        velocities = np.array(velocities,dtype=float).reshape(natoms,3)

    return data, index, molecule_tag, atom_type, charge, coordinates, velocities

def read_lammps_dump(
               fn_in = 'dump.water',
               fn_out = 'charges',
               len_header = 9,
               natoms = 301,
               timestep=100,
               nsteps=10000000,
               ):
   
    nframes = nsteps/timestep
    #nframes = 645
    nlines = (nframes+1)*(natoms+len_header)

    #charges = np.zeros((nframes,natoms))
    positions = np.zeros((nframes,natoms,3))

    skip_header = len_header
    #charges= np.genfromtxt(fn_in,dtype=float,skip_header=skip_header,usecols=7,invalid_raise=False)
    positions= np.genfromtxt(fn_in,dtype=float,skip_header=skip_header,usecols=(1,2,3),invalid_raise=False)
    #charges = charges[~np.isnan(charges)].reshape(-1,natoms)
    #
    #charges.dump(fn_out)
    positions.dump(fn_out)

    return 

def dump2xyz(
             symbols,
             fn_in = 'dump.water',
             fn_out = 'dump.xyz',
             len_header = 9,
             natoms = 901,
             timestep=100,
   #          nsteps=10000000,
    #         nframes = 645,
             ):
   

    skip_header = len_header
    coordinates= np.genfromtxt(fn_in,dtype=float,skip_header=skip_header,usecols=[1,2,3],invalid_raise=False)
    coordinates = coordinates[~np.isnan(coordinates)].reshape(-1,natoms,3)
    nframes = len(coordinates)
    nlines = (nframes+1)*(natoms+len_header)

    hs = open(fn_out,"w")
    for ind in range(nframes/timestep):
        hs.write('%d\n\n' % (natoms))
        for atom in range(natoms):
            hs.write('%s %3.3f %3.3f %3.3f\n' % (symbols[atom],coordinates[ind,atom,0],coordinates[ind,atom,1],coordinates[ind,atom,2]))
        print "step: %d/%d" % (ind,nframes)

    hs.close()
    return 

def write_lammps(
          fn_out,
          index,
          atom_type,
          charge,
          coordinates,
          velocities = 0,
          molecule_tag = [0],
          nx = 0,
          ny = 0,
          nz = 0,
          in_opts = {},
          ):

    options = {


        # Cell parameters
        "xhi" : 0.0,
        "yhi" : 0.0,
        "zhi" : 0.0,  
        "xlo" : 0.0,
        "ylo" : 0.0,
        "zlo" : 0.0,  

        # Atom parameters   
        "atoms" : 0,
        "bonds" : 0,
        "angles" : 0,
        "dihedrals" : 0,
        "impropers" : 0,
        "atom_types" : 2,
        "bond_types" : 0, 
        "angle_types" : 0,
        "dihedral_types" : 0,
        "improper_types" : 0,

        # Masses
        "masses" : [1.007940,15.999400],

        # Symbols (Order should correspond to order of Masses)
        "symbols" : ['H','O'],
        "ids" : [1,2],

    }

    # => Override Default Options <= #

    for key,val in in_opts.items():
       if key not in options.keys():
          raise ValueError('%s option unavailable, RTFM' % (key))
       #if options[key] != val:
       #   print "Default Option Overridden: %s = %s" % (key,str(val))
       options[key] = val
    print '\n'

    # => Declaring Variables from Options <= #

    xhi = options['xhi'] 
    yhi = options['yhi']
    zhi = options['zhi']
    xlo = options['xlo']
    ylo = options['ylo']
    zlo = options['zlo']
    
    natoms          = options["atoms"]
    bonds          = options["bonds"]
    angles         = options["angles"]
    dihedrals      = options["dihedrals"]
    impropers      = options["impropers"]
    atom_types     = options["atom_types"]
    bond_types     = options["bond_types"]
    angle_types    = options["angle_types"]
    dihedral_types = options["dihedral_types"]
    improper_types = options["improper_types"]
    masses = options["masses"]
    symbols = options["symbols"]
    ids = options["ids"]

    natoms = len(coordinates)

    if xlo == 0.0 and np.min(coordinates[:,0]) < 0.0:
        xlo = np.min(coordinates[:,0]) - 10.0
    if ylo == 0.0 and np.min(coordinates[:,1]) < 0.0:
        ylo = np.min(coordinates[:,1]) - 10.0
    if zlo == 0.0 and np.min(coordinates[:,2]) < 0.0:
        zlo = np.min(coordinates[:,2]) - 10.0

    #if xlo == 0.0 and np.min(coordinates[:,0]) < 10.0:
       # coordinates = man.displace(coordinates,disp=[10.0,0.0,0.0]) 
    #if ylo == 0.0 and np.min(coordinates[:,1]) < 10.0:
       # coordinates = man.displace(coordinates,disp=[0.0,10.0,0.0]) 
    #if zlo == 0.0 and np.min(coordinates[:,2]) < 10.0:
       # coordinates = man.displace(coordinates,disp=[0.0,0.0,10.0]) 

    if xhi == 0.0:
        xhi = np.max(coordinates[:,0]) + 10.0
    if yhi == 0.0:
        yhi = np.max(coordinates[:,1]) + 10.0
    if zhi == 0.0:
        zhi = np.max(coordinates[:,2]) + 10.0
        
    timestamp = "LAMMPS data file. Written: %s\n\n" % (str(date.today()))
    hs = open(fn_out,"w")
    hs.write(timestamp)
    hs.write("%d atoms\n" % (natoms))
    hs.write("%d bonds\n" % (bonds))
    hs.write("%d angles\n" % (angles))
    hs.write("%d dihedrals\n" % (dihedrals))
    hs.write("%d impropers\n" % (impropers))
    hs.write("%d atom types\n" % (atom_types))
    hs.write("%d bond types\n" % (bond_types))
    hs.write("%d angle types\n" % (angle_types))
    hs.write("%d dihedral types\n" % (dihedral_types))
    hs.write("%d improper types\n" % (improper_types))
    hs.write("%4.4f %4.4f xlo xhi\n" % (xlo,xhi))
    hs.write("%4.4f %4.4f ylo yhi\n" % (ylo,yhi))
    hs.write("%4.4f %4.4f zlo zhi\n\n" % (zlo,zhi))
    hs.write("Masses\n\n")
 
    for ind,mass in enumerate(masses):
        hs.write("%d %2.6f #%s\n" % (ind+1,mass,symbols[ind]))
    hs.write("\n\n")
    hs.write("Atoms\n\n")
       
    if len(molecule_tag) == 1:
        molecule_tag = np.ones(natoms) 

    if nx == 0:
        nx = np.zeros(natoms)
    if ny == 0:
        ny = np.zeros(natoms)
    if nz == 0:
        nz = np.zeros(natoms)

    #TODO: Separate the different components of the atoms so these indices do not depend on formatting!!!! 
    for i in range(natoms):
        hs.write("%d %1d %1d %.16e %.16e %.16e %.16e %1d %1d %1d\n" % (index[i], 
                                                                      molecule_tag[i],
                                                                      atom_type[i], 
                                                                      charge[i],
                                                                      coordinates[i,0], 
                                                                      coordinates[i,1], 
                                                                      coordinates[i,2], 
                                                                      nx[i], 
                                                                      ny[i], 
                                                                      nz[i]))
    hs.write("\n")
    if velocities != 0:
        hs.write("Velocities\n\n")
        for i in range(natoms):
            hs.write("%d %.16e %.16e %.16e\n") % (index[i], velocities[i,0], velocities[i,1],velocities[i,2]) 

    hs.close()

    return

def write_pdb(
              fn_out,
              coordinates,
              index,
              symbols,
              in_opts = {}
              ):

    options = {
        # Cell parameters
        "xhi" : 0.0,
        "yhi" : 0.0,
        "zhi" : 0.0,  
    }    
    # => Override Default Options <= #

    for key,val in in_opts.items():
       if key not in options.keys():
          raise ValueError('%s option unavailable, RTFM' % (key))
       #if options[key] != val:
       #   print "Default Option Overridden: %s = %s" % (key,str(val))
       options[key] = val
    print '\n'

    natoms = len(coordinates)
    res_name = 'MOL'
    chain_id = 'X'
    seq_no = '1'
    occupancy = 0.0
    temp = 0.0
    seg_id = 1

    alpha = 90.0
    beta = 90.0
    gamma = 90.0

    xhi = options['xhi'] 
    yhi = options['yhi']
    zhi = options['zhi']

    if xhi == 0.0:
        xhi = np.max(coordinates[:,0])
    if yhi == 0.0:
        yhi = np.max(coordinates[:,1])
    if zhi == 0.0:
        zhi = np.max(coordinates[:,2])

    fh = open(fn_out,"w")
    header = "CRYST1"
    header += "%9.3f" % (xhi)
    header += "%9.3f" % (yhi)
    header += "%9.3f" % (zhi)
    header += "%7.2f" % (alpha)
    header += "%7.2f" % (beta)
    header += "%7.2f" % (gamma)

    fh.write("%s\n\n" % (header))
    for n in range(natoms):
        s = "ATOM  " #1-6
        s += "%5d " % (index[n]) #7-12
        s += "%-4s " % (symbols[n]) #13-17
        s += "%3s " % (res_name) #18-21
        s += "%1s" % (chain_id) #22
        s += "%4s" % (seq_no) #23-26
        s += "    " #27-30
        s += "%8.3f" %(coordinates[n,0]) #31-38
        s += "%8.3f" %(coordinates[n,1]) 
        s += "%8.3f" %(coordinates[n,2])
        s += "%6.2f " %(occupancy)
        s += "%6.2f      " %(temp)
        s += "%-1s" %(seg_id)
        s += "%2s" %(symbols[n])
        s += "\n"
        fh.write(s)
    fh.close()


    return

def write_xyz(
    fn_out,
    atoms,
    symbols,
    append = False
    ):

    if append == False:
        fh = open(fn_out,"w")
    else:
        fh = open(fn_out,"a")

    if len(symbols) == 1:
        symbols *= atoms.shape[0]

    fh.write("%d\n\n" % (atoms.shape[0]))
    for n in range(atoms.shape[0]):
        fh.write("%s %14.6f %14.6f %14.6f\n" % (symbols[n],atoms[n,0],atoms[n,1],atoms[n,2]))
    fh.close()

    return
