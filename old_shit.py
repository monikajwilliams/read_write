import os, sys, re, math
import numpy as np
from numpy import linalg as LA
import com as CM
from operator import itemgetter
import mdtraj as mdt
import time
from datetime import date
import quotes as q
import copy as c

atomIDs = {
	1 : 'H',
	2 : 'O',
	}
revatomIDs = {
	'H' : 1,
	'O' : 2,
	}
boxAngle = 90.00


#<==============Reads Lammps data File================>#
def readFile(filename):
	fh = open(filename, "r")
	lines = fh.readlines()
	vals = []
	natoms = 0
	for line in lines:
		atoms1= re.match('^\s*(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s*$', line)
		if atoms1 != None:
			natoms += 1
			vals.append([
				float(atoms1.group(1)),
				float(atoms1.group(2)),
				float(atoms1.group(3)),
				float(atoms1.group(4)),
				float(atoms1.group(5)),
				float(atoms1.group(6)),
				float(atoms1.group(7)),
				float(atoms1.group(8)),
				float(atoms1.group(9)),
				float(atoms1.group(10)),
				])
	atomsList = sorted(vals, key=itemgetter(0))
	atoms = np.array(atomsList)
	q.h(0)
	return atoms
		
#<==============Reads PDB File================>#
def readPDB(filename):
	fh = open(filename, "r")
	lines = fh.readlines()
	vals = []
	natoms = 0
	for line in lines:
		atoms1= re.match('^\s*(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s*$', line)
		if atoms1 != None:
			natoms += 1
			vals.append([
				float(atoms1.group(2)),
				float(atoms1.group(6)),
				atoms1.group(3),
				float(atoms1.group(6)),
				float(atoms1.group(7)),
				float(atoms1.group(8)),
				float(atoms1.group(9)),
				float(atoms1.group(11)),
				float(atoms1.group(11)),
				float(atoms1.group(11)),
				])
	atomsList = sorted(vals, key=itemgetter(1))
	atoms = np.array(atomsList)
	atomsN = len(atoms)
	for i in range(atomsN):
		atoms[i,2] = (revatomIDs[atoms[i,2]])
	atoms1 = atoms.astype(np.float)
	q.h(3)
	return atoms1
		
#<==============Add Rows================>#
def addRow(
    atoms,
    nXRows=0.0, 
    nYRows=0.0,
    Xspace=12.0,
    Yspace=12.0, 
    atomsMol=3,
    ):
    initLen = len(atoms)
    newAtoms = np.zeros(((initLen)*(nXRows+1)*(nYRows+1),10))
    for i in range(nXRows+1):
	 newAtoms[initLen*(i):initLen*(i)+initLen,:] += c.deepcopy(atoms[:,:])
	 newAtoms[initLen*(i):initLen*(i)+initLen,4] += Xspace*i
	 newAtoms[initLen*(i):initLen*(i)+initLen,0] += initLen*i
    if nYRows != 0.0: 
  	 for i in range(nYRows):
		a = initLen*(nXRows+1)
		b = a + initLen*(nXRows+1)*i
		d = b + initLen*(nXRows+1)
		newAtoms[b:d,:] += c.deepcopy(newAtoms[:a,:])
	 	newAtoms[b:d,5] += Yspace*(i)
	 	newAtoms[b:d,0] += b

    q.dn(0) 
    return newAtoms 	

#<==============Add Molecules================>#
def addMol(atoms,nadd, atomsMol=3):
    atomsN = len(atoms)
    oDist = atoms[0,6] - atoms[atomsMol+1,6]
    offSet = max(atoms[:,6])
    molTemplate = c.deepcopy(atoms[0:atomsMol,:])
    newAtoms = np.zeros((nadd*atomsMol,10))
    m = 0
    for n in range(0,nadd*atomsMol,atomsMol):
	moleculeID = nadd*atomsMol - n
	k = 0
	m += 1
	for i in range(n, n+atomsMol):
		atomID = nadd*atomsMol - (i+1)
		newAtoms[atomID,:] += c.deepcopy(molTemplate[k,:])
		newAtoms[atomID,0] = atomID + 1
		newAtoms[atomID,6] += oDist*m + offSet
		k+=1
    atoms[:,0] += nadd*atomsMol
    finalAtoms = np.append(newAtoms, atoms, axis=0)
    return finalAtoms
     

#<==============Displace wire in x,y,z direction================>#
def displaceWire(atoms, xhi=0.0, yhi=0.0,displacementZ=0.0):
        if xhi != 0.0:
            centerX = xhi/2
	    displacementX = centerX-atoms[1,4]	
	    atoms[:,4] += displacementX
        if yhi != 0.0:
	    centerY = yhi/2
	    displacementY = centerY-atoms[1,5]	
	    atoms[:,5] += displacementY
	atoms[:,6] += displacementZ
	
	return atoms

#<==============Displace wire in x,y,z direction================>#
def displace(atoms, xdiff=0.0, ydiff=0.0,zdiff=0.0):
	atoms[:,4] += xdiff
	atoms[:,5] += ydiff
	atoms[:,6] += zdiff
	
	return atoms

#<==============Displace wire in x,y,z direction================>#
def new_displace(atoms, disp=[0.0,0.0,0.0]):

        new_atoms = np.copy(atoms)
	new_atoms[:,4] += disp[0]
	new_atoms[:,5] += disp[1]
	new_atoms[:,6] += disp[2]
	
	return new_atoms

def expand(atoms, axis, times,gap = 2.0):
    
    natoms = len(atoms)

    finalatoms = np.copy(atoms)
    for ind,n in enumerate(axis):
        disp = [0.0,0.0,0.0]
        newatoms = np.copy(finalatoms)
        for tind in range(times[ind]):
            disp[n-4] = (tind+1)*(np.max(newatoms[:,n]) + gap) 
            new_set = new_displace(newatoms,disp)
            finalatoms = np.vstack((finalatoms,new_set))

    finalatoms = index(finalatoms)

    return finalatoms
    

#<==============Alter the spacing between molecules in the z direction================>#
def constrict(atoms, density,ax_index=6,atomsMol=3):
      
        # Axis indices:
        # x = 4
        # y = 5
        # z = 6
 
        indsO = [ind for ind, val in enumerate(atoms[:,2]) if val == 2] 
        Os = np.copy(atoms[indsO,:])
        end1 = np.argmin(Os[:,ax_index])
        end2 = np.argmax(Os[:,ax_index])
        length = Os[end2,ax_index]- Os[end1,ax_index]
        
        atomsN = len(atoms)
        nWaters = int(atomsN/atomsMol)
        newLength = (nWaters/density)*10
        scale = newLength/length
        
        
        base = atoms[0,ax_index]
        newAtoms = np.copy(atoms)
        for n in range(atomsN):
            diff = atoms[n,ax_index] - atoms[0,ax_index]
            scaleDiff = scale*diff
            newAtoms[n,ax_index] = scaleDiff + atoms[0,ax_index]
        q.h(7)
        return newAtoms

# => New Constrict <= #            

def scale_density(fn_atoms,
                  axis,    # Along which axis to rescale
                  target_density, # Units must be appropriate for dimensionality
                  axis_scale,
                  ):

	# => Load Info <= #
	
	top = mdt.load(fn_atoms).topology
        atoms = [atom for atom in top.atoms]
        positions = mdt.load(fn_atoms)

        natoms = len(atoms)

        indsO = [ind for ind, val in enumerate(atoms) if val.name == 'O'] 
        indsH = [ind for ind, val in enumerate(atoms) if val.name == 'H'] 
	n_O = len(indsO)
	n_H = len(indsH)

        pairs = []
        for O in indsO:
            for H in indsH:
                pairs.append((H,O))
        
        distances = mdt.compute_distances(positions,pairs,opt=True).reshape((-1,n_O,n_H))

        # => Developing scaling factors <= #
        if target_density != 0.0:
            len_axis = []
            for ind,ax in enumerate(axis):
                length = np.max(positions.xyz[0,:,ax]) -  np.min(positions.xyz[0,:,ax])
                len_axis.append(length)
            
            space = np.prod(np.array(len_axis))
            init_density = space/n_O

            scale = target_density/init_density
            axis_scale = scale**(1.0/float(len(len_axis)))


        # Updating positions 
        new_positions = np.copy(positions.xyz)
        for ind,ax in enumerate(axis):
            base = 0.0
            for n in range(n_O):

                O = indsO[n]
            
                H1 = indsH[np.argsort(distances[0,n,:])[0]]
                H2 = indsH[np.argsort(distances[0,n,:])[1]] 
	        Hneighbors  = np.argmin(distances, axis=2)

                diff_O =  positions.xyz[0,O, ax] - base
                diff_H1 = positions.xyz[0,H1,ax] - positions.xyz[0,O,ax]
                diff_H2 = positions.xyz[0,H2,ax] - positions.xyz[0,O,ax]
                disp_O = axis_scale*diff_O 
                disp_H1= disp_O + diff_H1
                disp_H2= disp_O + diff_H2


                new_positions[0,O,ax] =  disp_O  + base
                new_positions[0,H1,ax] = base + disp_H1
                new_positions[0,H2,ax] = base + disp_H2


        q.h(9)

        return new_positions
        
       
#<==============Rotate molecules================>#
def rotate(atoms, theta, atomsMol=3, interval=2, nMols=0):
	atomsN = nMols*atomsMol
	if nMols!= 0:
		atomsN = len(atoms)
	v = []
	ID = [0]*atomsN
	for n in range(0,atomsN,atomsMol):
		for i in range(n, n+atomsMol,):
			ID[i]=(atomIDs[atoms[i,2]])
			v.append([
				ID[i],
				atoms[i,4],
				atoms[i,5],
				atoms[i,6]
				])
	R = np.matrix([[math.cos(theta), -math.sin(theta)],[math.sin(theta), math.cos(theta)]])
	for k in range(0,atomsN,atomsMol*interval):
		molecule = v[k:k+atomsMol]
		comZ = CM.comZ(molecule)
		comY = CM.comY(molecule)

		for i in range(atomsMol):
			molecule[i][3] -= comZ
			molecule[i][2] -= comY	
			zy = np.matrix(molecule[i][2:4])
			newZY = np.array(zy*R)
			atoms[i+k,6] = newZY[0,0] + comZ
			atoms[i+k,5] = newZY[0,1] + comY
	return atoms

#<==============Format for COM and like modules================>#
def reformat(atoms):
	atomsN = len(atoms)
	v = []
	ID = [0]*atomsN
	for i in range(atomsN):
		ID[i]=(atomIDs[atoms[i,2]])
		v.append([
			ID[i],
			atoms[i,4],
			atoms[i,5],
			atoms[i,6]
			])
	return v

#<==============Re-center in Box (z)================>#
def centerZ(atoms, zhi):
	atomsN = len(atoms)
	
	valsZ = [0]*atomsN
	for i in range(atomsN):
		valsZ[i] = atoms[i,6]
	minZ = np.amin(valsZ)
	maxZ = np.amax(valsZ)
        middleWire = (maxZ + minZ)/2
        middleBox = zhi/2
        disp = middleBox-middleWire
	print disp
	
	for n in range(atomsN):
		atoms[n,6] += disp
	
	return atoms

#<==============Re-assign atom charges================>#
def charges(atoms, atomID, charge):
	atomsN = len(atoms)
	for i in range(atomsN):
		if atoms[i,2] == atomID:
			atoms[i,3] = float(charge)
	return atoms		
		
#<==============Re-assign atom indices================>#
def index(atoms):
	atomsN = len(atoms)
	for i in range(atomsN):
		atoms[i,0] = i+1
	return atoms

#<==============Add Proton Defect================>#
def proton(atoms, wateratomID, h_up=True):
	atomsN = len(atoms)
	newatoms = np.zeros((atomsN+1,10))
	if h_up==True:
		Ox = atoms[wateratomID-1, 4] #subtract one since python indexing starts at 0, lammps starts at 1
		Oy = atoms[wateratomID-1, 5]
		Oz = atoms[wateratomID-1, 6]

		Hx1 = atoms[wateratomID-2, 4] #the h index precedes the O index
		Hy1 = atoms[wateratomID-2, 5]
		Hz1 = atoms[wateratomID-2, 6]
		vH1 = np.array([Hx1-Ox, Hy1-Oy, Hz1-Oz])
		nH1 = LA.norm(vH1)
	
		Hx2 = atoms[wateratomID-3, 4]
		Hy2 = atoms[wateratomID-3, 5]
		Hz2 = atoms[wateratomID-3, 6]
		vH2 = np.array([Hx2-Ox, Hy2-Oy, Hz2-Oz])
		nH2 = LA.norm(vH2)
	
	if h_up==False:
		Ox = atoms[wateratomID-1, 4]
		Oy = atoms[wateratomID-1, 5]
		Oz = atoms[wateratomID-1, 6]

		Hx1 = atoms[wateratomID, 4]
		Hy1 = atoms[wateratomID, 5]
		Hz1 = atoms[wateratomID, 6]
		vH1 = np.array([Hx1-Ox, Hy1-Oy, Hz1-Oz])
		nH1 = LA.norm(vH1)
	
		Hx2 = atoms[wateratomID+1,4]
		Hy2 = atoms[wateratomID+1,5]
		Hz2 = atoms[wateratomID+1,6]
		vH2 = np.array([Hx2-Ox, Hy2-Oy, Hz2-Oz])
		nH2 = LA.norm(vH2)

	HHd = math.sqrt((Hx1-Hx2)**2 + (Hy1-Hy2)**2 + (Hz1-Hz2)**2)
        b = np.array(vH1+vH2)
	b[:] = [x/2 for x in b]	
	c = vH2 - vH1
	
	
	m = [vH1, vH2, c]
	v = np.diag(m)

	c1 = ((HHd**2) - 2*nH1)/(-2)
	c2 = ((HHd**2) - 2*nH2)/(-2)
	c3 = np.dot(b,c)
		
	Hx3 = c1/v[0]
	Hy3 = c2/v[1]
	Hz3 = c3/v[2]

	H3 = np.array([Hx3, Hy3, Hz3])
	nH3 = LA.norm(H3)
	length = ((nH1+nH2)/2)/nH3
	H3[:] = [x*length for x in H3]

	H3[0] += Ox
	H3[1] += Oy
	H3[2] += Oz

	for i in range(atomsN):
		newatoms[i,:] += atoms[i,:]
#		newatoms[atomsN,:] = atoms[i,:]
		if atoms[i,2] == 1:
			newatoms[atomsN,:] = atoms[i,:]
	newatoms[atomsN, 4] = H3[0]
	newatoms[atomsN, 5] = H3[1]
	newatoms[atomsN, 6] = H3[2]
	
	newatoms = index(newatoms)

	return newatoms
		

# => Calculates pair-wise distances between all atoms <= #
def distances(atoms):

        xi, xj = np.meshgrid(atoms[:,4],atoms[:,4],indexing='ij')
        yi, yj = np.meshgrid(atoms[:,5],atoms[:,5],indexing='ij')
        zi, zj = np.meshgrid(atoms[:,6],atoms[:,6],indexing='ij')
        dx = xi - xj
        dy = yi - yj
        dz = zi - zj
        dr = np.sqrt(np.square(dx) + np.square(dy) + np.square(dz))
        
        return dr


# => Cuts Micelle Sphere out of Bulk based on Radius <= #
def r_sphere(atoms, radius, centerX=0.0, centerY=0.0, centerZ=0.0, atomsMol=3):
	atomsN = len(atoms)
	v = reformat(atoms)

	if centerX == 0.0:
		centerX += CM.comX(v)
	if centerY == 0.0:
		centerY += CM.comY(v)
	if centerZ == 0.0:
		centerZ += CM.comZ(v)

	newatoms = []	
	for atom in atoms:
            r = math.sqrt((atom[4]-centerX)**2 + (atom[5]-centerY)**2 + (atom[6]-centerZ)**2) 
            if r < radius:
                newatoms.append(atom)

        newatoms = np.array(newatoms)
        dr = distances(newatoms)
        cutoff = 1.05
        new2atoms = []
        diag_test = np.eye(len(dr))*cutoff*2.0
        dr += diag_test

        for k in range(len(newatoms)):
            if newatoms[k,2] == 2.0:
                if sum(float(num) < cutoff for num in dr[:,k]) < 2:
                    pass
                else:
                    new2atoms.append(newatoms[k])
            else:
                new2atoms.append(newatoms[k])

        new2atoms = np.array(new2atoms)
        finalAtoms = []
        dr2 = distances(new2atoms)
        diag_test2 = np.eye(len(dr2))*cutoff*2.0
        dr2 += diag_test2
        for k in range(len(new2atoms)):
            if sum(float(num) < cutoff for num in dr2[:,k]) < 1:
                pass
            else:
                finalAtoms.append(new2atoms[k])
        
	finalAtoms = index(np.array(finalAtoms))
        if len(finalAtoms) == 0.0:
            print "WARNING: No atoms survived! Modify cutoff distance!"
        elif len(finalAtoms)% 3 != 0.0:
            print "WARNING: Incomplete molecule present: %d atoms total!" % (len(finalAtoms))
       
        nOs = (sum(float(num) == 1.0 for num in finalAtoms[:,2]))
        nHs = (sum(float(num) == 2.0 for num in finalAtoms[:,2]))
        if nHs*2 != nOs:
            print "MISMATCH in hydrogens and oxygens!!!!"
            print "Number of hydrogens: %d" % (nOs)
            print "Number of oxygens: %d" % (nHs)

	return finalAtoms

# => Cuts Micelle Sphere out of Bulk based on Number of Water molecules <= #
def n_sphere(fn_pdb, nwaters, centerX=0.0, centerY=0.0, centerZ=0.0, atomsMol=3):

        frame = mdt.load(fn_pdb)
	top = mdt.load(fn_pdb).topology
        old_atoms = readPDB(fn_pdb)

        atoms = [atom for atom in top.atoms]
	atomsN = len(atoms)

        indsO = [ind for ind,val in enumerate(atoms) if val.name == 'O']
        indsH = [ind for ind,val in enumerate(atoms) if val.name == 'H']
        n_O = len(indsO)
        n_H = len(indsH)
    
	pairs = []
	for indH in indsH:
	    for indO in indsO:
	        pairs.append((indH,indO))

        distances = mdt.compute_distances(frame, pairs,opt=True).reshape((-1,n_H,n_O))

	newatoms = []	
        radii = []
	for ind,O in enumerate(indsO):
            atom = frame.xyz[0,O]
            r = math.sqrt((atom[0]-centerX)**2 + (atom[1]-centerY)**2 + (atom[2]-centerZ)**2) 
            radii.append(r)

        radii = np.array(radii)
        indsO = np.array(indsO)
        indsH = np.array(indsH)

        new_indsO = indsO[np.argsort(radii)[:nwaters]]
        local_indsO = np.argsort(radii)[:nwaters]
        new_indsH = indsH[np.argsort(distances[0,:,local_indsO])[:,:2]].reshape(-1)
       
        old_atoms =  
        H_atoms = frame.xyz[0,new_indsH]
        O_atoms = frame.xyz[0,new_indsO]

        exit()
        

	return finalAtoms
#<==============Reprint data.water with new coordinates================>#
def dataFile(atoms, fn_ref, fn_out="data.water", maxZ=0.0, distance=0.0, xhi=0.0, yhi=0.0, zhi=0.0,endData=False):
	fh = open(fn_ref, "r")
	lines = fh.readlines()
	
	l1 = 0
	index1 = 0
	positions1 = []
	for line in lines:
		l1 += 1
		natoms1 = re.search('^(\S+)\s+atoms\s*$', line)
		if natoms1 != None:
			index1 = l1

	box = []
	l2 = 0
	t1 = 0
	index2 = 0
	for line in lines:
		l2 += 1
		box1 = re.match('^\s*(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s*$', line)
		if box1 != None:
			t1 += 1
			box.append([
				float(box1.group(1)),
				float(box1.group(2)),
				box1.group(3),
				box1.group(4),
				])
			if t1 ==3:
				index2 = l2
	if maxZ != 0.0:
		box[2][1] = maxZ + (distance/2)
	if zhi != 0.0:
		box[2][1] = zhi
	
	for i in range(3):
		box[i][0] = 0.0

	if xhi != 0:
		box[0][1] = xhi
	if yhi != 0:
		box[1][1] = yhi

	indexN = 0
	l3 = 0	
	index3 = l3
	for line in lines:
		l3 += 1
		atoms1= re.match('^\s*(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s*$', line)
		if atoms1 != None:
			index3 = l3
			indexN += 1

	atomsN = len(atoms)

	fj = open(fn_out, "w")
	
	for i in range(index1-1):
		fj.write("%s" % lines[i])
	fj.write("%d %s\n" % (atomsN, "atoms"))
	for i in range(index1,index2-3):
		fj.write("%s" % lines[i])
	fj.write("%.16e %.16e xlo xhi\n" % (box[0][0], box[0][1]))
	fj.write("%.16e %.16e ylo yhi\n" % (box[1][0], box[1][1]))
	fj.write("%.16e %.16e zlo zhi\n" % (box[2][0], box[2][1]))

	for i in range(index2, index3-indexN):
		fj.write("%s" % lines[i])

	for i in range(atomsN):
		fj.write("%d %1d %1d %.16e %.16e %.16e %.16e %1d %1d %1d\n" % (atoms[i][0], atoms[i][1], atoms[i][2], atoms[i][3], atoms[i][4], atoms[i][5], atoms[i][6], atoms[i][7], atoms[i][8], atoms[i][9]))
	
	if endData==True:
		stop = len(lines)
		for i in range(index3,stop):
			fj.write("%s" % lines[i])	

#<==============Reprint data.water with new coordinates================>#
def pdbFile(atoms,filename, maxZ=0.0, distance=0.0, xhi=0, yhi=0, zhi=0):
	fh = open(filename, "r")
	lines = fh.readlines()

	atomsN = len(atoms)
	ID =[0]*atomsN
	for i in range(atomsN):
		ID[i]=(atomIDs[atoms[i,2]])

	box = []
	for line in lines:
		box1 = re.match('^\s*(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s*$', line)
		if box1 != None:
			box.append([
				float(box1.group(1)),
				float(box1.group(2)),
				box1.group(3),
				box1.group(4),
				])
	if maxZ != 0.0:
		box[2][1] = maxZ + (distance/2)
	if zhi != 0.0:
		box[2][1] = zhi

	for i in range(3):
		box[i][0] = 0.0

	if xhi != 0:
		box[0][1] = xhi
	if yhi != 0:
		box[1][1] = yhi
	fj = open("water2.pdb", "w")
	fj.write("CRYST1   %.2f   %.2f   %.3f   %.2f  %.2f  %.2f  1           1\n" % (box[0][1], box[1][1], box[2][1], boxAngle, boxAngle, boxAngle))
	
	for i in range(atomsN):
		fj.write("ATOM  %5d  %s   MOL X   1      %6.2f  %6.2f  %6.2f  1.00  0.00      %s\n" % (atoms[i][0], ID[i], atoms[i][4], atoms[i][5], atoms[i][6], ID[i]))
	fj.write("END")
	
	q.h(13)



def write_LAMMPS(
                 atoms,
                 in_opts = {},
                 ):

    options = {

        # File names 
        "fn_out" : "data.water",

        # Cell parameters
        "xhi" : 0.0,
        "yhi" : 0.0,
        "zhi" : 0.0,  
        "xlo" : 0.0,
        "ylo" : 0.0,
        "zlo" : 0.0,  

        # Atom parameters   
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

        # Symbols
        "symbols" : ['H','O'],

    }

    # => Override Default Options <= #

    for key,val in in_opts.items():
       if key not in options.keys():
          raise ValueError('%s option unavailable, RTFM' % (key))
       if options[key] != val:
          print "Default Option Overridden: %s = %s" % (key,str(val))
       options[key] = val
    print '\n'

    # => Declaring Variables from Options <= #

    fn_out = options['fn_out'] 

    xhi = options['xhi'] 
    yhi = options['yhi']
    zhi = options['zhi']
    xlo = options['xlo']
    ylo = options['ylo']
    zlo = options['zlo']
    
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

    natoms = len(atoms)

    if xlo == 0.0 and np.min(atoms[:,4]) < 0.0:
        xlo = np.min(atoms[:,4])
    if ylo == 0.0 and np.min(atoms[:,5]) < 0.0:
        ylo = np.min(atoms[:,5])
    if zlo == 0.0 and np.min(atoms[:,6]) < 0.0:
        zlo = np.min(atoms[:,6])

    if xhi == 0.0:
        xhi = np.max(atoms[:,4])
    if yhi == 0.0:
        yhi = np.max(atoms[:,5])
    if zhi == 0.0:
        zhi = np.max(atoms[:,6])
        
    timestamp = "LAMMPS data file. Written, date: %s\n\n" % (str(date.today()))
    hs = open(fn_out,"w")
    hs.write(timestamp)
    hs.write(" %d atoms\n" % (natoms))
    hs.write(" %d bonds\n" % (bonds))
    hs.write(" %d angles\n" % (angles))
    hs.write(" %d dihedrals\n" % (dihedrals))
    hs.write(" %d impropers\n" % (impropers))
    hs.write(" %d atom types\n" % (atom_types))
    hs.write(" %d bond types\n" % (bond_types))
    hs.write(" %d angle types\n" % (angle_types))
    hs.write(" %d dihedral types\n" % (dihedral_types))
    hs.write(" %d improper types\n" % (improper_types))
    hs.write(" %4.4f %4.4f xlo xhi\n" % (xlo,xhi))
    hs.write(" %4.4f %4.4f ylo yhi\n" % (ylo,yhi))
    hs.write(" %4.4f %4.4f zlo zhi\n\n" % (zlo,zhi))
    hs.write(" Masses\n\n")
   
    for ind,mass in enumerate(masses):
        hs.write(" %d %2.6f # %s\n" % (ind+1,mass,symbols[ind]))
    hs.write("\n\n")
    hs.write(" Atoms\n\n")
       
   
    for i in range(natoms):
        hs.write("%d %1d %1d %.16e %.16e %.16e %.16e %1d %1d %1d\n" % (atoms[i][0], atoms[i][1], atoms[i][2], atoms[i][3], atoms[i][4], atoms[i][5], atoms[i][6], atoms[i][7], atoms[i][8], atoms[i][9]))
    hs.write("\n")

    hs.close()
