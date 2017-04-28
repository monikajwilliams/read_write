#!/usr/bin/env python
import os, sys, re, math
import numpy as np


# NOTE: Following functions do NOT require this as the index for the positions
# in plumed start at unity, so only the number of atoms is required.
# => "r001: POSITION ATOM=  1"<= #
def position(
    targetAtoms,
    name='r'
    ):

    fh = open("plumed.dat", "w")
    natoms = len(targetAtoms)
    for i in range(natoms):
        fh.write("%s%05d: POSITION ATOM=%d NOPBC" % (name,i + 1, targetAtoms[i] ));
    	#print "%s%05d: POSITION ATOM=%d NOPBC" % (name,i + 1, targetAtoms[i] )
    	fh.write('\n')
    fh.write('\n')
    fh.close()
	
# calculates the distance of the specified atoms from the  z axis at the
# center specified by the user. 	
# => "rho001: MATHEVAL ARG=r001.x, r001.y FUNC=sqrt(x^2+y^2) PERIODIC=NO"<= #
def rho(
    natoms, 
    centerX=1.6, 
    centerY=1.6,
    arg='r',
    name='rho'
    ):

    fh = open("plumed.dat", "a")
    for i in range(natoms):
    	fh.write("%s%05d: MATHEVAL ARG=%s%05d.x,%s%05d.y FUNC=sqrt((x-%.3f)^2+(y-%.3f)^2) PERIODIC=NO" % (name,i+1,arg, i+1,arg, i+1,centerX,centerY));
    	fh.write('\n')
    fh.write('\n')

# calculates the 3D distance of the specified atoms from the
# center specified by the user. 	
# => "rad001: MATHEVAL ARG=r001.x, r001.y, r001.z, FUNC=sqrt(x^2+y^2+z^2) PERIODIC=NO"<= #
def rad(
    natoms, 
    centerX=1.6, 
    centerY=1.6,
    centerZ=1.6,
    name='rad',
    arg='r'
    ):

    fh = open("plumed.dat", "a")
    for i in range(natoms):
    	fh.write("%s%05d: MATHEVAL ARG=%s%05d.x,%s%05d.y,%s%05d.z FUNC=sqrt((x-%.3f)^2+(y-%.3f)^2+(z-%.3f)^2) PERIODIC=NO" % (name,i+1,arg,i+1,arg,i+1,arg,i+1,centerX,centerY,centerZ));
    	fh.write('\n')
    fh.write('\n')
	
# Applies an upper limit to the radial distance calculated in "rho"
# => "uwall: UPPER_WALLS ARG=rho001, rho002... AT=5.00 KAPPA=200 EXP=2"<= #
def uwall(
    natoms, 
    ARG="rho",
    at=0.1,
    kappa=1000.0, 
    exp=2.0, 
    name="uwall"
    ):

    fh = open("plumed.dat", "a")
    nmol = natoms

    fh.write("%s: UPPER_WALLS ARG=" % (name));
    for i in range(nmol-1):
    	fh.write("%s%05d," %(ARG,i+1));
    fh.write("%s%05d " % (ARG,nmol));

    fh.write("AT=");
    at_list = nmol*['{:.3f}'.format(at)]
    fh.write(','.join(at_list));
    fh.write(" KAPPA=");
    kappa_list = nmol*['{:.1f}'.format(kappa)]
    fh.write(','.join(kappa_list));
    fh.write(" EXP=");
    exp_list = nmol*['{:.1f}'.format(exp)]
    fh.write(','.join(exp_list));
    fh.write('\n')
    fh.write('\n')

# Applies an lower limit to the radial distance calculated in "rho"
# => "lwall: LOWER_WALLS ARG=rho001, rho002... AT=5.00 KAPPA=200 EXP=2"<= #
def lwall(
    natoms, 
    ARG="rho", 
    at=0.1,
    kappa=1000.0, 
    exp=2.0, 
    name="lwall"
    ):

    fh = open("plumed.dat", "a")
    nmol = natoms
    fh.write("%s: LOWER_WALLS ARG=" % (name));
    for i in range(nmol-1):
    	fh.write("%s%05d," %(ARG,i+1));
    fh.write("%s%05d " % (ARG,nmol));
    fh.write("AT=");
    at_list = nmol*['{:.3f}'.format(at)]
    fh.write(','.join(at_list));
    fh.write(" KAPPA=");
    kappa_list = nmol*['{:.1f}'.format(kappa)]
    fh.write(','.join(kappa_list));
    fh.write(" EXP=");
    exp_list = nmol*['{:.1f}'.format(exp)]
    fh.write(','.join(exp_list));
    fh.write('\n')
    fh.write('\n')

# => Applies general wall <= #
def wall(
    natoms, 
    at=0.1,
    kappa=1000.0, 
    exp=2.0, 
    name="uwally",
    d = 'y',
    upper=True,
    ):

    fh = open("plumed.dat", "a")
    nmol = natoms

    if upper == True:
        fh.write("%s: UPPER_WALLS ARG=" % (name));
    else:
        fh.write("%s: LOWER_WALLS ARG=" % (name));

    for i in range(nmol-1):
    	fh.write("r%05d.%s," %(i+1),d);
    fh.write("r%05d.%s " % (nmol),d);
    fh.write("AT=");
    at_list = nmol*['{:.3f}'.format(at)]
    fh.write(','.join(at_list));
    fh.write(" KAPPA=");
    kappa_list = nmol*['{:.1f}'.format(kappa)]
    fh.write(','.join(kappa_list));
    fh.write(" EXP=");
    exp_list = nmol*['{:.1f}'.format(exp)]
    fh.write(','.join(exp_list));
    fh.write('\n')
    fh.write('\n')

#Applies a custom potential using BIASVALUE
# => "Writes Custom Potential"<= #
def custom(
    natoms,
    arg='rho',
    ):

    fh = open("plumed.dat", "a")
    nmol = natoms
    
    cal2J = 4.184
    
    a2i = -0.2281 #kcal/(mol*A^2)
    a4i = 1.09    #kcal/(mol*A^4)
    a6i = 0.2341  #kcal/(mol*A^6)
    a8i = 0.3254  #kcal/(mol*A^8)
    
    a2 = a2i*cal2J*(10**2) #kJ/(mol*nm^2)
    a4 = a4i*cal2J*(10**4) #kJ/(mol*nm^4)
    a6 = a6i*cal2J*(10**6) #kJ/(mol*nm^6)
    a8 = a8i*cal2J*(10**8) #kJ/(mol*nm^8)
    
    for i in range(nmol):
    	fh.write("P%05d: MATHEVAL ARG=%s%05d FUNC=(%.3f*(x^2))+(%.3f*(x^4))+(%.3f*(x^6))+(%.3f*(x^8)) PERIODIC=NO" % (i+1,arg,i+1,a2,a4,a6,a8));
    	fh.write('\n')
    fh.write('\n')
    
    fh.write("custB: BIASVALUE ARG=");
    for i in range(nmol-1):
    	fh.write("P%05d," %(i+1));
    fh.write("P%05d " %(nmol));
    fh.write('\n')
    fh.write('\n')
    
#Another specific potential
# => "Writes Voth's Potential for a Given Radius"<= #
def Voth(
    natoms,  
    r_cyl=0.225,
    ):

    fh = open("plumed.dat", "a")
    nmol = natoms
    A = 38.72 #kcal/mol
    B = 4.18 #\AA^{-1}
    
    for i in range(nmol):
    	fh.write("P%05d: MATHEVAL ARG=rho%05d FUNC=(%.3f*exp(-%.3f*(%.3f-x))) PERIODIC=NO" % (i+1,i+1,A,B,r_cyl));
    	fh.write('\n')
    fh.write('\n')
    
    fh.write("Voth: BIASVALUE ARG=");
    for i in range(nmol-1):
    	fh.write("P%05d," %(i+1));
    fh.write("P%05d " %(nmol));
    fh.write('\n')
    fh.write('\n')
    
def philic(
    natoms,
    alpha,
    beta,
    sigma=0.3,
    epsilon=29.0,
    radius = 0.6,
    name_1="phil",
    name_2="philic",
    arg="rho",
    ):

    fh = open("plumed.dat", "a")
    eps_ls = np.sqrt(epsilon)
    nmol = natoms

    for i in range(natoms):
    	fh.write("%s%05d: MATHEVAL ARG=%s%05d FUNC=4.0*%.3f*%.3f*((%.3f/(%.3f-x))^12.0-%.3f*((%.3f/(%.3f-x))^6.0)) PERIODIC=NO" % (name_1,i+1,arg,i+1,alpha,eps_ls,sigma,radius,beta,sigma,radius));
    	fh.write('\n')
    fh.write('\n')
    
    fh.write("%s: BIASVALUE ARG="%(name_2));
    for i in range(nmol-1):
    	fh.write("%s%05d," %(name_1,i+1));
    fh.write("%s%05d " %(name_1,nmol));
    fh.write('\n')
    fh.write('\n')

#Applies a custom potential using BIASVALUE
# => "Writes New Spherical Potential"<= #
def buckingham(
    natoms,
    Rmax,
    sigma,
    epsilon,
    alpha,
    arg='rad',
    name='bias',
    ):

    fh = open("plumed.dat", "a")
    nmol = natoms

    A = epsilon*(6.0/(alpha-6.0))
    B = Rmax + sigma
   
    for i in range(nmol):
    	fh.write("P%05d: MATHEVAL ARG=%s%05d FUNC=(%.3f*(exp(%.3f*(1.0-((%.3f-x)/%.3f))))) PERIODIC=NO" % (i+1,arg,i+1,A,alpha,B,sigma));
    	fh.write('\n')
    fh.write('\n')
    
    fh.write("%s: BIASVALUE ARG=" % (name));
    for i in range(nmol-1):
    	fh.write("P%05d," %(i+1));
    fh.write("P%05d " %(nmol));
    fh.write('\n')
    fh.write('\n')

#Printing variables
# => "PRINT ARG= FILE="<= #
def output(
    nVals,
    ARG, 
    filename, 
    bias=False
    ):	

    fh = open("plumed.dat", "a")
    fh.write("PRINT  ARG=");

    if nVals > 1:
        if bias == True:
            for i in range(nVals-1):
                fh.write("%s%05d_bias," %(ARG,i+1));
	    fh.write("%s%05d_bias " %(ARG,nVals));
	if bias == False:
	    for i in range(nVals-1):
	        fh.write("%s%05d," %(ARG,i+1));
	    fh.write("%s%05d " %(ARG,nVals));
	if nVals == 1:
	    fh.write("%s " % (ARG));
	fh.write("FILE=%s" % (filename));
	fh.write('\n')
    fh.write('\n')

# => Combinations of Segments, Application Scope <= #	

def quadratic(
    filename,
    at, 
    kappa, 
    exp, 
    uzat, 
    uzkappa, 
    lzat, 
    lzkappa,
    zexp=2.0, 
    targetID=2,
    ):

    atomIDs = targetAtoms(filename, targetID)
    natoms = len(atomIDs)
    position(atomIDs)
    rho(natoms)
    uwall(natoms, at=at, kappa=kappa, exp=exp)
    uwall(natoms, at=uzat, kappa=uzkappa,exp=zexp, name="uwallz")	
    lwallz(natoms, molAtom, lzat, lzkappa, lzexp=zexp)
    #output(nVals=natoms, ARG="rho", filename="rhos.txt")
    #output(nVals=1, ARG="uwall.bias", filename="bias.txt")

def customBias(
    filename,
    uzat, 
    uzkappa, 
    uzexp,
    lzat, 
    lzkappa, 
    lzexp, 
    targetID=2,
    ):

    atomIDs = targetAtoms(filename, targetID)
    natoms = len(atomIDs)
    position(atomIDs)
    rho(natoms)
    custom(natoms)
    uwallz(natoms,uzat=uzat, uzkappa=uzkappa, uzexp=uzexp, name="uwallz")	
    lwallz(natoms,molAtom, lzat, lzkappa, lzexp=lzexp)
    #output(nVals=natoms, ARG="rho", filename="rhos.txt")
    #output(nVals=natoms, ARG="custB.P", filename="bias.txt", bias=True)

def litBias(
    filename,
    xhi, 
    yhi, 
    targetID=2,
    ):

    centerX = xhi/2
    centerY = yhi/2
    atomIDs = targetAtoms(filename, targetID)
    natoms = len(atomIDs)
    position(atomIDs)
    rho(natoms, centerX=centerX, centerY=centerY)
    custom(natoms)
    #output(nVals=natoms, ARG="rho", filename="rhos.txt")
    #output(nVals=natoms, ARG="custB.P", filename="bias.txt", bias=True)

def sphere(
    atomIDs,
    centerX, 
    centerY,
    centerZ,
    at,
    kappa, 
    exp, 
    targetID=2
    ):

    #atomIDs = targetAtoms(filename, targetID)
    natoms = len(atomIDs)
    position(atomIDs)
    rad(natoms,centerX,centerY,centerZ)
    uwall(natoms,ARG="rad",at=at,kappa=kappa, exp=exp,name="uwall")
    #output(nVals=natoms, ARG="rad", filename="rad.txt")
    #output(nVals=1, ARG="uwall.bias", filename="bias.txt", bias=True)

# This was a crazy idea which does not work yet. 
def ring(
    filename, 
    rad1, 
    rad2,
    xhi,
    yhi, 
    zhi, 
    kappa, 
    exp, 
    uzat,
    lzat,
    targetID=2
    ):

    centerZ = zhi/2
    centerX = xhi/2
    centerY = yhi/2
    atomsIDs = targetAtoms(filename, targetID)
    natoms = len(atomsIDs)
    position(atomsIDs)
    rho(natoms, centerX=centerX, centerY=centerY)
    uwallz(natoms,uzat=uzat, uzkappa=kappa, uzexp=exp, name="uwallz")	
    lwallz(natoms,molAtom, lzat, lzkappa=kappa, lzexp=exp)
    uwall(natoms,ARG="rho",at=rad1,kappa=kappa, exp=exp,name="uwall")
    lwall(natoms,ARG="rho",at=rad2,kappa=kappa, exp=exp,name="lwall")


def vothBias(
    filename,
    xhi, 
    yhi, 
    targetID=2,
    r_cyl=0.225,
    ):

    centerX = xhi/2
    centerY = yhi/2
    atomIDs = targetAtoms(filename, targetID)
    natoms = len(atomIDs)
    position(atomIDs)
    rho(natoms, centerX=centerX, centerY=centerY)
    Voth(natoms,r_cyl=r_cyl)

def vothBias2(
    filename,
    xhi, 
    yhi,
    uzat,
    uzkappa,
    uzexp,
    lzat,
    lzkappa,
    lzexp, 
    targetID=2,
    r_cyl=0.225,
    ):

    centerX = xhi/2
    centerY = yhi/2
    atomIDs = targetAtoms(filename, targetID)
    natoms = len(atomIDs)
    position(atomIDs)
    rho(natoms, centerX=centerX, centerY=centerY)
    Voth(natoms,r_cyl=r_cyl)
    uwallz(natoms,uzat=uzat, uzkappa=uzkappa, uzexp=uzexp, name="uwallz")	
    lwallz(natoms,molAtom, lzat, lzkappa, lzexp=lzexp)

def plane(atomIDs,
   xlat, 
   xuat, 
   ylat = 0, 
   yuat = 0, 
   zlat = 0, 
   zuat = 0, 
   kappa=200.0, 
   exp=2,
   targetID=2,
   ):
#	atomIDs = targetAtoms(filename, targetID)
   natoms = len(atomIDs)
   position(atomIDs)
   
   uwallx(natoms,at=xuat,kappa=kappa, exp=exp,name="uwallx")
   lwallx(natoms,at=xlat,kappa=kappa, exp=exp,name="lwallx")

   if ylat != 0 and yuat != 0:
      uwally(natoms,at=yuat,kappa=kappa, exp=exp,name="uwally")
      lwally(natoms,at=ylat,kappa=kappa, exp=exp,name="lwally")

   if zlat != 0 and zuat != 0:
      uwallz(natoms,uzat=zuat,uzkappa=kappa, uzexp=exp,name="uwallz")
      lwallz(natoms,lzat=zlat,lzkappa=kappa, lzexp=exp,name="lwallz")

def hydrophilic(
    O_IDs,
    H_IDs,
    xhi,
    yhi,
    alpha,
    beta,
    sigma=0.3,
    epsilon=29.0,
    radius = 0.6,
    ):

   centerX = xhi/2
   centerY = yhi/2
   nOs = len(O_IDs)
   nHs = len(H_IDs)

   position(O_IDs,name="r_O")
   rho(natoms=nOs, centerX=centerX, centerY=centerY,name="rho_O",arg="r_O")
   custom(natoms=nOs,arg="rho_O")

   position(H_IDs,name="r_H")
   rho(natoms=nHs, centerX=centerX, centerY=centerY,name="rho_H",arg="r_H")
   philic(natoms=nHs,
      alpha=alpha,
      beta=beta,
      sigma=sigma,
      epsilon=epsilon,
      name="phil",
      arg="rho_H",
      radius = radius,
      )

def hydro_grid(O_IDs,
      H_IDs,
      point_coords,
      xhi,
      yhi,
      alpha,
      beta,
      sigma=0.3,
      epsilon=29.0,
      radius = 0.6,
      ):

   centerX = xhi/2
   centerY = yhi/2
   nOs = len(O_IDs)
   nHs = len(H_IDs)

   position(O_IDs,name="r_O")
   rho(natoms=nOs, centerX=centerX, centerY=centerY,name="rho_O",arg="r_O")
   custom(natoms=nOs,arg="rho_O")

   position(H_IDs,name="r_H")
   for ind,point in enumerate(point_coords):
       rad(natoms=nHs, centerX=point[0], centerY=point[1],centerZ=point[2],name="rad_H%05d" % (ind),arg="r_H")
       philic(natoms=nHs,
     alpha=alpha,
     beta=beta,
     sigma=sigma,
     epsilon=epsilon,
     name_1="phil%05d" % (ind),
     name_2="philic%05d" % (ind),
     arg="rad_H%05d" % (ind),
     radius = radius,
     )

def micelle(
    atomIDs,
    centerX, 
    centerY,
    centerZ,
    Rmax,
    sigma,
    epsilon,
    alpha,
    arg='rad',
    name='bias',
    ):

    natoms = len(atomIDs)
    position(atomIDs)
    rad(natoms,centerX,centerY,centerZ)
    buckingham(
        natoms,
        Rmax,
        sigma,
        epsilon,
        alpha,
        arg='rad',
        name='bias',
        )









