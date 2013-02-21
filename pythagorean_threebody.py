#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 
#  test_bsint.py
#  glbsint
#  
#  Created by Alexander Rudy on 2013-02-12.
#  Copyright 2013 Alexander Rudy. All rights reserved.
# 

from __future__ import division
import numpy as np

# Make NUMPY print lots of things on a single line.
np.set_printoptions(linewidth=200,precision=2)

# Mass of our particles
mass = np.array([3.0,4.0,5.0])

# Integration Time for the 3-body system
t0 = 0
t1 = 65

# Number of bodies in use
nb = 3
Y0 = np.zeros((6*nb))

# Initial Conditions
Y0[0] = 1 # x1
Y0[1] = 0 # vx1
Y0[2] = 3 # y1
Y0[3] = 0 # vy1
Y0[4] = 0 # z1
Y0[5] = 0 # vz1
        
Y0[6] = -2 # x2
Y0[7] = 0  # vx2
Y0[8] = -1 # y2
Y0[9] = 0  # vy3
Y0[10] = 0 # z2
Y0[11] = 0 # vz2
        
Y0[12] = 1 # x3
Y0[13] = 0 # vx3
Y0[14] = -1 # y3
Y0[15] = 0 # vy3
Y0[16] = 0 # z3
Y0[17] = 0 # vz3
        
def force(t,Y,mass):
    """Force law for basic Grvaitation"""
    
    n = len(Y)
    nb = n//6 # Number of bodies, Use floor'd divison!
    
    denom = np.zeros((nb,nb)) # Denominator for gravitation (r-vectors)
    dYdt = np.zeros(Y.shape)  # Ouptut quantities, should match input ones.
            
    # Copy velocities from Y to dYdt
    for i in range(1,n,2):
        dYdt[i-1] = Y[i]
    
    # Calculate Separation Vectors
    for i in range(nb):
        for j in range(nb):
            if i != j:
                ib = i*6
                jb = j*6
                denom_tmp = (Y[jb]-Y[ib])**2.0 + (Y[jb+2]-Y[ib+2])**2.0 + (Y[jb+4]-Y[ib+4])**2.0
                denom[i,j] = np.power(denom_tmp,3/2)
                denom[j,i] = denom[i,j]
    
    # Apply force law
    for i in range(nb):
        ib = i*6
        for ic in range(1,4):
            dYdt[ib+(2*ic)-1] = 0
            for j in range(nb):
                jb = j*6
                if i != j:
                    dYdt[ib+(2*ic)-1] = dYdt[ib+(2*ic)-1] - mass[j] * (Y[ib+2*ic-2] - Y[jb+2*ic-2])/denom[i,j]
    return dYdt
            
import bsint
        
yout, tout = bsint.bsintegrate(force,Y0,t0,t1,mxstep=10000,args=(mass,))

# Plot the results!
import matplotlib.pyplot as plt
plt.plot(yout[:,0],yout[:,2],'.-',label="m=3")
plt.plot(yout[:,6],yout[:,8],'.-',label="m=4")
plt.plot(yout[:,12],yout[:,14],'.-',label="m=5")
plt.axis('equal')
plt.title("Pythagorean 3-body Problem")
plt.legend()
plt.savefig("output.png")
flatten = lambda l: [item for sublist in l for item in sublist]
header = "{:25s}".format("t") + " ".join(map("{:25s}".format,flatten(["{0:d}x {0:d}vx {0:d}y {0:d}vy {0:d}z {0:d}vz".format(i+1).split() for i in range(nb)])))
np.savetxt('output.dat',np.hstack((np.atleast_2d(tout).T,yout)),header=header,fmt="%25.18e")
