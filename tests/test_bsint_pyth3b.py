# -*- coding: utf-8 -*-
# 
#  test_bsint_pyth3b.py
#  glbsint
#  
#  Created by Alexander Rudy on 2013-02-12.
#  Copyright 2013 Alexander Rudy. All rights reserved.
# 
from __future__ import division
"""
Testing the 3-body problem
"""

import numpy as np

class test_pyth3body(object):
    """Test the pythagorean 3-body problem"""
    
    
    def test_pyth3body(self):
        """Test the pythagorean 3-body approach"""
        mass = np.array([3.0,4.0,5.0])
        t0 = 0
        t1 = 100
        nb = 3
        y = np.zeros((6*nb))
        
        y[0] = 1 # 1x
        y[1] = 0 # 1vx
        y[2] = 3 # 1y
        y[3] = 0 # 1vy
        y[4] = 0
        y[5] = 0
        
        y[6] = -2 #2x
        y[7] = 0
        y[8] = -1 #2y
        y[9] = 0
        y[10] = 0
        y[11] = 0
        
        y[12] = 1
        y[13] = 0
        y[14] = -1
        y[15] = 0
        y[16] = 0
        y[17] = 0
        
        def force(t,y):
            """Force law!"""
            n = len(y)
            nb = n//6
            denom = np.zeros((n,n))
            dydx = np.zeros(y.shape)
            
            for i in range(1,n,2):
                dydx[i-1] = y[i]
            
            for i in range(nb):
                for j in range(nb):
                    if i != j:
                        ib = (i-1)*6
                        jb = (j-1)*6
                        denom_tmp = (y[jb]-y[ib])**2.0 + (y[jb+2]-y[ib+2])**2.0 + (y[jb+4]-y[ib+4])**2.0
                        denom[i,j] = np.power(denom_tmp,3/2)
                        denom[j,i] = denom[i,j]
            
            for i in range(nb):
                ib = (i-1)*6
                for ic in range(3):
                    dydx[ib+(2*ic)] = 0
                    for j in range(nb):
                        jb = (i-1)*6
                        if i != j:
                            dydx[ib+(2*ic)] = dydx[ib+(2*ic)] - mass[j] * (y[ib+2*ic-1] - y[jb+2*ic-1])/denom[i,j]
                            
            print dydx
            return dydx
            
        import glbsint.bsint as bsint
        
        yout, tout = bsint.bsintegrate(force,y,t0,t1,0,0,1000)
        
        print tout.shape, yout.shape
        np.savetxt('youtput.txt',yout)
        np.savetxt('toutput.txt',tout)
                        
        