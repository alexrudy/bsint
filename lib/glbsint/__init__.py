# 
#  __init__.py
#  glbsint
#  
#  Created by Alexander Rudy on 2013-02-12.
#  Copyright 2012 Alexander Rudy. All rights reserved.
# 

__all__ = ['bsintegrate']

import bsint
import numpy
bsint.bsintegrate.__doc__ = """
A Bulirsch-Stoer Integrator. Based on Greg Laughlin's integration in fewbody.f

yout,tout = bsintegrate(derivs,y,t0,t1,[tacc,h0,mxstep,args])

Wrapper for ``bsintegrate``.

Parameters
----------
derivs : call-back function. Should take the current timestep and position. derivs(t,y)
y : input array of variables for integration at t0.
t0 : float, starting time for integration
t1 : float, ending point for integration

Other Parameters
----------------
args : input tuple, optional
    Default: ()
tacc : input float, optional
    Default: 1e-14
h0 : input float, optional
    Default: 1e-3
mxstep : input int, optional
    Default: 1e4

Returns
-------
yout : rank-2 array('d') with bounds (steps,nes)
tout : rank-1 array('d') with bounds (steps)

Notes
-----
Call-back functions::

  def derivs(x,y): return dydx
  Required arguments:
    x : input float
    y : input rank-1 array('d') with bounds (n)
  Return objects:
    dydx : rank-1 array('d') with bounds (n)

"""

def bsintegrate(derivs,y,t0,t1,tacc=1e-14,h0=1e-3,mxstep=1e4,args=()):
    t, y = bsint.bsintegrate(derivs,y,t0,t1,tacc,h0,mxstep,args)
    lindex = numpy.argmax(t)+1
    return t[:lindex], y[:lindex]
    
bsintegrate.__doc__ = bsint.bsintegrate.__doc__