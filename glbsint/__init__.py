# 
#  __init__.py
#  glbsint
#  
#  Created by Alexander Rudy on 2013-02-12.
#  Copyright 2012 Alexander Rudy. All rights reserved.
# 

__all__ = ['bsintegrate']

from bsint import bsintegrate
bsintegrate.__doc__ = """
A Bulirsch-Stoer Integrator. Based on Greg Laughlin's integration in fewbody.f

yout,tout = bsintegrate(derivs,y,t0,t1,[tacc,h0,mxstep,derivs_extra_args])

Wrapper for ``bsintegrate``.

Parameters
----------
derivs : call-back function. Should take the current timestep and position. derivs(t,y)
y : input array of variables for integration at t0.
t0 : float, starting time for integration
t1 : float, ending point for integration

Other Parameters
----------------
derivs_extra_args : input tuple, optional
    Default: ()
tacc : input float, optional
    Default: 1e-14
h0 : input float, optional
    Default: 0.001
mxstep : input int, optional
    Default: 1000

Returns
-------
yout : rank-2 array('d') with bounds (mxstep,nes)
tout : rank-1 array('d') with bounds (mxstep)

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