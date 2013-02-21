bsint - A Bulirsch-Stoer Integrator
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is a python adaptor for the fortan-based Bulirsch-Stoer integrator. The integrator is written in fortran, from Numerical Recipies. An example implementation is included in the file `pythagorean_threebody.py`. After installing, run this file from the command-line, and examine `output.dat` and `output.png` to see the results.

Install
=======

This module requires Numpy_. If you have Numpy_ installed and accessible, it should work. To install this module, use:

    $ sudo python setup.py install
    
If this works proprely, you will have a ``bsint`` module with a ``bsintegrate`` method that you can use to carry out Bulirsch-Stoer integrations.

Try it out!
===========

To run the example program, try

    $ ./pythagorean_threebody.py
    
and look at ``output.png``!

Use it!
=======

To use the ``bsintegrate()`` function, you'll need to write a function which calculates the derivatives of your system of equations. The signature of your derivatives function should be ``derivs(t,Ys,*args)`` where ``args`` are any extra arguments you wish to pass through the integrator to the derivatives function. Derivatives should return an array ``dYsdt``, the derivatives of ``Ys``, which should be the same shape as ``Ys``. ``Ys`` and ``dYsdt`` should be Numpy_ arrays.


The ``bsintegrate`` function itself:

    yout,tout = bsintegrate(derivs,y,t0,t1,[tacc,h0,mxstep,args])

Parameters
----------
``derivs`` : call-back function. Should take the current timestep and position. derivs(t,y)
``y`` : input array of variables for integration at t0.
``t0`` : float, starting time for integration
``t1`` : float, ending point for integration

Other Parameters
----------------
``args`` : input tuple, optional
    Default: ``()``
``tacc`` : input float, optional
    Default: ``1e-14``
``h0`` : input float, optional
    Default: ``1e-3``
``mxstep`` : input int, optional
    Default: ``1e4``

Returns
-------
``yout`` : rank-2 array() with bounds (steps,nes)
``tout`` : rank-1 array() with bounds (steps)

Notes
-----
Call-back functions::

  def derivs(x,y): return dydx
  Required arguments:
    x : input float
    y : input rank-1 array('d') with bounds (n)
  Return objects:
    dydx : rank-1 array('d') with bounds (n)

.. _Numpy: http://www.numpy.org
