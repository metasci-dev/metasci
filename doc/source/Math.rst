*******************
Mathematics Package
*******************
Other code packages offer great ways to do some very complicated mathematics.  
However, they sometimes do not have quite what is required for a certain task. 
This package is meant to fill in those gaps.  For instance, NumPy offers a very general 
`one dimensional polynomial class <http://docs.scipy.org/doc/numpy/reference/generated/numpy.poly1d.html>`_.
Unfortunately, they do not (yet) have 2, 3, 4 ... N-dimensional analogies.  MetaSci provides this here.

Additionally, some useful special functions are also offered that do not come stock in other packages.

=========================================
:mod:`mathematics` -- General Mathematics
=========================================

.. automodule:: metasci.mathematics
   :members:

==================================================================
:mod:`mathematics.polynomial` -- A More General Polynomial
==================================================================

.. automodule:: metasci.mathematics.polynomial
   :members:
   :inherited-members:
   :exclude-members: polyNd, polyNdcoef

   -----------------------
   Polynomial Coefficients
   -----------------------
   .. autoclass:: polyNdcoef
      :members:

   --------------------
   Polynomial Functions
   --------------------
   .. autoclass:: polyNd
      :members:

      .. automethod:: __call__


=================================================================
:mod:`mathematics.integrate` -- Some Integration Helper Functions
=================================================================

.. automodule:: metasci.mathematics.integrate
   :members:
   :inherited-members:
