*********************
Graphing and Plotting
*********************
The `MatPlotLib <http://matplotlib.sourceforge.net/>`_ package is one of the most 
useful tools for scientific computing.  However due to its generality, it is quite verbose.
Often the same commands are repetitively called to make figures that are only slightly different.
When this happens in a project, functions are created to wrap figure generation.
However, when the same or similar graphing commands are called across several projects, a central graphing 
package is warranted.  This is it.

.. warning::
   These functions are still under development and likely do not work!

==============================================================
:mod:`metasci.graph` -- MatPlotLib Helper Functions [unstable]
==============================================================

.. automodule:: metasci.graph
   :members:
   :inherited-members:
