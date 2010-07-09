******************
Contingency Tables
******************
A main part of statistics package is the ability to quickly and painlessly calculate contingency tables.
To date, there are two table classes: 2D and 3D.  These are formulated in a more restricted sense than a pure,
mathematical contingency table.

   * **Independent vs Response:** While they can be used more generally, the tables are meant to test associations
     between an independent variable(s) (*x*, *y*) and a response variable (*R*).  However,  in an information theoretic 
     sense, any two parameters may have a contingency table regardless of any implied dependencies.  Independent/Dependent 
     quantities are implied here because this is the natural way to study large, multivariate systems.
   * **Continuous (Non-Discrete) Variables:**  Pure contingency tables are often formulated for parameters that take one 
     discrete values. For instance, "Animal Diet" may be "Carnivore", "Omnivore", "Herbivore".  However, the *R*, *x*, and *y*
     parameters here are continuous and binned into *I*, *J*, or *K* bins.  This is done because most scientific measures 
     (eg "Age", "Temperature", "Pressure") are continuous parameters.

The contingency tables work via initialization.  When a new object is instantiated, all information is automatically calculated.  
The table is generated from lists (or NumPy arrays) of data that are subsequently binned.  The table is stored as an integer 
number of counts per bin.  These follow the naming convention that everything after an underscore ('_') is a subscript.  
Moreover, a 'dot' in the subscript indicates a summation over the dotted index. The total number of entries in a dataset is 
**N**.  These are stored as attributes of the table object.  Applying the naming conventions in 2D we see that:

   * :math:`N_{ij}` = *ContingencyTable.N_ij*, The contingency table itself.
   * :math:`N_{i\cdot}` = *ContingencyTable.N_idot*, Sum over the jth independent variable index.  Array of length *I*.
   * :math:`N_{\cdot j}` = *ContingencyTable.N_dotj*, Sum over the ith response variable index.  Array of length *J*.
   * :math:`N_{\cdot\cdot}` = :math:`N` = *ContingencyTable.N*, Total number of points in dataset.

Similar definitions exist for the associated probability table.  Simply replace *N* with *p* to access these.

Moreover upon initialization, a large number of statistics about the table are also calculated.  
These are also stored as attributes of the object.  They include, but are not limited to:

   * :math:`\chi^2`
   * Contingency Coefficient C
   * Cramer's V
   * G-test Statistic
   * R-value Statistic
   * Entropy
   * Mutual Information
   * Uncertainties

Entropies and mutual information have a separate naming convention from the table itself.  Everything after 
the first underscore is a list of the variables that the entropy is a function of.  The second underscore is the 
delimiter for the conditional "|" symbol.  Additionally, since order is not important in variable names here 
(:math:`H(R,x) = H(x,R)`), variables have a "*R* then *x* then *y*" default order applied to them. 
For example, 

   * :math:`H(R)` = *ContingencyTable.H_R*, Entropy of the repsonse variable.
   * :math:`H(x|R)` = *ContingencyTable.H_x_R*, Conditional entropy of *x* given a repsonse *R*.
   * :math:`I(R,x,y)` = *ContingencyTable.I_Rxy*, Mutual information across all parameters in a 3D system.
   * :math:`I(R,y|x)` = *ContingencyTable.I_Ry_x*, Conditional mutual information given *x* in a 3D system.

==============================================================
:mod:`ContingencyTable2D` -- Two Dimensional Contingency Table
==============================================================

.. automodule:: metasci.stats.ContingencyTable2D
   :members:
   :undoc-members:
   :inherited-members:

================================================================
:mod:`ContingencyTable3D` -- Three Dimensional Contingency Table
================================================================

.. automodule:: metasci.stats.ContingencyTable3D
   :members:
   :undoc-members:
   :inherited-members:
