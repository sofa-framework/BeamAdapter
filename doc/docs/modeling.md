This page presents the mechanical basis of the beam model and how the internal forces are computed.

# Modeling
Beam elements are used to model instruments (or anatomical structures) for which the length is greater than the other transverse dimensions. 
Examples of objects with this geometry are abundant in medical simulation. 
The particular nature of such objects generally leads to large geometric deformations, 
which is notorious for being computationally demanding. 

The approach that is used in the plugin is initially based on a linear beam analysis. 
But it extends this representation by a series of optimizations particularly suited for real-time animation. 
By using a corotational approach, our model can handle the important geometric non-linearity due to large changes in the shape of the object.

## Corotational beam model
To model the deformation of any solid body whose geometry and mechanical characteristics are similar to a wire, 
rod or beam, we use a representation based on three-dimensional beam theory *see Przemieniecki (1985)*, where the 
elementary stiffness matrix Ke is a 12 × 12 symmetric matrix that relates angular and spacial positions of each end 
of a beam element to the forces and torques applied to them:

<img src="../doc/docs/Matrix.jpg" align="left" width="700"/> 
<br/><br/><br/><br/><br/><br/><br/><br/><br/><br/><br/><br/>

parameters:

<img src="../doc/docs/Matrix_param.jpg" align="left" width="700"/> <br/>
<br/><br/><br/><br/><br/>

The assumption of the corotational model is that the deformations remain ”small” at the level of each element. 
Thus, a local frame is defined for each beam. Thus, the force  **f<sub>e</sub>** at the level of the element is equal to:

**f<sub>e</sub> = K<sub>e</sub> (u - u<sub>0</sub>)** (1)

where **u** and **u<sub>0</sub>** reflect respectively the actual and the initial configuration of the beam,
in the local frame of the beam.
Since the stiffness matrix **K<sub>e</sub>** is initially calculated in local coordinates, 
it is necessary to introduce transformation matrices changing the reference frame from a local to a global coordinate system. 
In order to determine the stiffness property of the complete structure, 
a common reference frame must be established for all unassembled structural elements so that all the displacements 
and their corresponding forces will be referred to a common (global) coordinate system. 
We need a matrix relationship between the variation of the position δ **u** in the local coordinate system and the variation of the node position 
δ **q** in the global coordinates.


<img src="../doc/docs/BeamFrames.jpg" align="left" width="700"/> 
<br/><br/><br/><br/><br/><br/><br/><br/><br/><br/><br/><br/><br/><br/><br/><br/>

This relationship is expressed by the matrix equation:


δ **u** = **Λ(q)**     δ **q** (2)

where **Λ(q)** is a matrix obtained from the direction cosines of angles between the local and global coordinate systems.
The linearization of the element force-displacement equation is obtained in global coordinates:

δ **F(q)** = [**Λ<sup>T</sup>(q)** **Ke** **Λ(q)**] δ **q**   (3)

The model can be combined with the computation of a mass matrix to obtain mechanical dynamic model. 
We then use the equation 3 to compute the system. 
One important feature of this model is that the interpolation is performed segment by segment between only two frames. 
Additionally, these frames are the independent Degrees of Freedom (DoFs) of the system.





## References
*Przemieniecki (1985) J. Przemieniecki. Theory of matrix structural analysis. McGraw-Hill, 1985.*