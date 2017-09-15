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
![Image of myMatrix](doc/docs/matrix.jpg)








*Przemieniecki (1985) J. Przemieniecki. Theory of matrix structural analysis. McGraw-Hill, 1985.*