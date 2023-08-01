# MATLAB CFD Code
Reza Sheikhi

University of Connecticut

This is a simple educational CFD code developed based on "Numerical Heat Transfer and Fluid Flow" book by S. Patankar (Chapter 6).
It incorporates an explicit spatial second order, temporal first order accurate numerical scheme along with staggered grid arrangement.
The boundary conditions are implemented using ghost cells. The pressure velocity link is handled using the Semi-Implicit Method for Pressure Linked Equations (SIMPLE) method.

Two cases are available:
* Duct-Flow: Incompressible flow in a 2D duct
* Mixing-layer: Incompressible flow in a spatial mixing layer
