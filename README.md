# Computational Fluid Dynamics
This fortran code is a implementation of a 3D pipe flow using Direct Forcing Immersed Boundary Method with FVM. Velocity and pressure fields are coupled to solve Navier Stokes equations. With the use of DFIB, geometry of a pipe can be created. Output data is in the form of binary file. Visualization can be achieve in "Plot3D" format with Paraview.

## Makefile
Use make to compile the program from source code and create executice file
'''bash
make
'''
to clean up object files and old data:
'''bash
make cleanall
'''

## Execution
Run in MPI
'''bash
mpirun -np (nprocs) ./main
'''

## Initial Author
Wei, Zi-Hsuan
