# Dynamic_Quad_lineshape
A program used to simulate static and MAS NMR spectra of quadrupolar nuclei experiencing motions

This is a simple C/C++ program utilizing the Eigen linear algebra library to simulate the
NMR spectra of half-integer quadrupolar nuclei in the precence of motions. Parameters
are written in a SIMPSON-style input file names input.txt, with an example present in 
the repository.

Different jump orientations are given with an orientation keywork followed by the three
Euler angles describing this orientation. kex is the hopping rate.

The program uses OpenMP parallelization 
