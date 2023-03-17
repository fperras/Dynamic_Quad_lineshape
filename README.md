# Dynamic_Quad_lineshape
A program used to simulate static and MAS NMR spectra of quadrupolar nuclei experiencing motions

This is a simple C/C++ program utilizing the Eigen linear algebra library to simulate the
NMR spectra of half-integer quadrupolar nuclei in the precence of motions. Parameters
are written in a SIMPSON-style input file names input.txt, with an example present in 
the repository.

Sites can be added with individual chemical shift, CQ, eta, and tensor orientation and these
can be dynamically linked by giving them the same site index. This enables for the simultaneous
simulation of multiple sites with different sets of motions, or none at all.

The linux version in the bin directory uses an argument for the name of the input file as it is
intended to be run remotely.

The program uses OpenMP parallelization 
