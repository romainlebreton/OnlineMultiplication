OnlineMultiplication
====================

Implementation of several online algorithm for polynomial multiplication

### Authors
Romain Lebreton, Université Montpellier II, Montpellier, France
Email: lebreton@lirmm.fr

Éric Schost, Western University, London, Ontario, Canada
Email: eschost@uwo.ca

### Note
This code corresponds to the implementation of the paper "A simple and fast online power series multiplication and its analysis" from the same authors 
(https://hal.archives-ouvertes.fr/lirmm-00867279/)

### How to compile 
We have a basic MAKEFILE that will compile the two programs 'test' and 'bench', provided that you have NTL and GMP installed on your computer.
You only have to run 'make' to compile both programs.

The program 'test' will check the correctness of computations. The program 'bench' launches a series of benchmark that illustrates the performance of the different algorithms.


