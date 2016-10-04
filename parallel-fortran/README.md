# parallel-fortran
Fortran 90 code written with Dr. Pinaki Sengupta to implement the Stochastic Series Expansion method for quantum Monte Carlo evaluation of the spin S=1 Heisenberg model. This code utilized MPI directives to perform parallel tempering on high performance computer clusters (multiple points in parameter space are evaluated in parallel with swapping of state configurations to reduce autocorrelation of the Markov chain data).

## random.f90
Module for pseudo-random number generation using the `mzran` generator of Marsaglia and Zaman (1994).

## main.f90
Main driver of the code. Reads input file, spawns threads, performs parallel tempering, and writes output.

## hb.f90
High-level Monte Carlo steps for each thread, along with setup of vertex weights (i.e. lookup tables for probability distributions).

## upd.f90
Low-level execution of Monte Carlo steps: cycles through operator string for diagonal updates, then creates a linked list of active operators for use during the off-diagonal loop update.
