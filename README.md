# cvrp_aco
Ant colony optimization(ranked aco) for the capacitated vehicle routing problem.

######################################################
######    Hybrid ACO algorithms for the CVRP ########
######################################################
=========
CONTENTS
=========


<h3> The main control routines, main: </h3>
* main.cpp

<h3> CVRP problem model, problem init/exit, solution check: </h3> 
problem.cpp
problem.h

<h3> Implementation of procedures for ants' behaviour: </h3> 
antColony.cpp
antColony.h

<h3> Parallel version of aco.Decompose the master problem into some subproblems, improve the time efficiency: </h3>
parallelAco.cpp
parallelAco.h

<h3> Simulated Annealing algorithm is used to improve the effectiveness of aco. </h3>
simulatedAnnealing.cpp
simulatedAnnealing.h

<h3> Performs three types of neighborhood search: sequence inversion, insertion, and exchange: </h3>
neighbourSearch.cpp
neighbourSearch.h
move.cpp
move.h

<h3> mainly input / output / statistic routines: </h3>
io.cpp
io.h

<h3> Local search procedures: </h3>
localSearch.cpp
localSearch.h

<h3> Additional useful / helping procedure: </h3>
utilities.c
utilities.h

<h3> VRP related procedures, distance computation, neighbour lists: </h3>
vrpHelper.cpp
vrpHelper.h

<h3> Time measurement:
timer.h 
unix_timer.c : in case you want to use rusage() instead, edit the
Makefile to use this one or compile with 'make TIMER=unix'

Makefile

=====
Code
=====

The program is developed under XCode, so it can be executed by XCode directly.

This program can also be compiled by the GNU g++. To run this program, you need to do the folllowing three steps:

* make all;
* g++ *.o -o solver -lpthread;
* ./solver filename


