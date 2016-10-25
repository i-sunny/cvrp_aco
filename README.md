# cvrp_aco
Ant colony optimization(ranked aco) for the capacitated vehicle routing problem.

######################################################
######    Hybrid ACO algorithms for the CVRP ########
######################################################

        Version: 1.0
        Author:  xiaoqi.sxq
        Copyright (c) xiaoqi.sxq, 2016


=========
CONTENTS
=========

The main control routines, main:
main.cpp

CVRP problem model, problem init/exit, solution check:
problem.cpp
problem.h

Implementation of procedures for ants' behaviour:
antColony.cpp
antColony.h

Parallel version of aco.Decompose the master problem into some subproblems, improve the time efficiency:
parallelAco.cpp
parallelAco.h

Simulated Annealing algorithm is used to improve the effectiveness of aco.
simulatedAnnealing.cpp
simulatedAnnealing.h

Performs three types of neighborhood search: sequence inversion, insertion, and exchange:
neighbourSearch.cpp
neighbourSearch.h
move.cpp
move.h

mainly input / output / statistic routines:
io.cpp
io.h

Local search procedures:
localSearch.cpp
localSearch.h

Additional useful / helping procedure:
utilities.c
utilities.h

VRP related procedures, distance computation, neighbour lists:
vrpHelper.cpp
vrpHelper.h

Time measurement:
timer.h 
unix_timer.c : in case you want to use rusage() instead, edit the
Makefile to use this one or compile with 'make TIMER=unix'

Makefile

=====
Code
=====

The program is developed under XCode, so it can be executed by XCode directly.

This program can also be compiled by the GNU g++. To run this program, you need to do the folllowing three steps:

make all;
g++ *.o -o solver -lpthread;
./solver filename


