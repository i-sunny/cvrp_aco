/*********************************
Ant Colony Optimization algorithms (AS, ACS, EAS, RAS, MMAS) for CVRP

Created by 孙晓奇 on 2016/10/8.
Copyright © 2016年 xiaoqi.sxq. All rights reserved.

Program's name: acovrp
Purpose: mainly input / output / statistic routines

email: sunxq1991@gmail.com

*********************************/

#ifndef inout_h
#define inout_h

#include <stdio.h>
#include "problem.h"


Point * read_instance_file(Problem *instance, const char *vrp_file_name);
const char* parse_commandline (long int argc, char *argv []);

void print_solution(Problem *instance, long int *tour, long int tour_size);
void print_single_route(Problem *instance, long int *route, long int route_size);
void print_probabilities(Problem *instance);
void print_solution_to_file(Problem *instance, FILE *file, long int *tour, long int tour_size);

void init_report(Problem *instance);
void exit_report(Problem *instance);
void write_best_so_far_report(Problem *instance);
void write_iter_report(Problem *instance);

#endif /* ants_h */
