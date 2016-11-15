/*********************************
 Ant Colony Optimization algorithms (AS, ACS, EAS, RAS, MMAS) for CVRP
 
 Created by 孙晓奇 on 2016/10/8.
 Copyright © 2016年 xiaoqi.sxq. All rights reserved.
 
 Program's name: acovrp
 Purpose: vrp related procedures, distance computation, neighbour lists
 
 email: sunxq1991@gmail.com
 
 *********************************/
#ifndef   vrpHelper_h
#define   vrpHelper_h

#include <stdio.h>
#include <vector>
#include "problem.h"

#define RRR            6378.388
#ifndef PI             /* as in stroustrup */
#define PI             3.14159265358979323846
#endif


double compute_tour_length(Problem *instance, int *t, int t_sz);
double compute_route_length(Problem *instance, int *route, int route_size);
double **compute_distances(Problem *instance);
int ** compute_nn_lists (Problem *instance);
void compute_route_centers(Problem *instance, int *tour, const vector<RouteCenter *>& centers);

#endif
