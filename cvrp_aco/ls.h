/*********************************
 
 Ant Colony Optimization algorithms (AS, ACS, EAS, RAS, MMAS, BWAS) for CVRP
 
 Created by 孙晓奇 on 2016/10/8.
 Copyright © 2016年 xiaoqi.sxq. All rights reserved.
 
 Program's name: acovrp
 Purpose: header file for local search routines
 
 email: sunxq1991@gmail.com
 
 *********************************/


extern long int ls_flag;

extern long int nn_ls; 

extern long int dlb_flag;

void two_opt_solution(long int *tour, long int tour_size);

void two_opt_single_route(long int *tour, long int route_begin, long int route_end,
                          long int *dlb, long int *route_node_map, long int *tour_node_pos);
