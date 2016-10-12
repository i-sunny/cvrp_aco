///*********************************
// 
// Ant Colony Optimization algorithms (AS, ACS, EAS, RAS, MMAS) for CVRP
// 
// Created by 孙晓奇 on 2016/10/8.
// Copyright © 2016年 xiaoqi.sxq. All rights reserved.
// 
// Program's name: acovrp
// Purpose: Decompose the master problem into some subproblems, improve the time efficiency.
// Comment: 参考文章 D-Ants: Savings Based Ants divide and conquer the vehicle routing problem
// 
// email: sunxq1991@gmail.com
// 
// *********************************/
//
//#include <stdio.h>
//#include <math.h>
//#include <limits.h>
//#include <assert.h>
//#include <string.h>
//#include <stdlib.h>
//#include <time.h>
//
////#include "ants.h"
//#include "utilities.h"
//#include "InOut.h"
//#include "vrp.h"
//#include "acovrp.h"
//#include "timer.h"
//#include "ls.h"
//
//long int parallel_flag;
//struct route_center *route_centers;            /* each route's center info */
//
//
//void init()
//{
//    long int max_route_num = num_node - 1;
//    if((route_centers = malloc(sizeof(struct route_center) * max_route_num +
//                         sizeof(struct route_center *) * max_route_num)) == NULL){
//        printf("Out of memory, exit.");
//        exit(1);
//    }
//    for (int i = 0; i < max_route_num; i++) {
//        route_centers[i].cp = malloc(sizeof(struct point));
//    }
//}
//
///*
// * compute the center of gravity for best so far solution's each route
// */
//void get_bsf_solution_centers(void)
//{
//    long int route_beg = 0;
//    long int route_cnt = 0;
//    for (int i = 0; i < best_so_far_ant->tour_size; i++) {
//        if (best_so_far_ant->tour[i] == 0) {
//            route_centers[route_cnt].beg = route_beg;
//            route_centers[route_cnt].end = i;
//            route_cnt++;
//            route_beg = i;
//        }
//        route_beg++;
//    }
//    compute_tour_centers(best_so_far_ant->tour, route_centers, route_cnt);
//}
//
//void decompose_master_problem()
//{
//    
//}
//
///*
// * 
// */
//void divide_and_conquer_runner()
//{
//    // computer master problem
//    for (int i = 0; i < master_problem_iteration_num; i++) {
//        aco_iteration_runner();
//    }
//    
//    // compute the center of gravity for each route
//    get_bsf_solution_centers();
//    
//    // decompose the best solution into some subproblems using Sweep Algorithm
//    
//    
//}


