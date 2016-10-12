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

#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#include "parallelAco.h"
#include "vrpHelper.h"
#include "utilities.h"


ParallelAco::~ParallelAco()
{
    for (int i = 0; i < route_centers.size(); i++) {
        delete route_centers[i]->coord;
        delete route_centers[i];
    }
}

void init()
{

}

/*
 * compute the center of gravity for best so far solution's each route
 */
void ParallelAco::get_solution_centers(AntStruct *ant)
{
    long int route_beg = 0;
    RouteCenter * center;
    for (int i = 1; i < ant->tour_size; i++) {
        if (ant->tour[i] == 0) {
            center = new RouteCenter();
            center->beg = route_beg;
            center->end = i;
            route_centers.push_back(center);
            route_beg = i;
        }
    }
    compute_route_centers(instance, ant->tour, route_centers);
}

bool cmp(const RouteCenter* ra, const RouteCenter* rb)
{
    return ra->angle < rb->angle;
}

/*
 * sort routes according to their polar angle with the depot
 */
void ParallelAco::sort_route_centers()
{
    long int x0 = nodeptr[0].x, y0 = nodeptr[0].y, dx, dy;
    
    for (int i = 0; i < route_centers.size(); i++) {
        dx = route_centers[i]->coord->x - x0;
        dy = route_centers[i]->coord->y - y0;
        if (dx > 0) {
            route_centers[i]->angle = (short)(atan(1.0*dy/dx)*180/PI);
        } else if (dx < 0){
            route_centers[i]->angle = (short)(atan(1.0*dy/dx)*180/PI) + 180;
        } else {
            route_centers[i]->angle = dy > 0 ? 90 : -90;
        }
    }
    sort(route_centers.begin(), route_centers.end(), cmp);
    
    for (int i = 0; i < route_centers.size(); i++) {
        printf("angle: %d (x,y): (%.2f,%.2f)\n",
               route_centers[i]->angle, route_centers[i]->coord->x, route_centers[i]->coord->y);
    }
    printf("\n");
}

/*
 * Decompose the master problem by these routes' centers
 * choose a starting node randomly and the 
 * remaining nodes are sorted according to their polar angle
 */
void ParallelAco::decompose_problem(void)
{
    // random start pos from [0, route_num)
    long int rnd_beg = (long int)(ran01(&seed) * (route_centers.size() - 1));
    /* 一个子问题包含的routes数目 */
    long int sub_problem_route_num = (long int)(route_centers.size() / g_num_sub_problems);
    /* 
     * 剩余的routes, 前remainder个子问题，每个子问题多分配一个，
     * 如此尽量保证分配均匀. 
     * 每个子问题分配的子routes数量为: sub_problem_route_num 或 sub_problem_route_num + 1
     */
    long int remainder = route_centers.size() % g_num_sub_problems;
    
    vector< vector<RouteCenter *> > sub_problem_routes;
    vector<RouteCenter *> tmp;

    sort_route_centers();
    for(int i = 0; i < rnd_beg; i++) {
        route_centers.insert(route_centers.begin(), *(route_centers.end()-1));
        route_centers.pop_back();
    }
    
    // 计算子问题分配的routes
    for (int i = 0; i < route_centers.size(); i++) {
        tmp.push_back(route_centers[i]);
        if (tmp.size() == sub_problem_route_num + (remainder > 0)) {
            if (remainder > 0) {
                remainder--;
            }
            sub_problem_routes.push_back(tmp);
            tmp.clear();
        }
    }
    
    build_sub_problem(sub_problem_routes);
}

void ParallelAco::build_sub_problem(const vector< vector<RouteCenter *> >& sub_problem_routes)
{
    Problem *sub_inst;
    
}

/*
 * 
 */
void ParallelAco::run_aco_iteration()
{
    // computer master problem
    for (int i = 0; i < g_master_problem_iteration_num; i++) {
        this->AntColony::run_aco_iteration();
    }
    
    // compute the center of gravity for each route
    get_solution_centers(instance->best_so_far_ant);
    
    // decompose the best solution into some subproblems using Sweep Algorithm
    decompose_problem();
}




