/*********************************
 
 Ant Colony Optimization algorithms (AS, ACS, EAS, RAS, MMAS) for CVRP
 
 Created by 孙晓奇 on 2016/10/8.
 Copyright © 2016年 xiaoqi.sxq. All rights reserved.
 
 Program's name: acovrp
 Purpose: Decompose the master problem into some subproblems, improve the time efficiency.
 Comment: 参考文章 D-Ants: Savings Based Ants divide and conquer the vehicle routing problem
 
 email: sunxq1991@gmail.com
 
 *********************************/

#ifndef parallelAco_h
#define parallelAco_h

#include <stdio.h>
#include <vector>
#include "problem.h"
#include "antColony.h"

using namespace std;

class ParallelAco: public AntColony {
public:
    ParallelAco(Problem *inst):AntColony(inst){}
    virtual ~ParallelAco();
    
    virtual void run_aco_iteration(void);
    void init_sub_pheromone(AntColony *sub_solver, Problem *master, Problem *sub, double ratio);
    void update_sub_best_pheromone(Problem *sub);
    void update_sub_to_master(Problem *master, Problem *sub, double ratio);
    
private:
    vector<Problem *> subs;            /* 多个子问题 */
    vector<RouteCenter *> route_centers;    /* each route's center info */
    
    void get_solution_centers(AntStruct *ant);
    void sort_route_centers();
    void decompose_problem(AntStruct *best_so_far_ant);
    void build_sub_problems(AntStruct *ant, const vector< vector<RouteCenter *> >& sub_problem_routes);
};

#endif /* parallelAco_h */
