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
#include "io.h"
#include "timer.h"


ParallelAco::~ParallelAco()
{
    for (int i = 0; i < route_centers.size(); i++) {
        delete route_centers[i]->coord;
        delete route_centers[i];
    }
}

/*
 * compute the center of gravity for best so far solution's each route
 */
void ParallelAco::get_solution_centers(AntStruct *ant)
{
    long int route_beg = 0;
    RouteCenter * center;
    route_centers.clear();   // 清空前一次迭代的数据
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
    
//    for (int i = 0; i < route_centers.size(); i++) {
//        printf("angle: %d (x,y): (%.2f,%.2f)\n",
//               route_centers[i]->angle, route_centers[i]->coord->x, route_centers[i]->coord->y);
//    }
//    printf("\n");
}

/*
 * Decompose the master problem by these routes' centers
 * choose a starting node randomly and the 
 * remaining nodes are sorted according to their polar angle
 */
void ParallelAco::decompose_problem(AntStruct *ant)
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
    
    build_sub_problems(ant, sub_problem_routes);
}

/*
 * 将主问题分解成多个独立的子问题，并初始化子问题结构体
 */
void ParallelAco::build_sub_problems(AntStruct *ant, const vector< vector<RouteCenter *> >& sub_problem_routes)
{
    subs.clear();   //清空上一次迭代数据
    
    Problem *sub, *master = instance;
    long int *tour = ant->tour;
    vector<RouteCenter *> routes;
    long int sub_best_length, sub_best_size;
    
    for (int p = 0; p < sub_problem_routes.size(); p++)
    {
        sub = new Problem(p+1);         // 子问题从1开始编号!!
        routes = sub_problem_routes[p];
        sub_best_length = 0;
        /*
         * 为了子问题能够直接使用串行蚁群算法，需要对nodes从0开始重新编号
         * i.e. nodes = {0,3,12,5,8} 重新编号后: {0,1,2,3,4}
         * 因此需要额外的数组记录编号前后的映射关系
         */
        sub->real_nodes.push_back(0);   // the depot
        sub->num_node = 1;
        for (long int  i = 0; i < routes.size(); i++)
        {
            for (long int j = routes[i]->beg; j < routes[i]->end; j++)
            {
                if (tour[j] != 0) {
                    sub->real_nodes.push_back(tour[j]);
                    sub->num_node++;
                }
                sub_best_length += master->distance[tour[j]][tour[j+1]];
            }
        }
        
        // 获取子问题的num_node之后便可以初始化
        init_sub_problem(master, sub);
        
        /* 
         * 由主问题计算出来的解作为子问题当前最优解
         * 需要在init_sub_problem()之后，best_so_far_ant才分配内存
         */
        sub->best_so_far_ant->tour_length = sub_best_length;
        sub->best_so_far_ant->tour_size = routes.size() + sub->num_node;
        
        subs.push_back(sub);
    }
}

/*
 * 子问题从主问题那里获取初始信息素
 * 由于信息素与tour_length相关，所以初始化子问题信息素时需要等比例转换
 */
void ParallelAco::init_sub_pheromone(Problem *master, Problem *sub, double ratio)
{
    long int ri, rj;
    
    for(long int i = 0; i < sub->num_node; i++) {
        ri = sub->real_nodes[i];
        for (long int j = 0; j < sub->num_node; j++) {
            rj = sub->real_nodes[j];
            sub->pheromone[i][j] = master->pheromone[ri][rj] / ratio;
        }
    }
    
//    print_pheromone(sub);
}

/*
 * 如果子问题获得更优解，则将子问题更新到主问题
 */
void ParallelAco::update_sub_to_master(Problem *master, Problem *sub, double ratio)
{
    printf("update master by sub-problem. (m-pid,s-pid,s-iter)=(%d,%d,%ld)\n", master->pid, sub->pid, sub->iteration);
    
    long int ri, rj, k;
    long int *sub_tour = sub->best_so_far_ant->tour;
    long int *master_tour = master->best_so_far_ant->tour;
    long int *tmp = new long int[2*master->num_node-1];
    bool *visited = new bool[master->num_node]();
    long int master_sz = master->best_so_far_ant->tour_size;
    long int sub_sz = sub->best_so_far_ant->tour_size;
    
    // 更新主问题信息素
    for(long int i = 0; i < sub->num_node; i++) {
        ri = sub->real_nodes[i];
        for (long int j = 0; j < sub->num_node; j++) {
            rj = sub->real_nodes[j];
            master->pheromone[ri][rj] = sub->pheromone[i][j] * ratio;
        }
    }
    
    // 更新主问题最优解
    for ( int i = 0; i < master_sz; i++) {
        tmp[i] = master_tour[i];
    }
    
    k = 0;
    for (int i = 0; i < sub_sz; i++) {
        ri = sub->real_nodes[sub_tour[i]];
        visited[ri] = TRUE;
        master_tour[k++] = ri;
    }
    
    for (int i = 0; i < master_sz; i++) {
        if (!visited[tmp[i]]) {
            while (tmp[i] != 0) {
                master_tour[k++] = tmp[i++];
            }
            master_tour[k++] = 0;
        }
    }
    master->best_so_far_ant->tour_length = compute_tour_length(master, master_tour, master_sz);
    
    g_best_so_far_time = elapsed_time(VIRTUAL);
    write_best_so_far_report(master);
    
    delete[] tmp;
    delete[] visited;
    
    check_solution(master, master_tour, master_sz);
//    print_solution(master, master_tour, master_sz);
    
//    print_pheromone(master);
}

/*
 * 
 */
void ParallelAco::run_aco_iteration()
{
    long int i, j;
    Problem *sub, *master = instance;
    AntColony *sub_solver;
    long int sub_best_length, master_best_length;
    double ratio;    // 主从问题信息素转换率
    
    // computer master problem
    for (i = 0; i < g_master_problem_iteration_num; i++) {
        this->AntColony::run_aco_iteration();
    }
    // 在子问题递归过程中可能会更新主问题最优解，该值用于计算主从问题信息素转换率ratio
    master_best_length = master->best_so_far_ant->tour_length;
    
    // compute the center of gravity for each route
    get_solution_centers(master->best_so_far_ant);
    
    // decompose the best solution into some subproblems using Sweep Algorithm
    decompose_problem(master->best_so_far_ant);
    
    for (i = 0; i < subs.size(); i++)
    {
        sub = subs[i];
        sub_solver = new AntColony(sub);
        
        sub_best_length = sub->best_so_far_ant->tour_length;
        // 主问题的信息素给子问题时，需要根据tour length做一步转化
        ratio = 1.0 * sub_best_length / master_best_length;
        
        // 根据主问题信息素初始化子问题信息素
        init_sub_pheromone(master, sub, 1.0);
        
        // 子问题递归
        for (j = 0; j < sub->max_iteration; j++)
        {
            sub_solver->AntColony::run_aco_iteration();
            
            // 子问题获得更优解, 则更新主问题
            if (sub->best_so_far_ant->tour_length < sub_best_length)
            {
                sub_best_length = sub->best_so_far_ant->tour_length;
                ratio =  1.0 * sub_best_length / master_best_length;
                update_sub_to_master(master, sub, 1.0);
                
//                print_solution(sub, sub->best_so_far_ant->tour, sub->best_so_far_ant->tour_size);
            }
            sub->iteration++;
        }
    }
}




