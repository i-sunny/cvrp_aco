/*********************************
 
 Ant Colony Optimization algorithms (AS, ACS, EAS, RAS, MMAS, BWAS) for CVRP
 
 Created by 孙晓奇 on 2016/10/8.
 Copyright © 2016年 xiaoqi.sxq. All rights reserved.
 
 Program's name: acovrp
 Purpose: header file for local search routines
 
 email: sunxq1991@gmail.com
 
 *********************************/

#ifndef localSearch_h
#define localSearch_h

#include "problem.h"

class LocalSearch {
public:
    
    LocalSearch(Problem *instance);
    ~LocalSearch(){}
    
    void do_local_search(void);
    void do_local_search(AntStruct *ant);
    
private:
    Problem *instance;
    
    // 直接引入，减少一次指针使用，加快执行速度
    AntStruct *ants;
    long int n_ants;
    bool ls_flag;
    double **distance;
    long int num_node;
    long int **nn_list;
    long int nn_ls;
    bool dlb_flag;
    
    void two_opt_solution(long int *tour, long int tour_size);
    void swap(long int *tour, long int tour_size, long int *path_load);
    long int * generate_random_permutation( long int n );
    void two_opt_single_route(long int *tour, long int rbeg, long int rend, long int *dlb,
                              long int *route_node_map, long int *tour_node_pos);
};

#endif /* localSearch_h */


