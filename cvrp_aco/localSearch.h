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
    int n_ants;
    bool ls_flag;
    double **distance;
    int num_node;
    int **nn_list;
    int nn_ls;
    bool dlb_flag;
    
    void two_opt_solution(int *tour, int tour_size);
    void swap(int *tour, int tour_size);
    int * generate_random_permutation( int n );
    void two_opt_single_route(int *tour, int rbeg, int rend, bool *dlb,
                              bool *route_node_map, int *tour_node_pos);
};

#endif /* localSearch_h */


