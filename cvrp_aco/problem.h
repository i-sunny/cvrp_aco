/*********************************
 Ant Colony Optimization algorithms (AS, ACS, EAS, RAS, MMAS) for CVRP
 
 Created by 孙晓奇 on 2016/10/8.
 Copyright © 2016年 xiaoqi.sxq. All rights reserved.
 
 Program's name: acovrp
 Purpose: vrp problem model, problem init/exit, solution check
 
 email: sunxq1991@gmail.com
 
 *********************************/

#ifndef Problem_h
#define Problem_h

#include <stdio.h>

#define LINE_BUF_LEN     255

/*************** global variables *********************/
extern long int seed;

extern double   g_max_runtime;                    /* maximal allowed run time of a try  */
extern double   g_best_so_far_time;               /* 最优解出现的时间 */
extern long int g_best_solution_iter;             /* iteration in which best solution is found */

extern long int g_master_problem_iteration_num;   /* 每次外循环，主问题蚁群的迭代的次数 */
extern long int g_sub_problem_iteration_num;      /* 每次外循环，子问题蚁群的迭代的次数 */
extern long int g_num_sub_problems;               /* 拆分子问题个数 */

extern double rho;           /* parameter for evaporation */
extern double alpha;         /* importance of trail */
extern double beta;          /* importance of heuristic evaluate */
extern long int ras_ranks;   /* additional parameter for rank-based version of ant system */

/****************** data struct ***********************/
enum DistanceTypeEnum {
    DIST_EUC_2D, DIST_CEIL_2D, DIST_GEO, DIST_ATT
};

struct Point {
    double x;
    double y;
    long int demand;     /* 每个配送点需求 */
};

struct RouteCenter {
    long int beg;                   /* route在tour中的开始位置 */
    long int end;                   /* route在tour中的结束位置 */
    Point    *cp;                   /* center pointroute */
};

/*
 * ant tour consists of some single route.
 * a single route begin with depot(0) and end with depot(0),
 * i.e. tour = [0,1,4,2,0,5,3,0] toute1 = [0,1,4,2,0] route2 = [0,5,3,0]
 */
struct AntStruct {
    long int  *tour;
    long int  tour_size;     /* 路程中经过的点个数（depot可能出现多次） */
    long int  tour_length;   /* 车辆路程行走长度 */
    char      *visited;
};


struct Problem {
    char          name[LINE_BUF_LEN];      	 /* instance name */
    char          edge_weight_type[LINE_BUF_LEN];  /* selfexplanatory */
    DistanceTypeEnum dis_type;               /* 用于决定使用哪种距离方式 */
    long int      optimum;                /* optimal tour length if known, otherwise a bound */
    long int      num_node;               /* number of nodes, depot included, numd_node = 1(depot) + num of target nodes*/
    Point         *nodeptr;               /* array of structs containing coordinates of nodes */
    long int      **distance;             /* distance matrix: distance[i][j] gives distance
                                           between node i und j */
    long int      **nn_list;              /* nearest neighbor list; contains for each node i a
                                           sorted list of n_near nearest neighbors */
    long int      vehicle_capacity;       /* 车辆最大装载量 */
    long int      *demand_meet_node_map;  /** 所有可配送的点(单次route中，目前车辆可以仍可配送的点) */
    
    
    /*----- local search -----*/
    long int ls_flag;          /* indicates whether and which local search is used */
    long int nn_ls;            /* maximal depth of nearest neighbour lists used in the
                                local search */
    long int dlb_flag;         /* flag indicating whether don't look bits are used. I recommend
                                to always use it if local search is applied */
    
    long int iteration;           /* counter of number iterations */
    long int max_iteration;       /* maximum number of iterations */
    
    /*----- ant info -----*/
    AntStruct *ants;                   /* this (array of) struct will hold the colony */
    AntStruct *best_so_far_ant;        /* struct that contains the best-so-far ant */
    AntStruct *iteration_best_ant;     /* 当前迭代表现最好的蚂蚁 */
    
    double   **pheromone;               /* pheromone matrix, one entry for each arc */
    double   **total_info;              /* combination of pheromone and heuristic information */
    
    double   *prob_of_selection;        /* 依概率选择下一个node */
    
    long int n_ants;                    /* number of ants */
    long int nn_ants;                   /* length of nearest neighbor lists for the ants'
                                         solution construction */
};

Problem * init_master_problem(const char *filename);
void exit_master_problem(Problem *instance);
void init_sub_problem(Problem *instance);
void exit_sub_problem(Problem *instance);
int check_solution(Problem *instance, const long int *tour, long int tour_size);
int check_route(Problem *instance, const long int *tour, long int rbeg, long int rend);

#endif /* Problem_h */
