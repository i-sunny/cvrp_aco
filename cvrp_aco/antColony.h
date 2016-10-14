/*********************************
 
 Ant Colony Optimization algorithms (AS, ACS, EAS, RAS, MMAS) for CVRP
 
 Created by 孙晓奇 on 2016/10/8.
 Copyright © 2016年 xiaoqi.sxq. All rights reserved.
 
 Program's name: acovrp
 Purpose: header of procedures for ants' behaviour
 
 email: sunxq1991@gmail.com
 
 *********************************/

#ifndef antColony_h
#define antColony_h

#define HEURISTIC(m,n)     (1.0 / ((double) distance[m][n] + 0.1))
/* add a small constant to avoid division by zero if a distance is
 zero */

#include "problem.h"
#include "localSearch.h"

#define EPSILON            0.00000000000000000000000000000001

#define MAX_ANTS       1024    /* max no. of ants */
#define MAX_NEIGHBOURS 512     /* max. no. of nearest neighbours in candidate set */

class AntColony {
public:
    
    Problem *instance;
    LocalSearch *local_search;
    
    // 直接引入，减少一次指针使用，加快执行速度
    AntStruct *ants;
    AntStruct *best_so_far_ant;
    
    long int **distance;
    bool     *demand_meet_node_map;
    double   *prob_of_selection;
    double   **pheromone;
    double   **total_info;
    long int num_node;
    long int n_ants;               /* number of ants */
    long int nn_ants;              /* length of nearest neighbor lists for the ants'
                                    solution construction */
    long int **nn_list;
    
    Point    *nodeptr;
    long int vehicle_capacity;
    
    bool ls_flag;               /* indicates whether and which local search is used */
    
    
    AntColony(Problem *instance);
    virtual ~AntColony();
    
    virtual void run_aco_iteration(void);
    void init_aco();
    void exit_aco();
    
    void construct_ant_solution(AntStruct *ant);
    void construct_solutions( void );
    void do_local_search( void );
    void ras_update( void );
    void pheromone_trail_update( void );
    void update_statistics( void );
    
    /* Pheromone manipulation etc. */
    void init_pheromone_trails ( double initial_trail );
    void evaporation ( void );
    void evaporation_nn_list ( void );
    void global_update_pheromone ( AntStruct *a );
    void global_update_pheromone_weighted ( AntStruct *a, long int weight );
    void compute_total_information( void );
    void compute_nn_list_total_information( void );
    
    /* Ants' solution construction */
    void ant_empty_memory( AntStruct *a );
    void init_ant_place( AntStruct *a , long int phase);
    long int choose_best_next( AntStruct *a, long int phase );
    long int neighbour_choose_best_next( AntStruct *a, long int phase );
    void choose_closest_next( AntStruct *a, long int phase );
    long int neighbour_choose_and_move_to_next( AntStruct *a, long int phase);
    
    /* Auxiliary procedures related to ants */
    long int find_best ( void );
    long int find_worst( void );
    void copy_solution_from_to(AntStruct *a1, AntStruct *a2);
};

#endif /* antColony_h */
