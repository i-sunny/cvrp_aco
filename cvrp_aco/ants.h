/*********************************
 
 Ant Colony Optimization algorithms (AS, ACS, EAS, RAS, MMAS) for CVRP
 
 Created by 孙晓奇 on 2016/10/8.
 Copyright © 2016年 xiaoqi.sxq. All rights reserved.
 
 Program's name: acovrp
 Purpose: header of procedures for ants' behaviour
 
 email: sunxq1991@gmail.com
 
 *********************************/


#define HEURISTIC(m,n)     (1.0 / ((double) instance.distance[m][n] + 0.1))
/* add a small constant to avoid division by zero if a distance is
 zero */

#define EPSILON            0.00000000000000000000000000000001

#define MAX_ANTS       1024    /* max no. of ants */
#define MAX_NEIGHBOURS 512     /* max. no. of nearest neighbours in candidate set */

/*
 * ant tour consists of some single route.
 * a single route begin with depot(0) and end with depot(0),
 * i.e. tour = [0,1,4,2,0,5,3,0] toute1 = [0,1,4,2,0] route2 = [0,5,3,0]
 */
typedef struct {
    long int  *tour;
    long int  tour_size;     /* 路程中经过的点个数（depot可能出现多次） */
    long int  tour_length;   /* 车辆路程行走长度 */
    char      *visited;
} ant_struct;

extern ant_struct *ant;      /* this (array of) struct will hold the colony */
extern ant_struct *best_so_far_ant;   /* struct that contains the best-so-far ant */
extern ant_struct *iteration_best_ant;     /* 当前迭代表现最好的蚂蚁 */

extern double   **pheromone;      /* pheromone matrix, one entry for each arc */
extern double   **total_info;     /* combination of pheromone and heuristic information */

extern double   *prob_of_selection;


extern long int n_ants;      /* number of ants */
extern long int nn_ants;     /* length of nearest neighbor lists for the ants'
                                solution construction */

extern double rho;           /* parameter for evaporation */
extern double alpha;         /* importance of trail */
extern double beta;          /* importance of heuristic evaluate */
extern double q_0;           /* probability of best choice in tour construction */

extern long int ras_flag;    /* = 1, run rank-based version of ant system */
extern long int ras_ranks;       /* additional parameter for rank-based version of ant
				    system */

extern long int u_gb;            /* every u_gb iterations update with best-so-far ant;
				    parameter used by MMAS for scheduling best-so-far update
				 */

/* Pheromone manipulation etc. */

void init_pheromone_trails ( double initial_trail );

void evaporation ( void );

void evaporation_nn_list ( void );

void global_update_pheromone ( ant_struct *a );

void global_update_pheromone_weighted ( ant_struct *a, long int weight );

void compute_total_information( void );

void compute_nn_list_total_information( void );

/* Ants' solution construction */

void ant_empty_memory( ant_struct *a );

void init_ant_place( ant_struct *a , long int phase);

long int choose_best_next( ant_struct *a, long int phase );

long int neighbour_choose_best_next( ant_struct *a, long int phase );

void choose_closest_next( ant_struct *a, long int phase );

long int neighbour_choose_and_move_to_next( ant_struct *a, long int phase);

/* Auxiliary procedures related to ants */

long int find_best ( void );

long int find_worst( void );

void copy_solution_from_to(ant_struct *a1, ant_struct *a2);

void allocate_ants ( void );

