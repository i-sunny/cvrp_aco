/*********************************
 Ant Colony Optimization algorithms (AS, ACS, EAS, RAS, MMAS) for CVRP
 
 Created by 孙晓奇 on 2016/10/8.
 Copyright © 2016年 xiaoqi.sxq. All rights reserved.
 
 Program's name: acovrp
 Purpose: vrp related procedures, distance computation, neighbour lists
 
 email: sunxq1991@gmail.com
 
 *********************************/
#ifndef   HEAD_VRP_H
#define   HEAD_VRP_H

#define RRR            6378.388
#ifndef PI             /* as in stroustrup */
#define PI             3.14159265358979323846
#endif

struct point {
  double x;
  double y;
  long int demand;     /* 每个配送点需求 */
};

struct problem{
  char          name[LINE_BUF_LEN];      	 /* instance name */
  char          edge_weight_type[LINE_BUF_LEN];  /* selfexplanatory */
  long int      optimum;                /* optimal tour length if known, otherwise a bound */
  long int      num_node;                      /* number of nodes, depot included */
  long int      n_near;                 /* number of nearest neighbors */
  struct point  *nodeptr;               /* array of structs containing coordinates of nodes */
  long int      **distance;	        /* distance matrix: distance[i][j] gives distance 
					   between node i und j */
  long int      **nn_list;              /* nearest neighbor list; contains for each node i a
                                           sorted list of n_near nearest neighbors */
  long int      vehicle_capacity;       /* 车辆最大装载量 */
};

extern struct problem instance;

extern long int num_node;                 /* number of nodes */

extern long int vehicle_capacity;         /** the max load of the vehicle */

extern long int  (*distance)(long int, long int);  /* pointer to function returning distance */


long int round_distance (long int i, long int j);

long int ceil_distance (long int i, long int j);

long int geo_distance (long int i, long int j);

long int att_distance (long int i, long int j);

long int compute_tour_length( long int *t, long int t_sz);

long int **compute_distances(void);

long int ** compute_nn_lists ( void );

int vrp_check_solution(const long int *tour, long int tour_size);

int vrp_check_route(const long int *tour, long int rbeg, long int rend);

#endif
