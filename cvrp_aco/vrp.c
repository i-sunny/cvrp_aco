/*********************************
Ant Colony Optimization algorithms (AS, ACS, EAS, RAS, MMAS) for CVRP

Created by 孙晓奇 on 2016/10/8.
Copyright © 2016年 xiaoqi.sxq. All rights reserved.

Program's name: acovrp
Purpose: vrp related procedures, distance computation, neighbour lists

email: sunxq1991@gmail.com

*********************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <assert.h>

#include "InOut.h"
#include "vrp.h"
#include "ants.h"
#include "ls.h"
#include "utilities.h"

#define M_PI 3.14159265358979323846264

/* number of nodes
 * note that: numd_node = 1(depot) + num of target nodes
 */
long int num_node;

long int vehicle_capacity;         /** the max load of the vehicle */

struct problem instance;

static double dtrunc (double x)
{
    int k;

    k = (int) x;
    x = (double) k;
    return x;
}

long int  (*distance)(long int, long int);  /* function pointer */

/*    
      FUNCTION: the following four functions implement different ways of 
                computing distances for VRPLIB instances
      INPUT:    two node indices
      OUTPUT:   distance between the two nodes
*/

long int round_distance (long int i, long int j) 
/*    
      FUNCTION: compute Euclidean distances between two nodes rounded to next 
                integer for VRPLIB instances
      INPUT:    two node indices
      OUTPUT:   distance between the two nodes
      COMMENTS: for the definition of how to compute this distance see VRPLIB
*/
{
    double xd = instance.nodeptr[i].x - instance.nodeptr[j].x;
    double yd = instance.nodeptr[i].y - instance.nodeptr[j].y;
    double r  = sqrt(xd*xd + yd*yd) + 0.5;

    return (long int) r;
}

long int ceil_distance (long int i, long int j) 
/*    
      FUNCTION: compute ceiling distance between two nodes rounded to next 
                integer for VRPLIB instances
      INPUT:    two node indices
      OUTPUT:   distance between the two nodes
      COMMENTS: for the definition of how to compute this distance see VRPLIB
*/
{
    double xd = instance.nodeptr[i].x - instance.nodeptr[j].x;
    double yd = instance.nodeptr[i].y - instance.nodeptr[j].y;
    double r  = sqrt(xd*xd + yd*yd);

    return (long int)(ceil (r));
}

long int geo_distance (long int i, long int j) 
/*    
      FUNCTION: compute geometric distance between two nodes rounded to next 
                integer for VRPLIB instances
      INPUT:    two node indices
      OUTPUT:   distance between the two nodes
      COMMENTS: adapted from concorde code
                for the definition of how to compute this distance see VRPLIB
*/
{
    double deg, min;
    double lati, latj, longi, longj;
    double q1, q2, q3;
    long int dd;
    double x1 = instance.nodeptr[i].x, x2 = instance.nodeptr[j].x, 
	y1 = instance.nodeptr[i].y, y2 = instance.nodeptr[j].y;

    deg = dtrunc (x1);
    min = x1 - deg;
    lati = M_PI * (deg + 5.0 * min / 3.0) / 180.0;
    deg = dtrunc (x2);
    min = x2 - deg;
    latj = M_PI * (deg + 5.0 * min / 3.0) / 180.0;

    deg = dtrunc (y1);
    min = y1 - deg;
    longi = M_PI * (deg + 5.0 * min / 3.0) / 180.0;
    deg = dtrunc (y2);
    min = y2 - deg;
    longj = M_PI * (deg + 5.0 * min / 3.0) / 180.0;

    q1 = cos (longi - longj);
    q2 = cos (lati - latj);
    q3 = cos (lati + latj);
    dd = (int) (6378.388 * acos (0.5 * ((1.0 + q1) * q2 - (1.0 - q1) * q3)) + 1.0);
    return dd;

}

long int att_distance (long int i, long int j) 
/*    
      FUNCTION: compute ATT distance between two nodes rounded to next 
                integer for VRPLIB instances
      INPUT:    two node indices
      OUTPUT:   distance between the two nodes
      COMMENTS: for the definition of how to compute this distance see VRPLIB
*/
{
    double xd = instance.nodeptr[i].x - instance.nodeptr[j].x;
    double yd = instance.nodeptr[i].y - instance.nodeptr[j].y;
    double rij = sqrt ((xd * xd + yd * yd) / 10.0);
    double tij = dtrunc (rij);
    long int dij;

    if (tij < rij)
        dij = (int) tij + 1;
    else
        dij = (int) tij;
    return dij;
}



long int ** compute_distances(void)
/*    
      FUNCTION: computes the matrix of all intercity distances
      INPUT:    none
      OUTPUT:   pointer to distance matrix, has to be freed when program stops
*/
{
    long int     i, j;
    long int     **matrix;

    if((matrix = malloc(sizeof(long int) * num_node * num_node +
			sizeof(long int *) * num_node	 )) == NULL){
        fprintf(stderr,"Out of memory, exit.");
        exit(1);
    }
    for ( i = 0 ; i < num_node ; i++ ) {
        matrix[i] = (long int *)(matrix + num_node) + i*num_node;
        for ( j = 0  ; j < num_node ; j++ ) {
            matrix[i][j] = distance(i, j);
        }
    }
    return matrix;
}



long int ** compute_nn_lists( void )
/*    
      FUNCTION: computes nearest neighbor lists of depth nn for each node
      INPUT:    none
      OUTPUT:   pointer to the nearest neighbor lists
*/
{
    long int i, node, nn;
    long int *distance_vector;
    long int *help_vector;
    long int **m_nnear;
 
    TRACE ( printf("\n computing nearest neighbor lists, "); )

    nn = MAX(nn_ls,nn_ants);
    if ( nn >= num_node)
        nn = num_node - 1;
    DEBUG ( assert( num_node > nn ); )
    
    TRACE ( printf("nn = %ld ... \n",nn); ) 

    if((m_nnear = malloc(sizeof(long int) * num_node * nn
			     + num_node * sizeof(long int *))) == NULL){
        exit(EXIT_FAILURE);
    }
    distance_vector = calloc(num_node, sizeof(long int));
    help_vector = calloc(num_node, sizeof(long int));
 
    for ( node = 0 ; node < num_node ; node++ ) {  /* compute cnd-sets for all node */
        m_nnear[node] = (long int *)(m_nnear + num_node) + node * nn;

        for ( i = 0 ; i < num_node ; i++ ) {  /* Copy distances from nodes to the others */
            distance_vector[i] = instance.distance[node][i];
            help_vector[i] = i;
        }
        distance_vector[node] = LONG_MAX;  /* node is not nearest neighbour */
        distance_vector[0] = LONG_MAX;     /* depot点需要排除在 nearest neighbour之外 */
        sort2(distance_vector, help_vector, 0, num_node-1);
        for ( i = 0 ; i < nn ; i++ ) {
            m_nnear[node][i] = help_vector[i];
        }
    }
    free(distance_vector);
    free(help_vector);
    TRACE ( printf("\n    .. done\n"); )
    return m_nnear;
}


/*
 FUNCTION: compute the tour length of tour tour
 INPUT:    pointer to tour tour, tour size tour_size
 OUTPUT:   tour length of tour t
 */
long int compute_tour_length( long int *tour, long int tour_size)
{
    int      i;
    long int tour_length = 0;
  
    for ( i = 0 ; i < tour_size-1; i++ ) {
        tour_length += instance.distance[tour[i]][tour[i+1]];
    }
    return tour_length;
}

/*
 * 检查 ant vrp solution 的有效性
 * i.e. tour = [0,1,4,2,0,5,3,0] toute1 = [0,1,4,2,0] route2 = [0,5,3,0]
 */
int vrp_check_solution(const long int *tour, long int tour_size)
{
    int i;
    int * used;
    long int route_beg;    /* 单条回路起点 */
//    const long int size = num_node;

    used = calloc (num_node, sizeof(int));

    if (tour == NULL) {
        fprintf (stderr,"\n%s:error: permutation is not initialized!\n", __FUNCTION__);
        exit(1);
    }

    route_beg = 0;
    used[0] = TRUE;
    for (i = 1; i < tour_size; i++) {
        
        if (tour[i] != 0) {
        // 非depot点只能路过一次
            if (used[tour[i]]) {
                fprintf(stderr,"\n%s:error: solution vector has two times the value %ld (last position: %d)\n", __FUNCTION__, tour[i], i);
                goto error;
            } else {
                used[tour[i]] = TRUE;
            }
        }
        
        if (tour[i] == 0) {
        // 形成单条回路
            if(!vrp_check_route(tour, route_beg, i)) {
                goto error;
            }
            route_beg = i;
        }
    }

    for (i = 0; i < num_node; i++) {
        if (!used[i]) {
            fprintf(stderr,"\n%s:error: vector position %d not occupied\n", __FUNCTION__, i);
            goto error;
        }
    }

    free (used);
    return TRUE;

error:
    fprintf(stderr,"\n%s:error: solution_vector:", __FUNCTION__);
    for (i = 0; i < tour_size; i++)
        fprintf(stderr, " %ld", tour[i]);
    fprintf(stderr,"\n");
    free(used);
    return FALSE;
}

/*
 * 检查单条回路有效性
 * 注意depot点出现两次（分别出现在首尾）
 * i.e. toute1 = [0,1,4,2,0]
 */
int vrp_check_route(const long int *tour, long int rbeg, long int rend)
{
    long int load = 0;
    
    if (tour[rbeg] != 0 || tour[rend] != 0) {
        fprintf(stderr,"\n%s:error: 车辆路径没有形成一条回路\n", __FUNCTION__);
        return FALSE;
    }
    if (rend - rbeg < 2) {
        fprintf(stderr,"\n%s:error: 单条回路长度不对. rbeg=%ld, rend=%ld\n", __FUNCTION__, rbeg, rend);
        return FALSE;
    }
    for (long int i = rbeg + 1; i < rend - 1; i++) {
        load += instance.nodeptr[tour[i]].demand;
    }
    if (load > vehicle_capacity) {
        fprintf(stderr,"\n%s:error: 单条回路超过车辆最大承载量 load = %ld, capacity = %ld rbeg = %ld rend = %ld\n",
                __FUNCTION__, load, vehicle_capacity, rbeg, rend);
        return FALSE;
    }
    return TRUE;
}
