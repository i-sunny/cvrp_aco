/*********************************
 
 Ant Colony Optimization algorithms (AS, ACS, EAS, RAS, MMAS, BWAS) for CVRP
 
 Created by 孙晓奇 on 2016/10/8.
 Copyright © 2016年 xiaoqi.sxq. All rights reserved.
 
 Program's name: acovrp
 Purpose: local search routines
 
 email: sunxq1991@gmail.com
 
 *********************************/

#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <limits.h>

#include "localSearch.h"
#include "utilities.h"
#include "vrpHelper.h"
#include "io.h"

LocalSearch::LocalSearch(Problem *instance) {
    this->instance = instance;
    
    ants = instance->ants;
    n_ants = instance->n_ants;
    ls_flag = instance->ls_flag;
    num_node = instance->num_node;
    nn_list = instance->nn_list;
    nn_ls = instance->nn_ls;
    dlb_flag = instance->dlb_flag;
    distance = instance->distance;
}

/*
 FUNCTION:       manage the local search phase; apply local search to ALL ants; in
 dependence of ls_flag one of 2-opt local search is chosen.
 INPUT:          none
 OUTPUT:         none
 (SIDE)EFFECTS:  all ants of the colony have locally optimal tours
 COMMENTS:       typically, best performance is obtained by applying local search
 to all ants. It is known that some improvements (e.g. convergence
 speed towards high quality solutions) may be obtained for some
 ACO algorithms by applying local search to only some of the ants.
 Overall best performance is typcially obtained by using 3-opt.
 */
void LocalSearch::do_local_search(void)
{
    long int k;
    
    TRACE ( printf("apply local search to all ants\n"); );
    
    for ( k = 0 ; k < n_ants ; k++ ) {
        //debug
//        printf("\n--Before local search:");
//        if(check_solution(instance, ants[k].tour, ants[k].tour_size)) {
//            print_solution(instance, ants[k].tour, ants[k].tour_size);
//        }
        
        if (ls_flag) {
            two_opt_solution(ants[k].tour, ants[k].tour_size);    /* 2-opt local search */
            ants[k].tour_length = compute_tour_length(instance, ants[k].tour, ants[k].tour_size);
        }
        
        //debug
//        printf("\n--After local search:");
//        if(check_solution(instance, ants[k].tour, ants[k].tour_size)) {
//            print_solution(instance, ants[k].tour, ants[k].tour_size);
//        }
    }
}

/*
 FUNCTION:       generate a random permutation of the integers 0 .. n-1
 INPUT:          length of the array
 OUTPUT:         pointer to the random permutation
 (SIDE)EFFECTS:  the array holding the random permutation is allocated in this
 function. Don't forget to free again the memory!
 COMMENTS:       only needed by the local search procedures
 */
long int * LocalSearch::generate_random_permutation( long int n )
{
   long int  i, help, node, tot_assigned = 0;
   double    rnd;
   long int  *r;

   r = (long int *)malloc(n * sizeof(long int));

   for ( i = 0 ; i < n; i++) 
     r[i] = i;

   for ( i = 0 ; i < n ; i++ ) {
     /* find (randomly) an index for a free unit */ 
     rnd  = ran01 ( &instance->rnd_seed );
     node = (long int) (rnd  * (n - tot_assigned)); 
     assert( i + node < n );
     help = r[i];
     r[i] = r[i+node];
     r[i+node] = help;
     tot_assigned++;
   }
   return r;
}

/*
 FUNCTION:       2-opt all routes of an ant's solution.
 INPUT:          tour, an ant's solution
                 depotId
 OUTPUT:         none
 */
void LocalSearch::two_opt_solution(long int *tour, long int tour_size)
{
    long int *dlb;               /* vector containing don't look bits */
    long int *route_node_map;    /* mark for all nodes in a single route */
    long int *tour_node_pos;     /* positions of nodes in tour */

    long int route_beg = 0;
    long int *path_load = new long int[num_node-1];  /* array of single route load */
    long int p = 0;
    long int load_tmp = 0;
    
    dlb = (long int *)malloc(num_node * sizeof(long int));
    for (int i = 0 ; i < num_node; i++) {
        dlb[i] = FALSE;
    }
    
    route_node_map =  (long int *)malloc(num_node * sizeof(long int));
    for (int j = 0; j < num_node; j++) {
        route_node_map[j] = FALSE;
    }
    
    tour_node_pos =  (long int *)malloc(num_node * sizeof(long int));
    for (int j = 0; j < tour_size; j++) {
        tour_node_pos[tour[j]] = j;
    }
    
    for (int i = 1; i < tour_size; i++) {
        // 2-opt a single route from tour
        if (tour[i] == 0) {
            tour_node_pos[0] = route_beg;
            two_opt_single_route(tour, route_beg, i-1, dlb, route_node_map, tour_node_pos);
            
            for (int j = 0; j < num_node; j++) {
                route_node_map[j] = FALSE;
            }
            route_beg = i;
            
            path_load[p++] = load_tmp;
            load_tmp = 0;
        } else {
            route_node_map[tour[i]] = TRUE;
            load_tmp += instance->nodeptr[tour[i]].demand;
        }
    }
    
    DEBUG(assert(p < num_node));
    swap(tour, tour_size, path_load);
    
    delete[] path_load;
    
    free( dlb );
    free(route_node_map);
    free(tour_node_pos);
}

/*
 FUNCTION:       2-opt a single route from an ant's solution.
                 This heuristic is applied separately to each
 of the vehicle routes built by an ant.
 INPUT:          rbeg, route的起始位置
                 rend, route的结束位置(包含rend处的点, rend处不为0)
 OUTPUT:         none
 COMMENTS:       the neighbourhood is scanned in random order
 */
void LocalSearch::two_opt_single_route(long int *tour, long int rbeg, long int rend,
                          long int *dlb, long int *route_node_map, long int *tour_node_pos)
{
    long int n1, n2;                            /* nodes considered for an exchange */
    long int s_n1, s_n2;                        /* successor nodes of n1 and n2     */
    long int p_n1, p_n2;                        /* predecessor nodes of n1 and n2   */
    long int pos_n1, pos_n2;                    /* positions of nodes n1, n2        */
    long int num_route_node = rend - rbeg + 1;  /* number of nodes in a single route(depot只计一次) */
    
    long int i, j, h, l;
    long int improvement_flag, help, n_improves = 0, n_exchanges = 0;
    long int h1=0, h2=0, h3=0, h4=0;
    long int radius;             /* radius of nn-search */
    long int gain = 0;
    long int *random_vector;

    // debug
//    print_single_route(instance, tour + rbeg, num_route_node+1);

    improvement_flag = TRUE;
    random_vector = generate_random_permutation(num_route_node);

    while ( improvement_flag ) {

        improvement_flag = FALSE;

        for (l = 0 ; l < num_route_node; l++) {

            /* the neighbourhood is scanned in random order */
            pos_n1 = rbeg + random_vector[l];
            n1 = tour[pos_n1];
            if (dlb_flag && dlb[n1])
                continue;
            
            s_n1 = pos_n1 == rend ? tour[rbeg] : tour[pos_n1+1];
            radius = distance[n1][s_n1];
            /* First search for c1's nearest neighbours, use successor of c1 */
            for ( h = 0 ; h < nn_ls ; h++ ) {
                n2 = nn_list[n1][h]; /* exchange partner, determine its position */
                if (route_node_map[n2] == FALSE) {
                    /* 该点不在本route中 */
                    continue;
                }
                if (radius > distance[n1][n2] ) {
                    pos_n2 = tour_node_pos[n2];
                    s_n2 = pos_n2 == rend ? tour[rbeg] : tour[pos_n2+1];
                    gain =  - radius + distance[n1][n2] +
                            distance[s_n1][s_n2] - distance[n2][s_n2];
                    if ( gain < 0 ) {
                        h1 = n1; h2 = s_n1; h3 = n2; h4 = s_n2;
                        goto exchange2opt;
                    }
                }
                else break;
            }
            
            /* Search one for next c1's h-nearest neighbours, use predecessor c1 */
            p_n1 = pos_n1 == rbeg ? tour[rend] : tour[pos_n1-1];
            radius = distance[p_n1][n1];
            for ( h = 0 ; h < nn_ls ; h++ ) {
                n2 = nn_list[n1][h];  /* exchange partner, determine its position */
                if (route_node_map[n2] == FALSE) {
                    /* 该点不在本route中 */
                    continue;
                }
                if ( radius > distance[n1][n2] ) {
                    pos_n2 = tour_node_pos[n2];
                    p_n2 = pos_n2 == rbeg ? tour[rend] : tour[pos_n2-1];
                    
                    if ( p_n2 == n1 || p_n1 == n2)
                        continue;
                    gain =  - radius + distance[n1][n2] +
                            distance[p_n1][p_n2] - distance[p_n2][n2];
                    if ( gain < 0 ) {
                        h1 = p_n1; h2 = n1; h3 = p_n2; h4 = n2;
                        goto exchange2opt;
                    }
                }
                else break;
            }
            /* No exchange */
            dlb[n1] = TRUE;
            continue;

exchange2opt:
            n_exchanges++;
            improvement_flag = TRUE;
            dlb[h1] = FALSE; dlb[h2] = FALSE;
            dlb[h3] = FALSE; dlb[h4] = FALSE;
            /* Now perform move */
            if ( tour_node_pos[h3] < tour_node_pos[h1] ) {
                help = h1; h1 = h3; h3 = help;
                help = h2; h2 = h4; h4 = help;
            }
            /* reverse inner part from pos[h2] to pos[h3] */
            i = tour_node_pos[h2]; j = tour_node_pos[h3];
            while (i < j) {
                n1 = tour[i];
                n2 = tour[j];
                tour[i] = n2;
                tour[j] = n1;
                tour_node_pos[n1] = j;
                tour_node_pos[n2] = i;
                i++; j--;
            }
            // debug
//            printf("after ls. pid %d", instance->pid);
//            print_single_route(instance, tour + rbeg, num_route_node+1);
        }
        if ( improvement_flag ) {
            n_improves++;
        }
    }
    free( random_vector );
}


/*
 * The swap operation selects two customers at random and 
 * then swaps these two customers in their positions.
 */
void LocalSearch::swap(long int *tour, long int tour_size, long int *path_load)
{
    long int i = 0, j = 0;
    long int gain = 0;
    long int n1, p_n1, s_n1, n2, p_n2, s_n2;
    long int p1 = 0, p2 = 0;     /* path idx of node n1 and n2 */
    long int load1 = 0, load2 = 0;
    
    for (i = 1; i < tour_size; i++) {
        n1 = tour[i];
        if (n1 == 0) {
            p1++;
            continue;
        }
        p_n1 = tour[i-1];
        s_n1 = tour[i+1];
        
        p2 = p1;
        for (j = i+1; j < tour_size; j++) {
            n2 = tour[j];
            if (n2 == 0) {
                p2++;
                continue;
            }
            p_n2 = tour[j-1];
            s_n2 = tour[j+1];
            
            // node n1 and n2 not in the same route
            if (p1 != p2) {
                load1 = path_load[p1] - instance->nodeptr[n1].demand + instance->nodeptr[n2].demand;
                load2 = path_load[p2] - instance->nodeptr[n2].demand + instance->nodeptr[n1].demand;
                if (load1 > instance->vehicle_capacity || load2 > instance->vehicle_capacity) {
                    continue;
                }
            }
            
            if (j == i+1) {
               gain = -(distance[p_n1][n1] + distance[n2][s_n2]) +(distance[p_n1][n2] + distance[n1][s_n2]);
            } else {
                gain = -(distance[p_n1][n1] + distance[n1][s_n1] + distance[p_n2][n2] + distance[n2][s_n2])
                +(distance[p_n1][n2] + distance[n2][s_n1] + distance[p_n2][n1] + distance[n1][s_n2]);
            }
            if (gain < 0) {
                tour[i] = n2;
                tour[j] = n1;

                if (p1 != p2) {
                    path_load[p1] = load1;
                    path_load[p2] = load2;
                }
                
                i--;
//                check_solution(instance, tour, tour_size);
//                print_solution(instance, tour, tour_size);
                break;
            }
        }
    }
}

