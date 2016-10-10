/*********************************
 
 Ant Colony Optimization algorithms (AS, ACS, EAS, RAS, MMAS) for CVRP
 
 Created by 孙晓奇 on 2016/10/8.
 Copyright © 2016年 xiaoqi.sxq. All rights reserved.
 
 Program's name: acovrp
 Purpose: implementation of procedures for ants' behaviour
 
 email: sunxq1991@gmail.com
 
*********************************/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <limits.h>
#include <time.h>

#include "InOut.h"
#include "vrp.h"
#include "ants.h"
#include "ls.h"
#include "utilities.h"
#include "timer.h"


ant_struct *ant;
ant_struct *best_so_far_ant;
ant_struct *iteration_best_ant;

double   **pheromone;
double   **total_info;

double   *prob_of_selection;

long int n_ants;      /* number of ants */
long int nn_ants;     /* length of nearest neighbor lists for the ants'
			 solution construction */

double rho;           /* parameter for evaporation */
double alpha;         /* importance of trail */
double beta;          /* importance of heuristic evaluate */
double q_0;           /* probability of best choice in tour construction */

long int ras_flag;    /* rank-based version of ant system */
long int ras_ranks;   /* additional parameter for rank-based version
                             of ant system */
long int u_gb;            /* every u_gb iterations update with best-so-far ant */


void allocate_ants ( void )
/*    
      FUNCTION:       allocate the memory for the ant colony, the best-so-far and 
                      the iteration best ant
      INPUT:          none
      OUTPUT:         none
      (SIDE)EFFECTS:  allocation of memory for the ant colony and two ants that 
                      store intermediate tours

*/
{
    long int i;
  
    if((ant = malloc(sizeof( ant_struct ) * n_ants +
		     sizeof(ant_struct *) * n_ants	 )) == NULL){
        printf("Out of memory, exit.");
        exit(1);
    }
    for ( i = 0 ; i < n_ants ; i++ ) {
        ant[i].tour        = calloc(2*num_node-1, sizeof(long int));   // tour最长为2 * num_node - 1
        ant[i].visited     = calloc(num_node, sizeof(char));
    }

    if((best_so_far_ant = malloc(sizeof( ant_struct ) )) == NULL){
        printf("Out of memory, exit.");
        exit(1);
    }
    best_so_far_ant->tour        = calloc(2*num_node-1, sizeof(long int));
    best_so_far_ant->visited     = calloc(num_node, sizeof(char));

    if ((prob_of_selection = malloc(sizeof(double) * (nn_ants + 1))) == NULL) {
        printf("Out of memory, exit.");
        exit(1);
    }
    /* Ensures that we do not run over the last element in the random wheel.  */
    prob_of_selection[nn_ants] = HUGE_VAL;
}



long int find_best( void ) 
/*    
      FUNCTION:       find the best ant of the current iteration
      INPUT:          none
      OUTPUT:         index of struct containing the iteration best ant
      (SIDE)EFFECTS:  none
*/
{
    long int   min;
    long int   k, k_min;

    min = ant[0].tour_length;
    k_min = 0;
    for( k = 1 ; k < n_ants ; k++ ) {
        if( ant[k].tour_length < min ) {
            min = ant[k].tour_length;
            k_min = k;
        }
    }
    return k_min;
}



long int find_worst( void ) 
/*    
      FUNCTION:       find the worst ant of the current iteration
      INPUT:          none
      OUTPUT:         pointer to struct containing iteration best ant
      (SIDE)EFFECTS:  none
*/
{
    long int   max;
    long int   k, k_max;

    max = ant[0].tour_length;
    k_max = 0;
    for( k = 1 ; k < n_ants ; k++ ) {
        if( ant[k].tour_length > max ) {
            max = ant[k].tour_length;
            k_max = k;
        }
    }
    return k_max;
}



/************************************************************
 ************************************************************
Procedures for pheromone manipulation 
 ************************************************************
 ************************************************************/



void init_pheromone_trails( double initial_trail )
/*    
      FUNCTION:      initialize pheromone trails
      INPUT:         initial value of pheromone trails "initial_trail"
      OUTPUT:        none
      (SIDE)EFFECTS: pheromone matrix is reinitialized
*/
{
    long int i, j;

    TRACE ( printf(" init trails with %.15f\n",initial_trail); );

    /* Initialize pheromone trails */
    for ( i = 0 ; i < num_node ; i++ ) {
        for ( j = 0 ; j < num_node ; j++ ) {
//            if(j == i) continue;
            pheromone[i][j] = initial_trail;
        }
    }
}



void evaporation( void )
/*    
      FUNCTION:      implements the pheromone trail evaporation
      INPUT:         none
      OUTPUT:        none
      (SIDE)EFFECTS: pheromones are reduced by factor rho
*/
{ 
    long int    i, j;

    TRACE ( printf("pheromone evaporation\n"); );

    for ( i = 0 ; i < num_node ; i++ ) {
        for ( j = 0 ; j < num_node ; j++ ) {
            if(j == i) continue;
            pheromone[i][j] = (1 - rho) * pheromone[i][j];
        }
    }
}



void evaporation_nn_list( void )
/*    
      FUNCTION:      simulation of the pheromone trail evaporation
      INPUT:         none
      OUTPUT:        none
      (SIDE)EFFECTS: pheromones are reduced by factor rho
      REMARKS:       if local search is used, this evaporation procedure 
                     only considers links between a node and those nodes
		     of its candidate list
*/
{ 
    long int    i, j, help_node;

    TRACE ( printf("pheromone evaporation nn_list\n"); );

    for ( i = 0 ; i < num_node ; i++ ) {
        for ( j = 0 ; j < nn_ants ; j++ ) {
            help_node = instance.nn_list[i][j];
            pheromone[i][help_node] = (1 - rho) * pheromone[i][help_node];
        }
    }
}



void global_update_pheromone( ant_struct *a )
/*    
      FUNCTION:      reinforces edges used in ant k's solution
      INPUT:         pointer to ant that updates the pheromone trail 
      OUTPUT:        none
      (SIDE)EFFECTS: pheromones of arcs in ant k's tour are increased
*/
{  
    long int i, j, h;
    double   d_tau;

    TRACE ( printf("global pheromone update\n"); );

    d_tau = 1.0 / (double) a->tour_length;
    for ( i = 0 ; i < a->tour_size ; i++ ) {
        j = a->tour[i];
        h = a->tour[i+1];
        pheromone[j][h] += d_tau;
    }
}



void global_update_pheromone_weighted( ant_struct *a, long int weight )
/*    
      FUNCTION:      reinforces edges of the ant's tour with weight "weight"
      INPUT:         pointer to ant that updates pheromones and its weight  
      OUTPUT:        none
      (SIDE)EFFECTS: pheromones of arcs in the ant's tour are increased
*/
{  
    long int      i, j, h;
    double        d_tau;

    TRACE ( printf("global pheromone update weighted\n"); );

    d_tau = (double) weight / (double) a->tour_length;
    for ( i = 0 ; i < a->tour_size ; i++ ) {
        j = a->tour[i];
        h = a->tour[i+1];
        pheromone[j][h] += d_tau;
    }       
}



void compute_total_information( void )
/*    
      FUNCTION: calculates heuristic info times pheromone for each arc
      INPUT:    none  
      OUTPUT:   none
*/
{
    long int     i, j;

    TRACE ( printf("compute total information\n"); );

    for ( i = 0 ; i < num_node ; i++ ) {
        for ( j = 0 ; j < num_node; j++ ) {
            if(j == i) continue;
            total_info[i][j] = pow(pheromone[i][j], alpha) * pow(HEURISTIC(i,j), beta);
        }
    }
}



void compute_nn_list_total_information( void )
/*    
      FUNCTION: calculates heuristic info times pheromone for arcs in nn_list
      INPUT:    none  
      OUTPUT:   none
*/
{ 
    long int    i, j, h;

    TRACE ( printf("compute total information nn_list\n"); );

    for ( i = 0 ; i < num_node ; i++ ) {
        for ( j = 0 ; j < nn_ants ; j++ ) {
            h = instance.nn_list[i][j];
            total_info[i][h] = pow(pheromone[i][h], alpha) * pow(HEURISTIC(i,h),beta);
        }
    }
}



/****************************************************************
 ****************************************************************
Procedures implementing solution construction and related things
 ****************************************************************
 ****************************************************************/



void ant_empty_memory( ant_struct *a ) 
/*    
      FUNCTION:       empty the ants's memory regarding visited nodes
      INPUT:          ant identifier
      OUTPUT:         none
      (SIDE)EFFECTS:  vector of visited nodes is reinitialized to FALSE
*/
{
    long int   i;

    for( i = 0 ; i < num_node ; i++ ) {
        a->visited[i]=FALSE;
    }
    a->tour_size = 0;
}



void init_ant_place( ant_struct *a , long int phase)
/*    
      FUNCTION:      place an ant on the single depot
      INPUT:         pointer to ant and the number of construction steps 
      OUTPUT:        none
      (SIDE)EFFECT:  ant is put on the depot
*/
{
    a->tour[phase] = 0;
    a->visited[0] = TRUE;
}



long int choose_best_next( ant_struct *a, long int phase )
/*    
      FUNCTION:      chooses for an ant as the next node the one with
                     maximal value of heuristic information times pheromone 
      INPUT:         pointer to ant and the construction step
      OUTPUT:        none 
      (SIDE)EFFECT:  ant moves to the chosen node
*/
{ 
    long int node, current_node, next_node;
    double   value_best;

    next_node = num_node;
    DEBUG( assert ( phase > 0 && phase < n ); );
    current_node = a->tour[phase-1];
    value_best = -1.;             /* values in total matrix are always >= 0.0 */    
    for ( node = 0 ; node < num_node ; node++ ) {
        if ( a->visited[node] ) {
            ; /* node already visited, do nothing */
        } else if(demand_meet_node_map[node] == FALSE) {
            ;  /* 该点不满足要求 */
        } else {
            if ( total_info[current_node][node] > value_best ) {
                next_node = node;
                value_best = total_info[current_node][node];
            }
        }
    }
    DEBUG( assert ( 0 <= next_node && next_node < n); );
    DEBUG( assert ( value_best > 0.0 ); )
    DEBUG( assert ( a->visited[next_node] == FALSE ); )
    a->tour[phase] = next_node;
    a->visited[next_node] = TRUE;
    
    return next_node;
}



long int neighbour_choose_best_next( ant_struct *a, long int phase )
/*    
      FUNCTION:      chooses for an ant as the next node the one with
                     maximal value of heuristic information times pheromone 
      INPUT:         pointer to ant and the construction step "phase" 
      OUTPUT:        none 
      (SIDE)EFFECT:  ant moves to the chosen node
*/
{ 
    long int i, current_node, next_node, help_node;
    double   value_best, help;
  
    next_node = num_node;
    DEBUG( assert ( phase > 0 && phase < num_node ); );
    current_node = a->tour[phase-1];
    DEBUG ( assert ( 0 <= current_node && current_node < num_node ); )
    value_best = -1.;             /* values in total matix are always >= 0.0 */    
    for ( i = 0 ; i < nn_ants ; i++ ) {
        help_node = instance.nn_list[current_node][i];
        if ( a->visited[help_node] ) {
            ;   /* node already visited, do nothing */
        } else if(demand_meet_node_map[help_node] == FALSE) {
            ;  /* 该点不满足要求 */
        } else {
            help = total_info[current_node][help_node];
            if ( help > value_best ) {
                value_best = help;
                next_node = help_node;
            }
        }
    }
    if ( next_node == num_node )
	/* all nodes in nearest neighbor list were already visited */
        choose_best_next( a, phase );
    else {
        DEBUG( assert ( 0 <= next_node && next_node < num_node); )
        DEBUG( assert ( value_best > 0.0 ); )
        DEBUG( assert ( a->visited[next_node] == FALSE ); )
        a->tour[phase] = next_node;
        a->visited[next_node] = TRUE;
    }
    return next_node;
}



void choose_closest_next( ant_struct *a, long int phase )
/*    
      FUNCTION:      Chooses for an ant the closest node as the next one 
      INPUT:         pointer to ant and the construction step "phase" 
      OUTPUT:        none 
      (SIDE)EFFECT:  ant moves to the chosen node
*/
{ 
    long int node, current_node, next_node, min_distance;
  
    next_node = num_node;
    DEBUG( assert ( phase > 0 && phase < num_node ); );
    current_node = a->tour[phase-1];
    min_distance = INFTY;             /* Search shortest edge */    
    for ( node = 0 ; node < num_node ; node++ ) {
        if ( a->visited[node] ) {
            ; /* node already visited */
        } else if(demand_meet_node_map[node] == FALSE) {
            ;  /* 该点不满足要求 */
        } else {
            if ( instance.distance[current_node][node] < min_distance) {
                next_node = node;
                min_distance = instance.distance[current_node][node];
            }
        } 
    }
    DEBUG( assert ( 0 <= next_node && next_node < num_node); );
    a->tour[phase] = next_node;
    a->visited[next_node] = TRUE;
}


long int neighbour_choose_and_move_to_next(ant_struct *a, long int phase)
/*    
     FUNCTION:      Choose for an ant probabilistically a next node among all
     unvisited and possible nodes in the current node's candidate list.
     If this is not possible, choose the closest next
     INPUT:         pointer to ant the construction step "phase"
     OUTPUT:        none
     (SIDE)EFFECT:  ant moves to the chosen node
*/
{
    long int i, help;
    long int current_node, neighbour_node;
    double   rnd, partial_sum = 0., sum_prob = 0.0;
    /*  double   *prob_of_selection; */ /* stores the selection probabilities 
	of the nearest neighbor nodes */
    double   *prob_ptr;



    if ( (q_0 > 0.0) && (ran01( &seed ) < q_0)  ) {
        /* with a probability q_0 make the best possible choice
           according to pheromone trails and heuristic information */
        /* we first check whether q_0 > 0.0, to avoid the very common case
           of q_0 = 0.0 to have to compute a random number, which is
           expensive computationally */
        return neighbour_choose_best_next(a, phase);
    }

    prob_ptr = prob_of_selection;

    current_node = a->tour[phase-1]; /* current_node node of ant k */
    DEBUG( assert ( current_node >= 0 && current_node < num_node ); )
    for ( i = 0 ; i < nn_ants ; i++ ) {
        neighbour_node = instance.nn_list[current_node][i];
        if ( a->visited[neighbour_node] ) {
            prob_ptr[i] = 0.0;   /* node already visited */
        } else if(demand_meet_node_map[neighbour_node] == FALSE) {
            prob_ptr[i] = 0.0;  /* 该点不满足要求 */
        } else {
            DEBUG( assert ( neighbour_node >= 0 && neighbour_node < num_node ); )
            prob_ptr[i] = total_info[current_node][neighbour_node];
            sum_prob += prob_ptr[i];
        } 
    }

    if (sum_prob <= 0.0) {
        /* All nodes from the candidate set are tabu */
        return choose_best_next( a, phase);
    } else {
        /* at least one neighbor is eligible, chose one according to the
           selection probabilities */
        rnd = ran01( &seed );
        rnd *= sum_prob;
        i = 0;
        partial_sum = prob_ptr[i];
        /* This loop always stops because prob_ptr[nn_ants] == HUGE_VAL  */
        while (partial_sum <= rnd) {
            i++;
            partial_sum += prob_ptr[i];
        }
        /* This may very rarely happen because of rounding if rnd is
           close to 1.  */
        if (i == nn_ants) {
            return neighbour_choose_best_next(a, phase);
        }
        DEBUG( assert ( 0 <= i && i < nn_ants); );
        DEBUG( assert ( prob_ptr[i] >= 0.0); );
        help = instance.nn_list[current_node][i];
        DEBUG( assert ( help >= 0 && help < n ); )
        DEBUG( assert ( a->visited[help] == FALSE ); )
        a->tour[phase] = help; /* instance.nn_list[current_node][i]; */
        a->visited[help] = TRUE;
        
        return help;
    }
}


/**************************************************************************
 **************************************************************************
Procedures specific to the ant's tour manipulation other than construction
***************************************************************************
 **************************************************************************/


/*
 FUNCTION:       copy solution from ant a1 into ant a2
 INPUT:          pointers to the two ants a1 and a2
 OUTPUT:         none
 (SIDE)EFFECTS:  a2 is copy of a1
 */
void copy_solution_from_to(ant_struct *a1, ant_struct *a2)
{
    int   i;
  
    a2->tour_length = a1->tour_length;
    a2->tour_size = a1->tour_size;
    for ( i = 0 ; i < a1->tour_size ; i++ ) {
        a2->tour[i] = a1->tour[i];
    }
}


