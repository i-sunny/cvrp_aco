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

#include "antColony.h"
#include "simulatedAnnealing.h"
#include "utilities.h"
#include "vrpHelper.h"
#include "problem.h"
#include "io.h"
#include "timer.h"

AntColony::AntColony(Problem *instance)
{
    this->instance = instance;
    local_search = new LocalSearch(instance);
    
    ants = instance->ants;
    best_so_far_ant = instance->best_so_far_ant;
    
    distance = instance->distance;
    prob_of_selection = instance->prob_of_selection;
    pheromone = instance->pheromone;
    total_info = instance->total_info;
    
    num_node = instance->num_node;
    n_ants = instance->n_ants;
    nn_ants = instance->nn_ants;
    nn_list = instance->nn_list;
    
    nodeptr = instance->nodeptr;
    vehicle_capacity = instance->vehicle_capacity;
    ls_flag = instance->ls_flag;
    
}

AntColony::~AntColony(){
    delete local_search;
}


/****************************************************************
 ****************************************************************
Procedures implementing init/exit aco and aco steps:
 1) solution construction
 2) pheromone update
 ****************************************************************
 ****************************************************************/

/*
 FUNCTION: initilialize variables appropriately when starting a trial
 INPUT:    trial number
 OUTPUT:   none
 COMMENTS: none
 */
void AntColony::init_aco()
{
    instance->best_so_far_time = elapsed_time(VIRTUAL);
    
    /* Initialize variables concerning statistics etc. */
    instance->iteration   = 0;
    best_so_far_ant->tour_length = INFTY;
    
    instance->iter_stagnate_cnt = 0;
    instance->best_stagnate_cnt = 0;
    
    /* Initialize the Pheromone trails */
    /* in the original papers on Ant System, Elitist Ant System, and
     Rank-based Ant System it is not exactly defined what the
     initial value of the pheromones is. Here we set it to some
     small constant, analogously as done in MAX-MIN Ant System.
     */
    double trail_0 = 0.5;
    init_pheromone_trails(trail_0);
    compute_total_information();
    
    // 第一次迭代用于设置一个合适的 pheromone init trail
    construct_solutions();
    if (ls_flag) {
        local_search->do_local_search();
    }
    update_statistics();
    trail_0 =  1.0 / ((rho) * best_so_far_ant->tour_length);
    init_pheromone_trails(trail_0);
    instance->iteration++;
    
    /* Calculate combined information pheromone times heuristic information */
    compute_total_information();
    
}

/*
 * exit of Ant Colony Optimization
 */
void AntColony::exit_aco()
{
    
}

/*
 * 蚁群算法单次迭代的执行
 */
void AntColony::run_aco_iteration()
{
    construct_solutions();
    
    if (ls_flag) {
        local_search->do_local_search();
    }
    
    update_statistics();
    
    pheromone_trail_update();
    
    if (sa_flag) {
        if ((instance == 0 && instance->best_stagnate_cnt >= instance->num_node)
            || (instance != 0 && instance->best_stagnate_cnt >= 30))
        {   // 2 * instance->num_node
            SimulatedAnnealing *annealer = new SimulatedAnnealing(instance, this, 5.0, 0.97, MAX(instance->num_node * 4, 250), 50);
            annealer->run();
            instance->best_stagnate_cnt = 0;
            
            delete annealer;
        }
    }
    
//    print_probabilities(instance);
}

/*
 FUNCTION:       manage the solution construction phase
 INPUT:          none
 OUTPUT:         none
 (SIDE)EFFECTS:  when finished, all ants of the colony have constructed a solution
 */
void AntColony::construct_solutions( void )
{
    long int k;
    
    TRACE ( printf("construct solutions for all ants\n"); );
    
    for(k = 0; k < n_ants; k++) {
        construct_ant_solution(&ants[k]);
    }
}

/*
 FUNCTION:       construct the solution for this ant
 INPUT:          this ant id
 OUTPUT:         none
 (SIDE)EFFECTS:  when finished, this ant has constructed a solution
 */

void AntColony::construct_ant_solution(AntStruct *ant)
{
    long int visited_node_cnt = 0;   /* count of visited node by this ant */
    
    long int path_load;          /* 单次从depot出发的送货量 */
    long int next_node, current_node;
    long int i, demand_meet_cnt, step;
    double path_distance;
    
    /* Mark all nodes as unvisited */
    ant_empty_memory(ant);
    
    path_load = 0;
    path_distance = 0;
    step = 0;
    init_ant_place(ant, step);
    
    while (visited_node_cnt < num_node - 1) {
        current_node = ant->tour[step];
        step++;
        
        /* 查看所有可以派送的点 */
        demand_meet_cnt = 0;
        for (i = 0; i < num_node; i++) {
            ant->demand_meet_node[i] = FALSE;
        }
        for(i = 0; i < num_node; i++) {
            if (ant->visited[i] == FALSE
                && path_load + nodeptr[i].demand <= vehicle_capacity
                && path_distance + (distance[current_node][i] + instance->service_time) + distance[i][0] <= instance->max_distance) {
                ant->demand_meet_node[i] = TRUE;
                demand_meet_cnt++;
            }
        }
        
        /*
         1)如果没有可行的配送点,则蚂蚁回到depot，重新开始新的路径
         2）否则，选择下一个配送点
         */
        if (demand_meet_cnt == 0) {
            path_load = 0;
            path_distance = 0;
            init_ant_place(ant, step);
        } else {
            next_node = neighbour_choose_and_move_to_next(ant, step);
            path_load += nodeptr[next_node].demand;
            path_distance += distance[current_node][next_node] + instance->service_time;
            visited_node_cnt++;
        }
    }
    
    // 最后回到depot
    step++;
    ant->tour[step] = ant->tour[0];
    ant->tour_size = step + 1;
    ant->tour_length = compute_tour_length(instance, ant->tour, ant->tour_size);
    
    // debug
    DEBUG(check_solution(instance, ant->tour, ant->tour_size));
}


/*
 FUNCTION:       manage global pheromone deposit for Rank-based Ant System
 INPUT:          none
 OUTPUT:         none
 (SIDE)EFFECTS:  the ras_ranks-1 best ants plus the best-so-far ant deposit pheromone
 on matrix "pheromone"
 COMMENTS:       this procedure could be implemented slightly faster, but it is
 anyway not critical w.r.t. CPU time given that ras_ranks is
 typically very small.
 */
void AntColony::ras_update( void )
{
    long int i, k, b, target;
    double *help_b;
    
    TRACE ( printf("Rank-based Ant System pheromone deposit\n"); );
    
    help_b = (double *)malloc( n_ants  * sizeof(double) );
    for ( k = 0 ; k < n_ants ; k++ )
        help_b[k] = ants[k].tour_length;
    
    for ( i = 0 ; i < ras_ranks-1 ; i++ ) {
        b = help_b[0]; target = 0;
        for ( k = 0 ; k < n_ants ; k++ ) {
            if ( help_b[k] < b ) {
                b = help_b[k];
                target = k;
            }
        }
        help_b[target] = LONG_MAX;
        global_update_pheromone_weighted(&ants[target], ras_ranks-i-1);
    }
    global_update_pheromone_weighted(best_so_far_ant, ras_ranks);
    free ( help_b );
}

/*
 * 蚁群停滞时，加入扰动跳出局部最优解
 */
void AntColony::pheromone_disturbance(void)
{
//    print_pheromone(instance);
    
    printf("pid %d start pheromone disturbance: iter %ld, best_stagnate %ld, iter_stagnate %ld\n",
           instance->pid, instance->iteration, instance->best_stagnate_cnt, instance->iter_stagnate_cnt);
    
    long int i, j;
    double sum_pheromone = 0, mean_pheromone;
    double delta = 0.7;
    
    for (i = 0; i < num_node; i++) {
        for (j = 0; j < num_node; j++) {
            sum_pheromone += pheromone[i][j];
        }
    }
    
    mean_pheromone = sum_pheromone / (num_node * num_node);
    
    for (i = 0; i < num_node; i++) {
        for (j = 0; j < num_node; j++) {
            pheromone[i][j] = (1- delta) * mean_pheromone + delta * pheromone[i][j];
        }
    }
    
//    print_pheromone(instance);
}


/*
 FUNCTION:       manage global pheromone trail update for the ACO algorithms
 INPUT:          none
 OUTPUT:         none
 (SIDE)EFFECTS:  pheromone trails are evaporated and pheromones are deposited
 according to the rules defined by the various ACO algorithms.
 */
void AntColony::pheromone_trail_update( void )
{
    /* Simulate the pheromone evaporation of all pheromones; this is not necessary
     for ACS (see also ACO Book) */
    if (ls_flag) {
        evaporation_nn_list();
        /* evaporate only pheromones on arcs of candidate list to make the
         pheromone evaporation faster for being able to tackle large TSP
         instances. For MMAS additionally check lower pheromone trail limits.
         */
    } else {
        /* if no local search is used, evaporate all pheromone trails */
        evaporation();
    }
    
    if (instance->iter_stagnate_cnt >= 5 ||
        instance->best_stagnate_cnt >= HUGE_VAL)
    {
        pheromone_disturbance();
        instance->iter_stagnate_cnt -= 2;
    } else {
        /* Next, apply the pheromone deposit for the various ACO algorithms */
        ras_update();
    }
    /* Compute combined information pheromone times heuristic info after
     the pheromone update for all ACO algorithms except ACS; in the ACS case
     this is already done in the pheromone update procedures of ACS */
    if (ls_flag) {
        compute_nn_list_total_information();
    } else {
        compute_total_information();
    }
}

/*
 FUNCTION:       manage some statistical information about the solution, especially
 if a new best solution (best-so-far or restart-best) is found
 INPUT:          none
 OUTPUT:         none
 (SIDE)EFFECTS:  best-so-far ant may be updated
 */
void AntColony::update_statistics()
{
    /* 本次迭代中结果最优的蚂蚁 */
    instance->iteration_best_ant = &ants[find_best()];
    
    if (instance->pid == 0) {
        write_iter_report(instance);
    }

    if (instance->iteration_best_ant->tour_length - best_so_far_ant->tour_length < -EPSILON) {
        // 获得更优解
        instance->best_stagnate_cnt = 0;
        
        instance->best_so_far_time = elapsed_time( VIRTUAL );
        copy_solution_from_to(instance->iteration_best_ant, best_so_far_ant );
        
        instance->best_solution_iter = instance->iteration;
        if (instance->pid == 0) {
            write_best_so_far_report(instance);
        }
    } else {
        instance->best_stagnate_cnt++;
        if (fabs(instance->last_iter_solution - instance->iteration_best_ant->tour_length) < EPSILON) {
            instance->iter_stagnate_cnt++;
        } else {
            instance->iter_stagnate_cnt = 0;
        }
    }
    
    instance->last_iter_solution = instance->iteration_best_ant->tour_length;
}


/****************************************************************
 ****************************************************************
 Procedures for finding best/worst tour in an iteration
 ****************************************************************
 ****************************************************************/

long int AntColony::find_best(void)
/*    
      FUNCTION:       find the best ant of the current iteration
      INPUT:          none
      OUTPUT:         index of struct containing the iteration best ant
      (SIDE)EFFECTS:  none
*/
{
    double   min;
    long int   k, k_min;

    min = ants[0].tour_length;
    k_min = 0;
    for( k = 1 ; k < n_ants ; k++ ) {
        if(ants[k].tour_length < min ) {
            min = ants[k].tour_length;
            k_min = k;
        }
    }
    return k_min;
}



long int AntColony::find_worst(void)
/*    
      FUNCTION:       find the worst ant of the current iteration
      INPUT:          none
      OUTPUT:         pointer to struct containing iteration best ant
      (SIDE)EFFECTS:  none
*/
{
    double   max;
    long int   k, k_max;

    max = ants[0].tour_length;
    k_max = 0;
    for( k = 1 ; k < n_ants ; k++ ) {
        if( ants[k].tour_length > max ) {
            max = ants[k].tour_length;
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



void AntColony::init_pheromone_trails( double initial_trail )
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



void AntColony::evaporation( void )
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



void AntColony::evaporation_nn_list( void )
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
            help_node = nn_list[i][j];
            pheromone[i][help_node] = (1 - rho) * pheromone[i][help_node];
        }
    }
}



void AntColony::global_update_pheromone( AntStruct *a )
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
    for ( i = 0 ; i < a->tour_size-1 ; i++ ) {
        j = a->tour[i];
        h = a->tour[i+1];
        pheromone[j][h] += d_tau;
    }
}



void AntColony::global_update_pheromone_weighted( AntStruct *a, long int weight )
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
    for ( i = 0 ; i < a->tour_size-1 ; i++ ) {
        j = a->tour[i];
        h = a->tour[i+1];
        pheromone[j][h] += d_tau;
    }       
}



void AntColony::compute_total_information( void )
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



void AntColony::compute_nn_list_total_information( void )
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
            h = nn_list[i][j];
            total_info[i][h] = pow(pheromone[i][h], alpha) * pow(HEURISTIC(i,h),beta);
        }
    }
}



/****************************************************************
 ****************************************************************
Procedures implementing solution construction and related things
 ****************************************************************
 ****************************************************************/



void AntColony::ant_empty_memory( AntStruct *a )
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



void AntColony::init_ant_place( AntStruct *a , long int phase)
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



long int AntColony::choose_best_next( AntStruct *a, long int phase )
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
    DEBUG( assert ( phase > 0 && phase < 2*num_node-2 ); );
    current_node = a->tour[phase-1];
    value_best = -1.;             /* values in total matrix are always >= 0.0 */    
    for ( node = 0 ; node < num_node ; node++ ) {
        if ( a->visited[node] ) {
            ; /* node already visited, do nothing */
        } else if(a->demand_meet_node[node] == FALSE) {
            ;  /* 该点不满足要求 */
        } else {
            if ( total_info[current_node][node] > value_best ) {
                next_node = node;
                value_best = total_info[current_node][node];
            }
        }
    }
    DEBUG( assert ( 0 <= next_node && next_node < num_node));
    DEBUG( assert ( value_best > 0.0 ); )
    DEBUG( assert ( a->visited[next_node] == FALSE ); )
    a->tour[phase] = next_node;
    a->visited[next_node] = TRUE;
    
    return next_node;
}



long int AntColony::neighbour_choose_best_next( AntStruct *a, long int phase )
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
    DEBUG( assert ( phase > 0 && phase < 2*num_node-2 ); );
    current_node = a->tour[phase-1];
    DEBUG ( assert ( 0 <= current_node && current_node < num_node ); )
    value_best = -1.;             /* values in total matix are always >= 0.0 */    
    for ( i = 0 ; i < nn_ants ; i++ ) {
        help_node = nn_list[current_node][i];
        if ( a->visited[help_node] ) {
            ;   /* node already visited, do nothing */
        } else if(a->demand_meet_node[help_node] == FALSE) {
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



void AntColony::choose_closest_next( AntStruct *a, long int phase )
/*    
      FUNCTION:      Chooses for an ant the closest node as the next one 
      INPUT:         pointer to ant and the construction step "phase" 
      OUTPUT:        none 
      (SIDE)EFFECT:  ant moves to the chosen node
*/
{ 
    long int node, current_node, next_node;
    double min_distance;
  
    next_node = num_node;
    DEBUG( assert ( phase > 0 && phase < 2*num_node-2 ); );
    current_node = a->tour[phase-1];
    min_distance = INFTY;             /* Search shortest edge */    
    for ( node = 0 ; node < num_node ; node++ ) {
        if ( a->visited[node] ) {
            ; /* node already visited */
        } else if(a->demand_meet_node[node] == FALSE) {
            ;  /* 该点不满足要求 */
        } else {
            if ( distance[current_node][node] < min_distance) {
                next_node = node;
                min_distance = distance[current_node][node];
            }
        } 
    }
    DEBUG( assert ( 0 <= next_node && next_node < num_node); );
    a->tour[phase] = next_node;
    a->visited[next_node] = TRUE;
}


long int AntColony::neighbour_choose_and_move_to_next(AntStruct *a, long int phase)
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

    prob_ptr = prob_of_selection;

    current_node = a->tour[phase-1]; /* current_node node of ant k */
    DEBUG( assert ( current_node >= 0 && current_node < num_node ); )
    for ( i = 0 ; i < nn_ants ; i++ ) {
        neighbour_node = nn_list[current_node][i];
        if ( a->visited[neighbour_node] ) {
            prob_ptr[i] = 0.0;   /* node already visited */
        } else if(a->demand_meet_node[neighbour_node] == FALSE) {
            prob_ptr[i] = 0.0;  /* 该点不满足要求 */
        } else {
            DEBUG( assert ( neighbour_node >= 0 && neighbour_node < num_node ); )
            prob_ptr[i] = total_info[current_node][neighbour_node];
            sum_prob += prob_ptr[i];
        } 
    }

    if (sum_prob <= 0.0) {
        /* All nodes from the candidate set are tabu */
        return choose_best_next(a, phase);
    } else {
        /* at least one neighbor is eligible, chose one according to the
           selection probabilities */
        rnd = ran01( &instance->rnd_seed );
        rnd *= sum_prob;
        DEBUG(assert ( rnd >= 0 && rnd <= sum_prob );)
        
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
        DEBUG( assert ( prob_ptr[i] > 0.0); );
        help = nn_list[current_node][i];
        DEBUG(assert ( help >= 0 && help < num_node );)
        DEBUG( assert ( a->visited[help] == FALSE ); )
        a->tour[phase] = help; /* nn_list[current_node][i]; */
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
void AntColony::copy_solution_from_to(AntStruct *a1, AntStruct *a2)
{
    int   i;
  
    a2->tour_length = a1->tour_length;
    a2->tour_size = a1->tour_size;
    for ( i = 0 ; i < a1->tour_size ; i++ ) {
        a2->tour[i] = a1->tour[i];
    }
}


