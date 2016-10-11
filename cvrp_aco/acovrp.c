/*********************************
 
 Ant Colony Optimization algorithms (AS, ACS, EAS, RAS, MMAS) for CVRP
 
 Created by 孙晓奇 on 2016/10/8.
 Copyright © 2016年 xiaoqi.sxq. All rights reserved.
 
 Program's name: acovrp
 Purpose: main routines and control for the ACO algorithms
 
 email: sunxq1991@gmail.com
 
 *********************************/


#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#include "ants.h"
#include "utilities.h"
#include "InOut.h"
#include "vrp.h"
#include "timer.h"
#include "ls.h"


/*
 FUNCTION:       checks whether termination condition is met
 INPUT:          none
 OUTPUT:         0 if condition is not met, number neq 0 otherwise
 (SIDE)EFFECTS:  none
 */
long int termination_condition( void )
{
  return ((iteration >= max_iteration) ||
          (elapsed_time( VIRTUAL ) >= max_time) ||
          (best_so_far_ant->tour_length <= optimal));
}

/*
 FUNCTION:       construct the solution for this ant
 INPUT:          this ant id
 OUTPUT:         none
 (SIDE)EFFECTS:  when finished, this ant has constructed a solution
 */

void construct_solution(long int ant_id)
{
    ant_struct *antk = &ant[ant_id];
    long int visited_node_cnt = 0;   /* count of visited node by this ant */
    
    long int path_load;          /* 单次从depot出发的送货量 */
    long int next_node;
    long int i, demand_meet_cnt, step;
    
    /* Mark all nodes as unvisited */
    ant_empty_memory(antk);
    
    path_load = 0;
    step = 0;
    init_ant_place(antk, step);
    
    while (visited_node_cnt < num_node - 1) {
        
        step++;
        /* 查看所有可以派送的点 */
        demand_meet_cnt = 0;
        for (i = 0; i < num_node; i++) {
            demand_meet_node_map[i] = FALSE;
        }
        for(i = 0; i < num_node; i++) {
            if (antk->visited[i] == FALSE && path_load + instance.nodeptr[i].demand <= vehicle_capacity) {
                demand_meet_node_map[i] = TRUE;
                demand_meet_cnt++;
            }
        }
        
        /*
         1)如果没有可行的配送点,则蚂蚁回到depot，重新开始新的路径
         2）否则，选择下一个配送点
         */
        if (demand_meet_cnt == 0) {
            path_load = 0;
            init_ant_place(antk, step);
        } else {
            next_node = neighbour_choose_and_move_to_next(antk, step);
            path_load += instance.nodeptr[next_node].demand;
            visited_node_cnt++;
        }
    }
    
    // 最后回到depot
    step++;
    antk->tour[step] = antk->tour[0];
    antk->tour_size = step + 1;
    antk->tour_length = compute_tour_length(antk->tour, antk->tour_size);
    
    // debug
//    if(vrp_check_solution(antk->tour, antk->tour_size)) {
//        print_solution(ant_id, antk->tour, antk->tour_size);
//    }
}


/*
 FUNCTION:       manage the solution construction phase
 INPUT:          none
 OUTPUT:         none
 (SIDE)EFFECTS:  when finished, all ants of the colony have constructed a solution
 */
void construct_solutions( void )
{
    long int k;

    TRACE ( printf("construct solutions for all ants\n"); );

    for(k = 0; k < n_ants; k++) {
        construct_solution(k);
    }
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
void local_search( void )
{
    long int k;

    TRACE ( printf("apply local search to all ants\n"); );

    for ( k = 0 ; k < n_ants ; k++ ) {
        switch (ls_flag) {
            case 1:
                two_opt_solution(ant[k].tour, ant[k].tour_size);    /* 2-opt local search */
                break;
            default:
                fprintf(stderr,"type of local search procedure not correctly specified\n");
                exit(1);
        }
        ant[k].tour_length = compute_tour_length(ant[k].tour, ant[k].tour_size);
        
        // debug
//        printf("\n--After local search:--\n");
//        if(vrp_check_solution(ant[k].tour, ant[k].tour_size)) {
//            print_solution(k, ant[k].tour, ant[k].tour_size);
//        }
        
        if (termination_condition()) return;
    }
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
void ras_update( void )
{
    long int i, k, b, target;
    long int *help_b;

    TRACE ( printf("Rank-based Ant System pheromone deposit\n"); );

    help_b = malloc( n_ants  * sizeof(long int) );
    for ( k = 0 ; k < n_ants ; k++ )
        help_b[k] = ant[k].tour_length;

    for ( i = 0 ; i < ras_ranks-1 ; i++ ) {
        b = help_b[0]; target = 0;
        for ( k = 0 ; k < n_ants ; k++ ) {
            if ( help_b[k] < b ) {
                b = help_b[k];
                target = k;
            }
        }
        help_b[target] = LONG_MAX;
        global_update_pheromone_weighted( &ant[target], ras_ranks-i-1 );
    }
    global_update_pheromone_weighted( best_so_far_ant, ras_ranks );
    free ( help_b );
}

/*
 FUNCTION:       manage global pheromone trail update for the ACO algorithms
 INPUT:          none
 OUTPUT:         none
 (SIDE)EFFECTS:  pheromone trails are evaporated and pheromones are deposited
 according to the rules defined by the various ACO algorithms.
 */
void pheromone_trail_update( void )
{
    /* Simulate the pheromone evaporation of all pheromones; this is not necessary 
       for ACS (see also ACO Book) */
    if ( ls_flag ) {
        evaporation_nn_list();
        /* evaporate only pheromones on arcs of candidate list to make the 
           pheromone evaporation faster for being able to tackle large TSP 
           instances. For MMAS additionally check lower pheromone trail limits.
        */
    } else {
        /* if no local search is used, evaporate all pheromone trails */
        evaporation();
    }

    /* Next, apply the pheromone deposit for the various ACO algorithms */
    ras_update();

  /* Compute combined information pheromone times heuristic info after
     the pheromone update for all ACO algorithms except ACS; in the ACS case 
     this is already done in the pheromone update procedures of ACS */
    if ( ls_flag ) {
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
void update_statistics( void )
{
    
    /* 本次迭代中结果最优的蚂蚁 */
    iteration_best_ant = &ant[find_best()];
    write_iter_report();
    
    if (iteration_best_ant->tour_length < best_so_far_ant->tour_length) {
        
        best_so_far_time = elapsed_time( VIRTUAL );
        copy_solution_from_to(iteration_best_ant, best_so_far_ant );
        
        best_solution_iter = iteration;
        write_best_so_far_report();
    }
}

void exit_try()
{
    
}

void init_try()
/*
 FUNCTION: initilialize variables appropriately when starting a trial
 INPUT:    trial number
 OUTPUT:   none
 COMMENTS: none
 */
{
    
    TRACE ( printf("INITIALIZE TRIAL\n"); );
    
    start_timers();
    best_so_far_time = elapsed_time( VIRTUAL );
    
    
    /* Initialize variables concerning statistics etc. */
    iteration   = 0;
    best_so_far_ant->tour_length = INFTY;
    
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
    if (ls_flag > 0) {
        local_search();
    }
    update_statistics();
    trail_0 =  1.0 / ((rho) * iteration_best_ant->tour_length);
    init_pheromone_trails(trail_0);
    iteration++;
    
    /* Calculate combined information pheromone times heuristic information */
    compute_total_information();
    
}

/*
 * 蚁群算法单次迭代的执行
 */
void aco_iteration_runner(void)
{
    
    construct_solutions();
    
    if (ls_flag > 0) {
        local_search();
    }
    
    update_statistics();
    
    pheromone_trail_update();
}

/* --- main program ------------------------------------------------------ */
/*
 FUNCTION:       main control for running the ACO algorithms
 INPUT:          none
 OUTPUT:         none
 (SIDE)EFFECTS:  none
 COMMENTS:       this function controls the run of "max_tries" independent trials
 
 */
int main(int argc, char *argv[]) {

    long int i;

    start_timers();

    init_program(argc, argv);

    best_so_far_time = elapsed_time( VIRTUAL );
    printf("Initialization took %.10f seconds\n",best_so_far_time);

	init_try();

	while ( !termination_condition() ) {
        
        aco_iteration_runner();
	    iteration++;
	}
	exit_try();
    
    exit_program();

    free( instance.distance );
    free( instance.nn_list );
    free( pheromone );
    free( total_info );
    for ( i = 0 ; i < n_ants ; i++ ) {
        free( ant[i].tour );
        free( ant[i].visited );
    }
    free( ant );
    free( best_so_far_ant->tour );
    free( best_so_far_ant->visited );
    free( prob_of_selection );
    free(demand_meet_node_map);
    return(0);
}
