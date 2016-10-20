//
//  main.cpp
//  cvrp_aco
//
//  Created by 孙晓奇 on 2016/10/11.
//  Copyright © 2016年 xiaoqi.sxq. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "utilities.h"
#include "antColony.h"
#include "parallelAco.h"
#include "simulatedAnnealing.h"
#include "problem.h"
#include "timer.h"
#include "io.h"

static bool parallel_flag  = false;  /* 是否使用并行算法 */
/*
 FUNCTION:       checks whether termination condition is met
 INPUT:          none
 OUTPUT:         0 if condition is not met, number neq 0 otherwise
 (SIDE)EFFECTS:  none
 */
bool termination_condition(Problem *instance)
{
    return ((instance->iteration >= instance->max_iteration) ||
            (elapsed_time( VIRTUAL ) >= instance->max_runtime) ||
            (instance->best_so_far_ant->tour_length <= instance->optimum));
}

/*
 * 解析命令行，获取文件名
 */
const char* parse_commandline (long int argc, char *argv [])
{
    const char *filename;
    
    if (argc <= 1) {
        fprintf (stderr,"Error: No vrp instance file.\n");
        exit(1);
    }
    
    filename = argv[1];
    return filename;
}

/* --- main program ------------------------------------------------------ */
/*
 FUNCTION:       main control for running the ACO algorithms
 INPUT:          none
 OUTPUT:         none
 (SIDE)EFFECTS:  none
 COMMENTS:       this function controls the run of "max_tries" independent trials
 
 */
int main(int argc, char *argv[])
{
    Problem *instance = new Problem(0);
    AntColony *solver;
    SimulatedAnnealing *annealer;
    
    start_timers();
    
    const char *filename = parse_commandline(argc, argv);
    read_instance_file(instance, filename);
    init_problem(instance);
    init_report(instance);
    
    printf("Initialization took %.10f seconds\n", elapsed_time(VIRTUAL));

    if (parallel_flag && instance->num_subs > 1) {
        solver = new ParallelAco(instance);
    } else {
        solver = new AntColony(instance);
    }
    
    solver->init_aco();
    
    while (!termination_condition(instance)) {
    
        solver->run_aco_iteration();
        instance->iteration++;

        if (instance->best_stagnate_cnt >= 2 * instance->num_node) {
            annealer = new SimulatedAnnealing(instance, solver, 2.0, 0.97, MAX(instance->num_node * 4, 250), 30);
            annealer->run();
            instance->best_stagnate_cnt = 0;
            
            delete annealer;
        }
    }

    solver->exit_aco();
    
    delete solver;
    exit_report(instance);
    exit_problem(instance);
    
    return(0);
}
