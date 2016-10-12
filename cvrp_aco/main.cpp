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
#include "problem.h"
#include "timer.h"
#include "io.h"

/*
 FUNCTION:       checks whether termination condition is met
 INPUT:          none
 OUTPUT:         0 if condition is not met, number neq 0 otherwise
 (SIDE)EFFECTS:  none
 */
bool termination_condition(Problem *instance)
{
    return ((instance->iteration >= instance->max_iteration) ||
            (elapsed_time( VIRTUAL ) >= g_max_runtime) ||
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
int main(int argc, char *argv[]) {
    
    Problem *instance;
    AntColony *solver;
    
    start_timers();
    
    const char *filename = parse_commandline(argc, argv);
    
    instance = init_master_problem(filename);
    
    g_best_so_far_time = elapsed_time( VIRTUAL );
    printf("Initialization took %.10f seconds\n", g_best_so_far_time);
    
    solver = new AntColony(instance);
    solver->init_aco();
    
    while (!termination_condition(instance)) {
        
        solver->run_aco_iteration();
        instance->iteration++;
    }
    solver->exit_aco();
    
    exit_master_problem(instance);
    
    return(0);
}
