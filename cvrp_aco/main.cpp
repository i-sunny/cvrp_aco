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
#include <math.h>

#include "utilities.h"
#include "antColony.h"
#include "parallelAco.h"
#include "simulatedAnnealing.h"
#include "problem.h"
#include "timer.h"
#include "io.h"

static bool parallel_flag  = true;  /* 是否使用并行算法 */
static int tries = 15;

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
            (fabs(instance->best_so_far_ant->tour_length - instance->optimum) < 10 * EPSILON));
}

/*
 * 解析命令行，获取文件名
 */
char* parse_commandline (int argc, char *argv [])
{
    char *filename;
    
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
    for (int i = 1; i <= 2; i++) {
    for (int ntry = 0 ; ntry < tries; ntry++)
    {
        Problem *instance = new Problem(0);
        AntColony *solver;
        
        start_timers();
        
        char *filename = parse_commandline(argc, argv);
//        sprintf(filename, "../dataset/CMT/CMT%d.vrp", i);
        sprintf(filename, "../dataset/Golden/Golden_%d.vrp", i);
        
        read_instance_file(instance, filename);
        init_problem(instance);
        init_report(instance, ntry);
        
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
        }

        solver->exit_aco();
        
        delete solver;
        
        exit_report(instance, ntry);
        exit_problem(instance);
    }
    }
    
    return(0);
}
