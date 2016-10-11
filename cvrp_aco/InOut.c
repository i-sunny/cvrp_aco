/*********************************
 
 Ant Colony Optimization algorithms (AS, ACS, EAS, RAS, MMAS) for CVRP
 
 Created by 孙晓奇 on 2016/10/8.
 Copyright © 2016年 xiaoqi.sxq. All rights reserved.
 
 Program's name: acovrp
 Purpose: mainly input / output / statistic routines
 
 email: sunxq1991@gmail.com
 
 *********************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <limits.h>
#include <time.h>

#include "InOut.h"
#include "vrp.h"
#include "timer.h"
#include "utilities.h"
#include "ants.h"
#include "ls.h"
#include "parallel_acovrp.h"

long int iteration;           /* counter of number iterations */
long int max_iteration;         /* maximum number of iterations */
long int seed;

double   max_time;          /* maximal allowed run time of a try  */
double   best_so_far_time;        /* 当前最优解出现的时间 */

long int optimal;                 /* optimal solution or bound to find */
long int best_solution_iter;      /* iteration in which best solution is found */
long int *demand_meet_node_map;   /** 所有可配送的点(单次route中，目前车辆可以仍可配送的点) */

/**** 用于parallel aco参数 ***/
long int master_problem_iteration_num;   /* 每次外循环，主问题蚁群的迭代的次数 */
long int sub_problem_iteration_num;      /* 每次外循环，子问题蚁群的迭代的次数 */
long int num_sub_problems;               /* 拆分子问题个数 */

/* ------------------------------------------------------------------------ */

FILE *report, *best_so_far_report, *iter_report;

long int report_flag;    /* 结果是否打印输出 */

/*
 FUNCTION: parse and read instance file
 INPUT:    instance name
 OUTPUT:   list of coordinates for all nodes
 COMMENTS: Instance files have to be in vrpLIB format, otherwise procedure fails
 */
struct point * read_instance_file(const char *vrp_file_name)
{
    FILE         *vrp_file;
    char         buf[LINE_BUF_LEN];
    long int     i, j;
    struct point *nodeptr;

    vrp_file = fopen(vrp_file_name, "r");
    if ( vrp_file == NULL ) {
        fprintf(stderr,"No instance file specified, abort\n");
        exit(1);
    }
    assert(vrp_file != NULL);
    printf("\nreading vrp-file %s ... \n\n", vrp_file_name);

    fscanf(vrp_file,"%s", buf);
    while ( strcmp("NODE_COORD_SECTION", buf) != 0 ) {
        if ( strcmp("NAME", buf) == 0 ) {
            fscanf(vrp_file, "%s", buf);
            TRACE ( printf("%s ", buf); )
            fscanf(vrp_file, "%s", buf);
            strcpy(instance.name, buf);
            TRACE ( printf("%s \n", instance.name); )
            buf[0]=0;
        }
        else if ( strcmp("NAME:", buf) == 0 ) {
            fscanf(vrp_file, "%s", buf);
            strcpy(instance.name, buf);
            TRACE ( printf("%s \n", instance.name); )
            buf[0]=0;
        }
        else if ( strcmp("COMMENT", buf) == 0 ){
            fgets(buf, LINE_BUF_LEN, vrp_file);
            TRACE ( printf("%s", buf); )
            buf[0]=0;
        }
        else if ( strcmp("COMMENT:", buf) == 0 ){
            fgets(buf, LINE_BUF_LEN, vrp_file);
            TRACE ( printf("%s", buf); )
            buf[0]=0;
        }
        else if ( strcmp("TYPE", buf) == 0 ) {
            fscanf(vrp_file, "%s", buf);
            TRACE ( printf("%s ", buf); )
            fscanf(vrp_file, "%s", buf);
            TRACE ( printf("%s\n", buf); )
            if( strcmp("CVRP", buf) != 0 ) {
                fprintf(stderr,"\n Not a vrp instance in vrpLIB format !!\n");
            exit(1);
            }
            buf[0]=0;
        }
        else if ( strcmp("TYPE:", buf) == 0 ) {
            fscanf(vrp_file, "%s", buf);
            TRACE ( printf("%s\n", buf); )
            if( strcmp("CVRP", buf) != 0 ) {
                fprintf(stderr,"\n Not a vrp instance in vrpLIB format !!\n");
            exit(1);
            }
            buf[0]=0;
        }
        else if( strcmp("DIMENSION", buf) == 0 ){
            fscanf(vrp_file, "%s", buf);
            TRACE ( printf("%s ", buf); );
            fscanf(vrp_file, "%ld", &num_node);
            instance.num_node = num_node;
            TRACE ( printf("%ld\n", num_node); );
            buf[0]=0;
        }
        else if ( strcmp("DIMENSION:", buf) == 0 ) {
            fscanf(vrp_file, "%ld", &num_node);
            instance.num_node = num_node;
            TRACE ( printf("%ld\n", num_node); );
            buf[0]=0;
        }
        else if( strcmp("DISPLAY_DATA_TYPE", buf) == 0 ){
            fgets(buf, LINE_BUF_LEN, vrp_file);
            TRACE ( printf("%s", buf); );
            buf[0]=0;
        }
        else if ( strcmp("DISPLAY_DATA_TYPE:", buf) == 0 ) {
            fgets(buf, LINE_BUF_LEN, vrp_file);
            TRACE ( printf("%s", buf); );
            buf[0]=0;
        }
        else if( strcmp("EDGE_WEIGHT_TYPE", buf) == 0 ){
            buf[0]=0;
            fscanf(vrp_file, "%s", buf);
            TRACE ( printf("%s ", buf); );
            buf[0]=0;
            fscanf(vrp_file, "%s", buf);
            TRACE ( printf("%s\n", buf); );
            if ( strcmp("EUC_2D", buf) == 0 ) {
                distance = round_distance;
            }
            else if ( strcmp("CEIL_2D", buf) == 0 ) {
                distance = ceil_distance;
            }
            else if ( strcmp("GEO", buf) == 0 ) {
                distance = geo_distance;
            }
            else if ( strcmp("ATT", buf) == 0 ) {
                distance = att_distance;
            }
            else
                fprintf(stderr,"EDGE_WEIGHT_TYPE %s not implemented\n",buf);
            strcpy(instance.edge_weight_type, buf);
            buf[0]=0;
        }
        else if( strcmp("EDGE_WEIGHT_TYPE:", buf) == 0 ){
            /* set pointer to appropriate distance function; has to be one of 
               EUC_2D, CEIL_2D, GEO, or ATT. Everything else fails */
            buf[0]=0;
            fscanf(vrp_file, "%s", buf);
            TRACE ( printf("%s\n", buf); )
            printf("%s\n", buf);
            printf("%s\n", buf);
            if ( strcmp("EUC_2D", buf) == 0 ) {
                distance = round_distance;
            }
            else if ( strcmp("CEIL_2D", buf) == 0 ) {
                distance = ceil_distance;
            }
            else if ( strcmp("GEO", buf) == 0 ) {
                distance = geo_distance;
            }
            else if ( strcmp("ATT", buf) == 0 ) {
                distance = att_distance;
            }
            else {
                fprintf(stderr,"EDGE_WEIGHT_TYPE %s not implemented\n",buf);
                exit(1);
            }
            strcpy(instance.edge_weight_type, buf);
            buf[0]=0;
        }
        else if( strcmp("CAPACITY", buf) == 0 ){
            fscanf(vrp_file, "%s", buf);
            TRACE ( printf("%s ", buf); );
            fscanf(vrp_file, "%ld", &vehicle_capacity);
            instance.vehicle_capacity = vehicle_capacity;
            TRACE ( printf("%ld\n", vehicle_capacity); );
            buf[0]=0;
        }
        else if ( strcmp("CAPACITY:", buf) == 0 ) {
            fscanf(vrp_file, "%ld", &vehicle_capacity);
            instance.vehicle_capacity = vehicle_capacity;
            TRACE ( printf("%ld\n", vehicle_capacity); );
            buf[0]=0;
        }
        
        buf[0]=0;
        fscanf(vrp_file,"%s", buf);
    }


    // read coordinates
    if( strcmp("NODE_COORD_SECTION", buf) == 0 ){
        TRACE ( printf("found section contaning the node coordinates\n"); )
    }
    else{
        fprintf(stderr,"\n\nSome error ocurred finding start of coordinates from vrp file !!\n");
        exit(1);
    }

    if( (nodeptr = malloc(sizeof(struct point) * num_node)) == NULL )
        exit(EXIT_FAILURE);
    else {
        for ( i = 0 ; i < num_node ; i++ ) {
            fscanf(vrp_file,"%ld %lf %lf", &j, &nodeptr[i].x, &nodeptr[i].y );
        }
    }
    
    // read demand
    while(strcmp("DEMAND_SECTION", buf) != 0) {
        fscanf(vrp_file,"%s", buf);
    }
    if(strcmp("DEMAND_SECTION", buf) == 0) {
        TRACE ( printf("found section contaning the node demand\n"); )
    }
    for ( i = 0 ; i < num_node ; i++ ) {
        fscanf(vrp_file,"%ld %ld", &j, &nodeptr[i].demand);
    }
    
    TRACE ( printf("number of cities is %ld\n",n); )
    TRACE ( printf("\n... done\n"); )
	return (nodeptr);
}

void set_default_parameters ()
{
    n_ants         = num_node;    /* number of ants*/
    nn_ants        = MIN(MAX(n_ants>>2, 25), num_node - 2);    /* number of nearest neighbours in tour construction(neighbor不应该包括depot和自身) */
    nn_ls          = MIN(nn_ants, 20);                         /* use fixed radius search in the 20 nearest neighbours */
    
    alpha          = 1.0;
    beta           = 2.0;
    rho            = 0.1;
    q_0            = 0.0;
    ras_ranks      = 6;          /* number of ranked ants, top-{ras_ranks} ants */
    
    // local search
    ras_flag       = TRUE;
    ls_flag        = 1;     /* apply local search */
    dlb_flag       = TRUE;  /* apply don't look bits in local search */

    seed           = (long int) time(NULL);
    max_time       = 20.0;
    max_iteration  = 5000;
    optimal        = 1;
    iteration      = 0;
    
    report_flag    = 1;
    
    // parallel aco
    parallel_flag                   = TRUE;
    master_problem_iteration_num    = 1;
    sub_problem_iteration_num       = 75;
    num_sub_problems                = num_node/50;
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

/*
 FUNCTION:       initialize the program,
 INPUT:          program arguments, needed for parsing commandline
 OUTPUT:         none
 COMMENTS:
 */
void init_program( long int argc, char *argv[] )

{
    char temp_buffer[LINE_BUF_LEN];
    
    const char *filename = parse_commandline(argc, argv);
    instance.nodeptr = read_instance_file(filename);
    
    set_default_parameters();

    instance.distance = compute_distances();
    instance.nn_list = compute_nn_lists();
    
    pheromone = generate_double_matrix( num_node, num_node );
    total_info = generate_double_matrix( num_node, num_node );
    
    demand_meet_node_map = calloc(num_node, sizeof(long int));
    
    allocate_ants();
    
    if (report_flag) {
        sprintf(temp_buffer,"best.%s",instance.name);
        report = fopen(temp_buffer, "w");
        
        sprintf(temp_buffer,"best_so_far.%s",instance.name);
        best_so_far_report = fopen(temp_buffer, "w");
        
        sprintf(temp_buffer,"iter.%s",instance.name);
        iter_report = fopen(temp_buffer, "w");
    } else {
        best_so_far_report = NULL;
        iter_report = NULL;
    }
}

/*
 * 结束程序
 */
void exit_program(void)
{
    if (report) {
        fprintf(report, "Best Length: %ld\t Iterations: %ld\t At time %.2f\t Tot.time %.2f\n",
                best_so_far_ant->tour_length, iteration, best_so_far_time, elapsed_time( VIRTUAL ));
        fflush(report);
    }
    
    if (best_so_far_report){
        vrp_check_solution(best_so_far_ant->tour, best_so_far_ant->tour_size);
        printSolutionToFile(best_so_far_report, best_so_far_ant->tour, best_so_far_ant->tour_size);
        fflush(best_so_far_report);
    }
    if (iter_report) {
        fflush(iter_report);
    }
}


/**************************************************************************
 **************************************************************************
                    print
 ***************************************************************************
 **************************************************************************/

/*
 FUNCTION:       print ant soulution *tour
 INPUT:          pointer to a tour
 OUTPUT:         none
 */
void print_solution(long int ant_id, long int *tour, long int tour_size)
{
    long int   i;
    
    printf("\n--------\n");
    
    printf("Ant %ld soulution:", ant_id);
    for( i = 0 ; i < tour_size ; i++ ) {
        if (!i%25) {
            printf("\n");
        }
        printf("%ld ", tour[i]);
    }
    printf("\n");
    printf("Tour length = %ld",compute_tour_length(tour, tour_size));
    
    printf("\n--------\n\n");
}

/*
 FUNCTION:       print a single route from an ant's solution
 INPUT:          pointer to a route
 OUTPUT:         none
 */
void print_single_route(long int route_id, long int *route, long int route_size)
{
    long int   i;
    
    printf("\n--------\n");
    
    printf("Route %ld : ", route_id);
    for( i = 0 ; i < route_size ; i++ ) {
        if (i!= 0 && !i%25) {
            printf("\n");
        }
        printf("%ld ", route[i]);
    }
    printf("\n");
    printf("Route length = %ld",compute_tour_length(route, route_size));
    
    printf("\n--------\n");
}

/*
 FUNCTION:       print the solution *t to best_so_far.vrplibfile
 INPUT:          pointer to a tour
 OUTPUT:         none
 */
void printSolutionToFile(FILE *file, long int *tour, long int tour_size)
{
    if (file){
        fprintf(file,"Begin Solution\n");
        for(int i = 0 ; i < tour_size ; i++ ) {
            fprintf(file, "%ld ", tour[i]);
        }
        fprintf(file,"\n");
        fprintf(file,"Solution Length %ld\n", compute_tour_length(tour, tour_size));
        fprintf(file,"End Solution\n");
    }
}

/*
 FUNCTION: output some info about best-so-far solution quality, and its time
 INPUT:    none
 OUTPUT:   none
 COMMENTS: none
 */
void write_best_so_far_report( void )
{
    printf("best so far length %ld, iteration: %ld, time %.2f\n", best_so_far_ant->tour_length, iteration, elapsed_time( VIRTUAL));
    if (best_so_far_report) {
        fprintf(best_so_far_report, "best %ld\t iteration %ld\t time %.3f\n",
                best_so_far_ant->tour_length, iteration, best_so_far_time);
    }
}

/*
 FUNCTION: output some info about best-so-far solution quality, and its time
 INPUT:    none
 OUTPUT:   none
 COMMENTS: none
 */
void write_iter_report(void)
{
    printf("iteration: %ld, iter best length %ld, time %.2f\n",iteration, iteration_best_ant->tour_length, elapsed_time( VIRTUAL));
    if (iter_report) {
        fprintf(iter_report, "iteration %ld\t best %ld\t time %.3f\n",
                iteration, iteration_best_ant->tour_length, elapsed_time( VIRTUAL));
    }
}











