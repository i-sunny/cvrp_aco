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

#include "io.h"
#include "timer.h"
#include "utilities.h"
#include "antColony.h"
#include "vrpHelper.h"


static bool report_flag = TRUE;   /* 结果是否输出文件 */
static FILE *report, *best_so_far_report, *iter_report, *anneal_report;

void write_params(Problem *instance);
static void fprintf_parameters (FILE *stream, Problem *instance);

/*
 FUNCTION:       print ant soulution *tour
 INPUT:          pointer to a tour
 OUTPUT:         none
 */
void print_solution(Problem *instance, long int *tour, long int tour_size)
{
    long int   i;
    
    printf("\n--------\n");
    
    printf("Ant soulution:");
    for( i = 0 ; i < tour_size ; i++ ) {
        if (!i%25) {
            printf("\n");
        }
        printf("%ld ", tour[i]);
    }
    printf("\n");
    printf("Tour length = %f, pid = %d", compute_tour_length(instance, tour, tour_size), instance->pid);
    
    printf("\n--------\n\n");
}

/*
 FUNCTION:       print a single route from an ant's solution
 INPUT:          pointer to a route
 OUTPUT:         none
 */
void print_single_route(Problem *instance, long int *route, long int route_size)
{
    long int   i;
    
    printf("\n--------\n");
    
    printf("Route: ");
    for( i = 0 ; i < route_size ; i++ ) {
        if (i!= 0 && !i%25) {
            printf("\n");
        }
        printf("%ld ", route[i]);
    }
    printf("\n");
    printf("Route length = %f, pid = %d", compute_tour_length(instance, route, route_size), instance->pid);
    
    printf("\n--------\n");
}

void print_problem_decompositon(const vector<Problem *>& subs)
{
    Problem *sub;
    printf("\n--------\n");
    printf("Decompositons:\n");
    for (int i = 0; i < subs.size(); i++) {
        sub = subs[i];
        printf("sub %d, length %f:", i+1, sub->best_so_far_ant->tour_length);
        for (int j = 0; j < sub->real_nodes.size(); j++) {
            printf("%ld ", sub->real_nodes[j]);
        }
        printf("\n");
    }
   printf("--------\n");
}

/*
 FUNCTION:       prints the selection probabilities as encountered by an ant
 INPUT:          none
 OUTPUT:         none
 COMMENTS:       this computation assumes that no choice has been made yet.
 */
void print_probabilities(Problem *instance)
{
    long int i, j;
    double   *p;
    double   sum_prob;
    long int num_node = instance->num_node;
    
    printf("Selection Probabilities, iteration: %ld\n",instance->iteration);
    p = (double *)calloc(num_node, sizeof(double) );
    
    for (i=0; i < num_node; i++) {
        printf("From %ld:  ",i);
        sum_prob = 0.;
        for ( j = 0 ; j < num_node ; j++) {
            if ( i == j )
                p[j] = 0.;
            else
                p[j] = instance->total_info[i][j];
            sum_prob += p[j];
        }
        for ( j = 0 ; j < num_node ; j++) {
            p[j] = p[j] / sum_prob;
        }
        for ( j = 0 ; j < num_node-1 ; j++) {
            printf(" %.5f ", p[j]);
        }
        printf(" %.5f\n", p[num_node-1]);
        if (!(j % 26)) {
            printf("\n");
        }
        printf("\n");
    }
    printf("\n");
    free ( p );
}

/*
 FUNCTION:       print distance matrix
 INPUT:          none
 OUTPUT:         none
 */
void print_distance(Problem *instance)
{
    long int i,j;
    
    printf("Distance Matrix:\n");
    for ( i = 0 ; i < instance->num_node ; i++) {
        printf("From %ld:  ",i);
        for ( j = 0 ; j < instance->num_node - 1 ; j++ ) {
            printf(" %f", instance->distance[i][j]);
        }
        printf(" %f\n", instance->distance[i][instance->num_node-1]);
        printf("\n");
    }
    printf("\n");
}

/*
 FUNCTION:       print pheromone trail values
 INPUT:          none
 OUTPUT:         none
 */
void print_pheromone(Problem *instance)
{
    long int i,j;
    
    printf("pheromone Trail matrix, iteration: %ld\n\n",instance->iteration);
    for ( i = 0 ; i < instance->num_node ; i++) {
        printf("From %ld:  ",i);
        for ( j = 0 ; j < instance->num_node ; j++ ) {
            printf(" %.10f ", instance->pheromone[i][j]);
            if (instance->pheromone[i][j] > 1.0)
                printf("XXXXX\n");
        }
        printf("\n");
    }
    printf("\n");
}

/*
 FUNCTION:       print values of pheromone times heuristic information
 INPUT:          none
 OUTPUT:         none
 */
void print_total_info(Problem *instance)
{
    long int i, j, num_node =  instance->num_node;
    double **total = instance->total_info;
    
    printf("combined pheromone and heuristic info\n\n");
    for (i=0; i < num_node; i++) {
        for (j = 0; j < num_node - 1 ; j++) {
            printf(" %.15f &", total[i][j]);
            if ( total[i][j] > 1.0 )
                printf("XXXXX\n");
        }
        printf(" %.15f\n", total[i][num_node-1]);
        if ( total[i][num_node-1] > 1.0 )
            printf("XXXXX\n");
    }
    printf("\n");
}

/*
 * 问题开始时
 */
void init_report(Problem *instance, long int ntry)
{
    printf("\n############### start try %ld ###############\n", ntry);
    
    char temp_buffer[LINE_BUF_LEN];
    
    if (report_flag) {
        sprintf(temp_buffer,"./report/best.%s",instance->name);
        report = fopen(temp_buffer, "w");
        
        sprintf(temp_buffer,"./report/best_so_far.%s",instance->name);
        best_so_far_report = fopen(temp_buffer, "a");
        
        sprintf(temp_buffer,"./report/iter.%s",instance->name);
        iter_report = fopen(temp_buffer, "w");
        
        sprintf(temp_buffer,"./report/anneal.%s",instance->name);
        anneal_report = fopen(temp_buffer, "w");
    } else {
        report = NULL;
        anneal_report = NULL;
        best_so_far_report = NULL;
        iter_report = NULL;
    }
    
    if (best_so_far_report) {
        fprintf(best_so_far_report,"\n############### start try %ld ###############\n", ntry);
    }
    write_params(instance);
}

/*
 * 问题结束时
 */
void exit_report(Problem *instance, long int ntry) {
    
    if(!check_solution(instance, instance->best_so_far_ant->tour, instance->best_so_far_ant->tour_size)) {
        exit(EXIT_FAILURE);
    }
    
    printf("\n\nBest Length: %f\t Iterations: %ld\t At time %.2f\t Tot.time %.2f\n",
            instance->best_so_far_ant->tour_length, instance->iteration, instance->best_so_far_time, elapsed_time(VIRTUAL));
    printf("############### end try %ld ###############\n\n", ntry);
    
    if (report) {
        fprintf(report, "Best Length: %f\t Iterations: %ld\t At time %.2f\t Tot.time %.2f\n",
                instance->best_so_far_ant->tour_length, instance->iteration, instance->best_so_far_time, elapsed_time(VIRTUAL));
        fflush(report);
    }

    if (best_so_far_report){
        print_solution_to_file(instance, best_so_far_report, instance->best_so_far_ant->tour, instance->best_so_far_ant->tour_size);
        fprintf(best_so_far_report,"############### end try %ld ###############\n\n",ntry);
        fflush(best_so_far_report);
    }
    if (iter_report) {
        fflush(iter_report);
    }
    if (anneal_report) {
        fflush(anneal_report);
    }
}

/*
 FUNCTION:       print the solution *t to best_so_far.vrplibfile
 INPUT:          pointer to a tour
 OUTPUT:         none
 */
void print_solution_to_file(Problem *instance, FILE *file, long int *tour, long int tour_size)
{
    if (file){
        fprintf(file,"Begin Solution\n");
        for(int i = 0 ; i < tour_size ; i++ ) {
            fprintf(file, "%ld ", tour[i]);
        }
        fprintf(file,"\n");
        fprintf(file,"Solution Length %f\n", compute_tour_length(instance, tour, tour_size));
        fprintf(file,"End Solution\n");
    }
}

/*
 FUNCTION: output some info about best-so-far solution quality, and its time
 INPUT:    none
 OUTPUT:   none
 COMMENTS: none
 */
void write_best_so_far_report(Problem *instance)
{
    printf("best so far length %f, iteration: %ld, time %.2f\n",
           instance->best_so_far_ant->tour_length, instance->iteration, instance->best_so_far_time);
    if (best_so_far_report) {
        fprintf(best_so_far_report, "%f\t %ld\t %.3f\n",
                instance->best_so_far_ant->tour_length, instance->iteration, instance->best_so_far_time);
    }
}

/*
 FUNCTION: output some info about best-so-far solution quality, and its time
 INPUT:    none
 OUTPUT:   none
 COMMENTS: none
 */
void write_iter_report(Problem *instance)
{
    DEBUG(printf("iteration: %ld, iter best length %f, time %.2f\n",
           instance->iteration, instance->iteration_best_ant->tour_length, elapsed_time( VIRTUAL));)
    if (iter_report) {
        fprintf(iter_report, "%f\t %ld\t %.3f\n",
                instance->iteration_best_ant->tour_length, instance->iteration, elapsed_time( VIRTUAL));
    }
}

/*
 * 退火算法report
 */
void write_anneal_report(Problem *instance, AntStruct *ant, Move *move)
{
    if (InversionMove *p = dynamic_cast<InversionMove *>(move)) {
        DEBUG(printf("[Inversion Move] moved length %f, gain:%f, pos_n1:%ld, pos_n2:%ld, time %.2f\n",
               ant->tour_length, p->gain, p->pos_n1, p->pos_n2, elapsed_time( VIRTUAL));)
        if (anneal_report) {
            fprintf(anneal_report, "[Inversion Move]: best length %f, gain:%f, pos_n1:%ld, pos_n2:%ld, time %.2f\n",
                    ant->tour_length, p->gain, p->pos_n1, p->pos_n2, elapsed_time( VIRTUAL));
        }
    } else if (InsertionMove *p = dynamic_cast<InsertionMove *>(move)) {
        DEBUG(printf("[Insertion Move] moved length %f, gain:%f, pos_n1:%ld, pos_n2:%ld, load_r1:%ld, load_r2:%ld, time %.2f\n",
               ant->tour_length, p->gain, p->pos_n1, p->pos_n2, p->load_r1, p->load_r2, elapsed_time( VIRTUAL));)
        if (anneal_report) {
            fprintf(anneal_report, "[Insertion Move]: best length %f, gain:%f, pos_n1:%ld, pos_n2:%ld, time %.2f\n",
                    ant->tour_length, p->gain, p->pos_n1, p->pos_n2, elapsed_time( VIRTUAL));
        }
    } else if (ExchangeMove *p = dynamic_cast<ExchangeMove *>(move)) {
        DEBUG(printf("[Exchange Move] moved length %f, gain:%f, pos_n1:%ld, pos_n2:%ld, load_r1:%ld, load_r2:%ld, time %.2f\n",
               ant->tour_length, p->gain, p->pos_n1, p->pos_n2, p->load_r1, p->load_r2, elapsed_time( VIRTUAL));)
        if (anneal_report) {
            fprintf(anneal_report, "[Exchange Move] moved length %f, gain:%f, pos_n1:%ld, pos_n2:%ld, load_r1:%ld, load_r2:%ld, time %.2f\n",
                    ant->tour_length, p->gain, p->pos_n1, p->pos_n2, p->load_r1, p->load_r2, elapsed_time( VIRTUAL));
        }
    }
//    print_solution(instance, ant->tour, ant->tour_size);
    print_solution_to_file(instance, anneal_report, ant->tour, ant->tour_size);
    
}

void write_params(Problem *instance)
/*
 FUNCTION:       writes chosen parameter settings in standard output and in
 report files
 INPUT:          none
 OUTPUT:         none
 */
{
    fprintf(stdout, "\nParameter-settings: \n\n");
    fprintf_parameters (stdout, instance);
    fprintf(stdout, "\n");
    
    if (report) {
        fprintf(report,"\nParameter-settings: \n\n");
        fprintf_parameters (report, instance);
        fprintf(report,"\n");
    }
}

static void fprintf_parameters (FILE *stream, Problem *instance)
{
    fprintf(stream,"max_time\t\t %.2f\n", instance->max_runtime);
    fprintf(stream,"seed\t\t %ld\n", instance->rnd_seed);
    fprintf(stream,"optimum\t\t\t %f\n", instance->optimum);
    fprintf(stream,"n_ants\t\t\t %ld\n", instance->n_ants);
    fprintf(stream,"nn_ants\t\t\t %ld\n", instance->nn_ants);
    fprintf(stream,"alpha\t\t\t %.2f\n", alpha);
    fprintf(stream,"beta\t\t\t %.2f\n", beta);
    fprintf(stream,"rho\t\t\t %.2f\n", rho);
    fprintf(stream,"ras_ranks\t\t %ld\n", ras_ranks);
    fprintf(stream,"ls_flag\t\t\t %d\n", instance->ls_flag);
    fprintf(stream,"nn_ls\t\t\t %ld\n", instance->nn_ls);
    fprintf(stream,"dlb_flag\t\t %d\n", instance->dlb_flag);
}

/*
 FUNCTION: parse and read instance file
 INPUT:    instance name
 OUTPUT:   none
 COMMENTS: Instance files have to be in vrpLIB format, otherwise procedure fails
 */
void read_instance_file(Problem *instance, const char *vrp_file_name)
{
    FILE         *vrp_file;
    char         buf[LINE_BUF_LEN];
    long int     i, j;
    Point *nodeptr;
    
    vrp_file = fopen(vrp_file_name, "r");
    if ( vrp_file == NULL ) {
        fprintf(stderr,"No instance file specified, abort\n");
        exit(1);
    }
    printf("\nreading vrp-file %s ... \n\n", vrp_file_name);
    
    instance->max_distance = LONG_MAX;
    instance->service_time = 0;
    instance->optimum = 0;
    
    fscanf(vrp_file,"%s", buf);
    while ( strcmp("NODE_COORD_SECTION", buf) != 0 ) {
        if ( strcmp("NAME", buf) == 0 ) {
            fscanf(vrp_file, "%s", buf);
            fscanf(vrp_file, "%s", buf);
            strcpy(instance->name, buf);
            buf[0]=0;
        }
        else if ( strcmp("NAME:", buf) == 0 ) {
            fscanf(vrp_file, "%s", buf);
            strcpy(instance->name, buf);
            buf[0]=0;
        }
        else if ( strcmp("COMMENT", buf) == 0 ){
            fscanf(vrp_file, "%lf", &instance->optimum);
            buf[0]=0;
        }
        else if ( strcmp("COMMENT:", buf) == 0 ){
            fscanf(vrp_file, "%lf", &instance->optimum);
            buf[0]=0;
        }
        else if ( strcmp("TYPE", buf) == 0 ) {
            fscanf(vrp_file, "%s", buf);
            fscanf(vrp_file, "%s", buf);
            if( strcmp("CVRP", buf) != 0 ) {
                fprintf(stderr,"\n Not a vrp instance in vrpLIB format !!\n");
                exit(1);
            }
            buf[0]=0;
        }
        else if ( strcmp("TYPE:", buf) == 0 ) {
            fscanf(vrp_file, "%s", buf);
            if( strcmp("CVRP", buf) != 0 ) {
                fprintf(stderr,"\n Not a vrp instance in vrpLIB format !!\n");
                exit(1);
            }
            buf[0]=0;
        }
        else if( strcmp("DIMENSION", buf) == 0 ){
            fscanf(vrp_file, "%s", buf);
            fscanf(vrp_file, "%ld", &instance->num_node);
            buf[0]=0;
        }
        else if ( strcmp("DIMENSION:", buf) == 0 ) {
            fscanf(vrp_file, "%ld", &instance->num_node);
            buf[0]=0;
        }
        else if( strcmp("DISPLAY_DATA_TYPE", buf) == 0 ){
            fgets(buf, LINE_BUF_LEN, vrp_file);
            buf[0]=0;
        }
        else if ( strcmp("DISPLAY_DATA_TYPE:", buf) == 0 ) {
            fgets(buf, LINE_BUF_LEN, vrp_file);
            buf[0]=0;
        }
        else if( strcmp("EDGE_WEIGHT_TYPE", buf) == 0 ){
            buf[0]=0;
            fscanf(vrp_file, "%s", buf);
            buf[0]=0;
            fscanf(vrp_file, "%s", buf);
            if ( strcmp("EUC_2D", buf) == 0 ) {
                instance->dis_type = DIST_EUC_2D;
            }
            else if ( strcmp("CEIL_2D", buf) == 0 ) {
                instance->dis_type = DIST_CEIL_2D;
            }
            else if ( strcmp("GEO", buf) == 0 ) {
                instance->dis_type = DIST_GEO;
            }
            else if ( strcmp("ATT", buf) == 0 ) {
                instance->dis_type = DIST_ATT;
            }
            else
                fprintf(stderr,"EDGE_WEIGHT_TYPE %s not implemented\n",buf);
            strcpy(instance->edge_weight_type, buf);
            buf[0]=0;
        }
        else if( strcmp("EDGE_WEIGHT_TYPE:", buf) == 0 ){
            /* set pointer to appropriate distance function; has to be one of
             EUC_2D, CEIL_2D, GEO, or ATT. Everything else fails */
            buf[0]=0;
            fscanf(vrp_file, "%s", buf);
            printf("%s\n", buf);
            printf("%s\n", buf);
            if ( strcmp("EUC_2D", buf) == 0 ) {
                instance->dis_type = DIST_EUC_2D;
            }
            else if ( strcmp("CEIL_2D", buf) == 0 ) {
                instance->dis_type = DIST_CEIL_2D;
            }
            else if ( strcmp("GEO", buf) == 0 ) {
                instance->dis_type = DIST_GEO;
            }
            else if ( strcmp("ATT", buf) == 0 ) {
                instance->dis_type = DIST_ATT;
            }
            else {
                fprintf(stderr,"EDGE_WEIGHT_TYPE %s not implemented\n",buf);
                exit(1);
            }
            strcpy(instance->edge_weight_type, buf);
            buf[0]=0;
        }
        else if( strcmp("CAPACITY", buf) == 0 ){
            fscanf(vrp_file, "%s", buf);
            fscanf(vrp_file, "%ld", &instance->vehicle_capacity);
            buf[0]=0;
        }
        else if ( strcmp("CAPACITY:", buf) == 0 ) {
            fscanf(vrp_file, "%ld", &instance->vehicle_capacity);
            buf[0]=0;
        }
        else if ( strcmp("DISTANCE", buf) == 0 ) {
            fscanf(vrp_file, "%lf", &instance->max_distance);
            buf[0]=0;
        }
        else if ( strcmp("DISTANCE:", buf) == 0 ) {
            fscanf(vrp_file, "%lf", &instance->max_distance);
            buf[0]=0;
        }
        else if ( strcmp("SERVICE_TIME", buf) == 0 ) {
            fscanf(vrp_file, "%lf", &instance->service_time);
            buf[0]=0;
        }
        else if ( strcmp("SERVICE_TIME:", buf) == 0 ) {
            fscanf(vrp_file, "%lf", &instance->service_time);
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
    
    if((nodeptr = (Point *)malloc(sizeof(Point) * instance->num_node)) == NULL)
        exit(EXIT_FAILURE);
    else {
        for ( i = 0 ; i < instance->num_node ; i++ ) {
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
    for ( i = 0 ; i < instance->num_node ; i++ ) {
        fscanf(vrp_file,"%ld %ld", &j, &nodeptr[i].demand);
    }
    
    instance->nodeptr = nodeptr;
    
    TRACE ( printf("\n... done\n"); )
}


