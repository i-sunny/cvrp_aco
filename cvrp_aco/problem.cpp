/*********************************
 Ant Colony Optimization algorithms (AS, ACS, EAS, RAS, MMAS) for CVRP
 
 Created by 孙晓奇 on 2016/10/8.
 Copyright © 2016年 xiaoqi.sxq. All rights reserved.
 
 Program's name: acovrp
 Purpose: vrp problem model, problem init/exit, check
 
 email: sunxq1991@gmail.com
 
 *********************************/

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

#include "problem.h"
#include "io.h"
#include "utilities.h"
#include "vrpHelper.h"
#include "timer.h"


/**** 用于parallel aco参数 ***/
int g_master_problem_iteration_num;   /* 每次外循环，主问题蚁群的迭代的次数 */
int g_sub_problem_iteration_num;      /* 每次外循环，子问题蚁群的迭代次数 */
bool sa_flag = true;                        /* 是否使用sa */

double rho;           /* parameter for evaporation */
double alpha;         /* importance of trail */
double beta;          /* importance of heuristic evaluate */
int ras_ranks;   /* additional parameter for rank-based version of ant system */

/* ------------------------------------------------------------------------ */

void set_default_parameters (Problem *instance);
void allocate_ants (Problem *instance);

/*
 * 初始化问题
 */
void init_problem(Problem *instance)
{
    set_default_parameters(instance);
    
    // 为 problem 实例的成员分配内存
    // 只有主问题需要计算distance矩阵
    if (instance->pid == 0) {
        instance->distance = compute_distances(instance);
    }
    instance->nn_list = compute_nn_lists(instance);
    instance->pheromone = generate_double_matrix(instance->num_node, instance->num_node);
    instance->total_info = generate_double_matrix(instance->num_node, instance->num_node );
    allocate_ants(instance);
}

void exit_problem(Problem *instance)
{
    // 释放内存
    free( instance->distance );
    free(instance->nodeptr);
    free( instance->nn_list );
    free( instance->pheromone );
    free( instance->total_info );
    for (int i = 0 ; i < instance->n_ants ; i++ ) {
        free( instance->ants[i].tour );
        free( instance->ants[i].visited );
        free( instance->ants[i].candidate );
    }
    free( instance->ants );
    free( instance->best_so_far_ant->tour );
    free( instance->best_so_far_ant->visited );
    free( instance->prob_of_selection );
    free(instance);
}


/*
 * 初始子问题
 * Note: 子问题的distance可以从主问题直接赋值，减少重复计算
 */
void init_sub_problem(Problem *master, Problem *sub)
{
    double **sub_dis;
    int ri, rj;
    Point *nodeptr, *m_node;
    
    // 初始化 sub-problem distance矩阵
    if((sub_dis = (double **)malloc(sizeof(double) * sub->num_node * sub->num_node +
                                      sizeof(double *) * sub->num_node)) == NULL) {
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < sub->num_node; i++ ) {
        sub_dis[i] = (double *)(sub_dis + sub->num_node) + i * sub->num_node;
        ri = sub->real_nodes[i];
        for (int j = 0; j < sub->num_node; j++ ) {
            rj = sub->real_nodes[j];
            sub_dis[i][j] = master->distance[ri][rj];
        }
    }
    sub->distance = sub_dis;
    
    // 初始化nodeptr
    if((nodeptr = (Point *)malloc(sizeof(Point) * sub->num_node)) == NULL) {
        exit(EXIT_FAILURE);
    }
    for (int i = 0 ; i < sub->num_node ; i++ ) {
        ri = sub->real_nodes[i];
        m_node = &master->nodeptr[ri];
        nodeptr[i].x = m_node->x;
        nodeptr[i].y = m_node->y;
        nodeptr[i].demand = m_node->demand;
    }
    sub->nodeptr = nodeptr;
    
    // sub 需要额外的结构存储当前最优解所对应的信息素
    sub->best_pheromone = generate_double_matrix(sub->num_node, sub->num_node);
//    print_distance(sub);
    
    init_problem(sub);
    
    sub->max_iteration = g_sub_problem_iteration_num;
    sub->dis_type = master->dis_type;
    sub->vehicle_capacity = master->vehicle_capacity;
    
    sub->max_distance = master->max_distance;
    sub->service_time = master->service_time;
}

void exit_sub_problem(Problem *sub)
{
    free(sub->best_pheromone);
    exit_problem(sub);
}

/*
 FUNCTION:       allocate the memory for the ant colony, the best-so-far and
 the iteration best ant
 INPUT:          none
 OUTPUT:         none
 (SIDE)EFFECTS:  allocation of memory for the ant colony and two ants that
 store intermediate tours
 
 */
void allocate_ants (Problem *instance)
{
    int i;
    AntStruct *ants, *best_so_far_ant;
    double   *prob_of_selection;
    
    if((ants = (AntStruct *)malloc(sizeof( AntStruct ) * instance->n_ants +
                                   sizeof(AntStruct *) * instance->n_ants)) == NULL){
        printf("Out of memory, exit.");
        exit(1);
    }
    for (i = 0 ; i < instance->n_ants ; i++) {
        ants[i].tour        = (int *)calloc(2*instance->num_node-1, sizeof(int));   // tour最长为2 * num_node - 1
        ants[i].visited     = (bool *)calloc(instance->num_node, sizeof(bool));
        ants[i].candidate = (bool *)calloc(instance->num_node, sizeof(bool));
    }
    
    if((best_so_far_ant = (AntStruct *)malloc(sizeof(AntStruct))) == NULL){
        printf("Out of memory, exit.");
        exit(1);
    }
    best_so_far_ant->tour        = (int *)calloc(2*instance->num_node-1, sizeof(int));
    best_so_far_ant->visited     = (bool *)calloc(instance->num_node, sizeof(bool));
    
    if ((prob_of_selection = (double *)malloc(sizeof(double) * (instance->nn_ants + 1))) == NULL) {
        printf("Out of memory, exit.");
        exit(1);
    }
    /* Ensures that we do not run over the last element in the random wheel.  */
    prob_of_selection[instance->nn_ants] = HUGE_VAL;
    
    instance->ants = ants;
    instance->best_so_far_ant = best_so_far_ant;
    instance->prob_of_selection = prob_of_selection;
}

/*
 * 参数设置
 */
void set_default_parameters (Problem *instance)
{
    /* number of ants */
    instance->n_ants         = instance->num_node;
    /* number of nearest neighbours in tour construction(neighbor不应该包括depot和自身) */
    instance->nn_ants        = instance->num_node - 2;//MIN(MAX(instance->n_ants>>2, 25), instance->num_node - 2);
    /* use fixed radius search in the 20 nearest neighbours */
    instance->nn_ls          = instance->nn_ants;//MIN(instance->nn_ants, 25);
    
    /* maximum number of iterations */
    instance->max_iteration  = 100;
    /* optimal tour length if known, otherwise a bound */
//    instance->optimum        = 1;
    /* counter of number iterations */
    instance->iteration      = 0;
    
    /* apply local search */
    instance->ls_flag        = TRUE;
    /* apply don't look bits in local search */
    instance->dlb_flag       = TRUE;
    
    alpha          = 1.0;
    beta           = 2.0;
    rho            = 0.1;
    ras_ranks      = 6;          /* number of ranked ants, top-{ras_ranks} ants */
    
    instance->rnd_seed       = (int) time(NULL);
    instance->max_runtime    = 600.0;
    
    // parallel aco
    g_master_problem_iteration_num    = 1;      /* 每次外循环，主问题蚁群的迭代次数 */
    g_sub_problem_iteration_num       = 75;     /* 每次外循环，子问题蚁群的迭代次数 */
    instance->num_subs                = instance->num_node/50;
    
}

/*
 * 检查 ant vrp solution 的有效性
 * i.e. tour = [0,1,4,2,0,5,3,0] toute1 = [0,1,4,2,0] route2 = [0,5,3,0]
 */
bool check_solution(Problem *instance, int *tour, int tour_size)
{
    int i;
    int * used;
    int num_node = instance->num_node;
    int route_beg;    /* 单条回路起点 */
    
    used = (int *)calloc (num_node, sizeof(int));
    
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
                fprintf(stderr,"\n%s:error: solution vector has two times the value %d (last position: %d)\n", __FUNCTION__, tour[i], i);
                goto error;
            } else {
                used[tour[i]] = TRUE;
            }
        }
        
        if (tour[i] == 0) {
            // 形成单条回路
            if(!check_route(instance, tour, route_beg, i)) {
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
        fprintf(stderr, " %d", tour[i]);
    fprintf(stderr,"\n");
    free(used);
    return FALSE;
}

/*
 * 检查单条回路有效性
 * 注意depot点出现两次（分别出现在首尾）
 * i.e. toute = [0,1,4,2,0]
 */
bool check_route(Problem *instance, int *tour, int rbeg, int rend)
{
    int load = 0;
    double distance = 0;
    
    if (tour[rbeg] != 0 || tour[rend] != 0) {
        fprintf(stderr,"\n%s:error: 车辆路径没有形成一条回路\n", __FUNCTION__);
        return FALSE;
    }
    if (rend - rbeg < 2) {
        fprintf(stderr,"\n%s:error: 单条回路长度不对. rbeg=%d, rend=%d\n", __FUNCTION__, rbeg, rend);
        return FALSE;
    }
    for (int i = rbeg + 1; i < rend; i++) {
        load += instance->nodeptr[tour[i]].demand;
    }
    if (load > instance->vehicle_capacity) {
        fprintf(stderr,"\n%s:error: 单条回路超过车辆最大承载量 load = %d, capacity = %d rbeg = %d rend = %d\n",
                __FUNCTION__, load, instance->vehicle_capacity, rbeg, rend);
        return FALSE;
    }
    
    distance = compute_route_length(instance, tour + rbeg, rend - rbeg + 1);
    distance += instance->service_time * (rend - rbeg - 1);
    
    if (distance > instance->max_distance) {
        fprintf(stderr,"\n%s:error: 单条回路超过车辆最大路程 distance = %f, max_distance = %f rbeg = %d rend = %d\n",
                __FUNCTION__, distance, instance->max_distance, rbeg, rend);
        return FALSE;
    }
    
    return TRUE;
}

