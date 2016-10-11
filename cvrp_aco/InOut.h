/*********************************
Ant Colony Optimization algorithms (AS, ACS, EAS, RAS, MMAS) for CVRP

Created by 孙晓奇 on 2016/10/8.
Copyright © 2016年 xiaoqi.sxq. All rights reserved.

Program's name: acovrp
Purpose: mainly input / output / statistic routines

email: sunxq1991@gmail.com

*********************************/


#define LINE_BUF_LEN     255

extern long int iteration;           /* counter of number iterations */
extern long int max_iteration;         /* maximum number of iterations */
extern long int seed;

extern double   max_time;                   /* maximal allowed run time of a try  */
extern double   best_so_far_time;         /* 最优解出现的时间 */

extern long int optimal;                 /* optimal solution or bound to find */
extern long int best_solution_iter;   /* iteration in which best solution is found */
extern long int *demand_meet_node_map;   /** 所有可配送的点(单次route中，目前车辆可以仍可配送的点) */

extern long int master_problem_iteration_num;   /* 每次外循环，主问题蚁群的迭代的次数 */
extern long int sub_problem_iteration_num;      /* 每次外循环，子问题蚁群的迭代的次数 */
extern long int num_sub_problems;               /* 拆分子问题个数 */


void init_program( long int argc, char *argv[] );
void exit_program(void);

struct point * read_instance_file(const char *vrp_file_name);

const char* parse_commandline (long int argc, char *argv []);

void set_default_parameters (void);

void print_solution(long int ant_id, long int *tour, long int tour_size);

void print_single_route(long int route_id, long int *route, long int route_size);

void printSolutionToFile(FILE *file, long int *tour, long int tour_size);

void write_best_so_far_report(void);

void write_iter_report(void);
