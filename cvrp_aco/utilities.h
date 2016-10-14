/*********************************
 Ant Colony Optimization algorithms (AS, ACS, EAS, RAS, MMAS) for CVRP
 
 Created by 孙晓奇 on 2016/10/8.
 Copyright © 2016年 xiaoqi.sxq. All rights reserved.
 
 Program's name: acovrp
 Purpose: some additional useful procedures
 
 email: sunxq1991@gmail.com
 
 *********************************/

#define INFTY                 LONG_MAX
#define MAXIMUM_NO_TRIES      100

#define TRUE  1
#define FALSE 0

/* general macros */

#define MAX(x,y)        ((x)>=(y)?(x):(y))
#define MIN(x,y)        ((x)<=(y)?(x):(y))

#define DEBUG( x )  x

#define TRACE( x )

/* constants for a random number generator, for details see numerical recipes in C */

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876

double mean ( long int *values, long int max);

double meanr ( double *values, long int max );

double std_deviation ( long int *values, long int i, double mean );

double std_deviationr ( double *values, long int i, double mean );

long int best_of_vector ( long int *values, long int i );

long int worst_of_vector ( long int *values, long int i );

void swap ( long int v[], long int i, long int j );

void sort ( long int v[], long int left, long int right );

double quantil ( long int vector[], double q, long int numbers );

void swap2(long int v[], long int v2[], long int i, long int j);

void sort2(long int v[], long int v2[], long int left, long int right);

double ran01 ( long *idum );

long int random_number ( long *idum );

long int ** generate_int_matrix( long int n, long int m);

double ** generate_double_matrix( long int n, long int m);
