/*********************************
 Ant Colony Optimization algorithms (AS, ACS, EAS, RAS, MMAS) for CVRP
 
 Created by 孙晓奇 on 2016/10/8.
 Copyright © 2016年 xiaoqi.sxq. All rights reserved.
 
 Program's name: acovrp
 Purpose: some additional useful procedures
 
 email: sunxq1991@gmail.com
 
 *********************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "utilities.h"
#include "timer.h"


double mean( int *values, int max ) 
/*    
      FUNCTION:       compute the average value of an integer array of length max 
      INPUT:          pointer to array, length of array
      OUTPUT:         average 
      (SIDE)EFFECTS:  none
*/
{
  int j;
  double   m;

  m = 0.;
  for ( j = 0 ; j < max ; j++ ) {
    m += (double)values[j];
  }
  m = m / (double)max;
  return m;
}



double meanr( double *values, int max ) 
/*    
      FUNCTION:       compute the average value of a floating number array of length max 
      INPUT:          pointer to array, length of array
      OUTPUT:         average 
      (SIDE)EFFECTS:  none
*/
{
  int j;
  double   m;

  m = 0.;
  for ( j = 0 ; j < max ; j++ ) {
    m += values[j];
  }
  m = m / (double)max;
  return m;
}



double std_deviation( int *values, int max, double mean ) 
/*    
      FUNCTION:       compute the standard deviation of an integer array  
      INPUT:          pointer to array, length of array, mean 
      OUTPUT:         standard deviation
      (SIDE)EFFECTS:  none
*/
{
  int j;
  double   dev = 0.;

  if (max <= 1)
    return 0.;
  for ( j = 0 ; j < max; j++ ) {
    dev += ((double)values[j] - mean) * ((double)values[j] - mean);
  }
  return sqrt(dev/(double)(max - 1));
}



double std_deviationr( double *values, int max, double mean ) 
/*    
      FUNCTION:       compute the standard deviation of a floating number array  
      INPUT:          pointer to array, length of array, mean 
      OUTPUT:         standard deviation
      (SIDE)EFFECTS:  none
*/
{
  int j;
  double   dev;

  if (max <= 1)
    return 0.;
  dev = 0.;
  for ( j = 0 ; j < max ; j++ ) {
    dev += ((double)values[j] - mean) * ((double)values[j] - mean);
  }
  return sqrt(dev/(double)(max - 1));
}



int best_of_vector( int *values, int l ) 
/*    
      FUNCTION:       return the minimum value in an integer value  
      INPUT:          pointer to array, length of array
      OUTPUT:         smallest number in the array
      (SIDE)EFFECTS:  none
*/
{
  int min, k;

  k = 0;
  min = values[k];
  for( k = 1 ; k < l ; k++ ) {
    if( values[k] < min ) {
      min = values[k];
    }
  }
  return min;
}



int worst_of_vector( int *values, int l ) 
/*    
      FUNCTION:       return the maximum value in an integer value  
      INPUT:          pointer to array, length of array
      OUTPUT:         largest number in the array
      (SIDE)EFFECTS:  none
*/
{
  int max, k;

  k = 0;
  max = values[k];
  for( k = 1 ; k < l ; k++ ) {
    if( values[k] > max ){
      max = values[k];
    }
  }
  return max;
}



double quantil(int v[], double q, int l)
/*    
      FUNCTION:       return the q-quantil of an ordered integer array  
      INPUT:          one array, desired quantil q, length of array
      OUTPUT:         q-quantil of array
      (SIDE)EFFECTS:  none
*/
{
  int i,j;
  double tmp;

  tmp = q * (double)l;
  if ((double)((int)tmp) == tmp) {  
    i = (int)tmp;
    j = (int)(tmp + 1.);
    return ((double)v[i-1] + (double)v[j-1]) / 2.;
  } else {
    i = (int)(tmp +1.);
    return v[i-1];
  }
}


void swap(int *i, int *j)
{
    int tmp = *i;
    *i = *j;
    *j = tmp;
}

void swap(int v[], int i, int j)
/*    
      FUNCTION:       auxiliary routine for sorting an integer array  
      INPUT:          array, two indices
      OUTPUT:         none
      (SIDE)EFFECTS:  elements at position i and j of array are swapped
*/
{
  int tmp;

  tmp = v[i];
  v[i] = v[j];
  v[j] = tmp;
}




void sort(int v[], int left, int right)
/*    
      FUNCTION:       recursive routine (quicksort) for sorting an array  
      INPUT:          one array, two indices
      OUTPUT:         none
      (SIDE)EFFECTS:  elements at position i and j of the two arrays are swapped
*/
{
  int k, last;

  if (left >= right) 
    return;
  swap(v, left, (left + right)/2);
  last = left;
  for (k=left+1; k <= right; k++)
    if (v[k] < v[left])
      swap(v, ++last, k);
  swap(v, left, last);
  sort(v, left, last);
  sort(v, last+1, right);
}



void swap2(double v[], int v2[], int i, int j)
/*    
      FUNCTION:       auxiliary routine for sorting an integer array  
      INPUT:          two arraya, two indices
      OUTPUT:         none
      (SIDE)EFFECTS:  elements at position i and j of the two arrays are swapped
*/
{
    double tmp1;
    int tmp2;

    tmp1 = v[i];
    v[i] = v[j];
    v[j] = tmp1;

    tmp2 = v2[i];
    v2[i] = v2[j];
    v2[j] = tmp2;
}



void sort2(double v[], int v2[], int left, int right)
/*    
      FUNCTION:       recursive routine (quicksort) for sorting one array; second 
                      arrays does the same sequence of swaps  
      INPUT:          two arrays, two indices
      OUTPUT:         none
      (SIDE)EFFECTS:  elements at position i and j of the two arrays are swapped
*/
{
  int k, last;

  if (left >= right) 
    return;
  swap2(v, v2, left, (left + right)/2);
  last = left;
  for (k=left+1; k <= right; k++)
    if (v[k] < v[left])
      swap2(v, v2, ++last, k);
  swap2(v, v2, left, last);
  sort2(v, v2, left, last);
  sort2(v, v2, last+1, right);
}



double ran01( int *idum )
/*    
      FUNCTION:       generate a random number that is uniformly distributed in [0,1]
      INPUT:          pointer to variable with the current seed
      OUTPUT:         random number uniformly distributed in [0,1]
      (SIDE)EFFECTS:  random number seed is modified (important, this has to be done!)
      ORIGIN:         numerical recipes in C
*/
{
  int k;
  double ans;

  k =(*idum)/IQ;
  *idum = IA * (*idum - k * IQ) - IR * k;
  if (*idum < 0 ) *idum += IM;
  ans = AM * (*idum);
  return ans;
}



int random_number( int *idum )
/*    
      FUNCTION:       generate an integer random number
      INPUT:          pointer to variable containing random number seed
      OUTPUT:         integer random number uniformly distributed in {0,2147483647}
      (SIDE)EFFECTS:  random number seed is modified (important, has to be done!)
      ORIGIN:         numerical recipes in C
*/
{
  int k;

  k =(*idum)/IQ;
  *idum = IA * (*idum - k * IQ) - IR * k;
  if (*idum < 0 ) *idum += IM;
  return *idum;
}



int ** generate_int_matrix( int n, int m)
/*    
      FUNCTION:       malloc a matrix and return pointer to it
      INPUT:          size of matrix as n x m 
      OUTPUT:         pointer to matrix
      (SIDE)EFFECTS:  
*/
{
  int i;
  int **matrix;

  if((matrix = (int **)malloc(sizeof(int) * n * m +
                                   sizeof(int *) * n	 )) == NULL){
    printf("Out of memory, exit.");
    exit(1);
  }
  for ( i = 0 ; i < n ; i++ ) {
    matrix[i] = (int *)(matrix + n) + i*m;
  }

  return matrix;
}



double ** generate_double_matrix( int n, int m)
/*    
      FUNCTION:       malloc a matrix and return pointer to it
      INPUT:          size of matrix as n x m 
      OUTPUT:         pointer to matrix
      (SIDE)EFFECTS:  
*/
{

  int i;
  double **matrix;

  if((matrix = (double **)malloc(sizeof(double) * n * m +
                                 sizeof(double *) * n	 )) == NULL){
    printf("Out of memory, exit.");
    exit(1);
  }
  for ( i = 0 ; i < n ; i++ ) {
    matrix[i] = (double *)(matrix + n) + i*m;
  }
  return matrix;
}
