//
//  move.cpp
//  cvrp_aco
//
//  Created by 孙晓奇 on 2016/10/19.
//  Copyright © 2016年 xiaoqi.sxq. All rights reserved.
//

#include <assert.h>
#include "utilities.h"
#include "move.h"

/*
 * exchange node[pos_n1] and node[pos_n2]
 * pos_n1 < pos_n2
 */
void ExchangeMove::apply()
{
    long int *tour = ant->tour;
    
    long int tmp = tour[pos_n1];
    tour[pos_n1] = tour[pos_n2];
    tour[pos_n2] = tmp;
    
    ant->tour_length += gain;
}

/*
 * if pos_n1 < pos_n2, then insert node[pos_n1] after node[pos_n2]
 * if pos_n1 > pos_n2, then insert node[pos_n1] before node[pos_n1]
 */
void InsertionMove::apply()
{
    long int *tour = ant->tour;
    long int tmp = tour[pos_n1];
    long int i;
    
    if (pos_n1 < pos_n2) {
        for (i = pos_n1; i < pos_n2; i++) {
            tour[i] = tour[i+1];
        }
        tour[pos_n2] = tmp;
    } else {
        for (i = pos_n1; i > pos_n2; i--) {
            tour[i] = tour[i-1];
        }
        tour[pos_n2] = tmp;
    }
    
    ant->tour_length += gain;
}

/*
 * reverse nodes(in the same route) between pos_n1 and pos_n2,
 * pos_n1 and pos_n2 included.
 */
void InversionMove::apply()
{
    long int *tour = ant->tour;
    long int i, j;
    long int tmp;
    
    if (pos_n1 > pos_n2) {
        tmp = pos_n1;
        pos_n1 = pos_n2;
        pos_n2 = tmp;
    }
    
    i = pos_n1;
    j = pos_n2;
    DEBUG(assert(pos_n1 < pos_n2);)
    while (i < j) {
        tmp = tour[i];
        tour[i] = tour[j];
        tour[j] = tmp;
        i++;
        j--;
    }
    
    ant->tour_length += gain;
}