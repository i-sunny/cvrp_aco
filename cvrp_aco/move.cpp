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

void ExchangeMove::apply()
{
    long int *tour = ant->tour;
    
    long int tmp = tour[pos_n1];
    tour[pos_n1] = tour[pos_n2];
    tour[pos_n2] = tmp;
    
    ant->tour_length += gain;
}

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
        for (i = pos_n1; i > pos_n2+1; i--) {
            tour[i] = tour[i-1];
        }
        tour[pos_n2+1] = tmp;
    }
    ant->tour_length += gain;
}

void InversionMove::apply()
{
    long int *tour = ant->tour;
    long int i, j;
    long int tmp;
    
    i = pos_n1 + 1;
    j = pos_n2 - 1;
    DEBUG(assert(pos_n1 + 1 < pos_n2);)
    while (i < j) {
        tmp = tour[i];
        tour[i] = tour[j];
        tour[j] = tmp;
        i++;
        j--;
    }
    
    ant->tour_length += gain;
}