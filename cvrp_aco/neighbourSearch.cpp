///*********************************
//
// Ant Colony Optimization algorithms RAS for CVRP
//
// Created by 孙晓奇 on 2016/10/19.
// Copyright © 2016年 xiaoqi.sxq. All rights reserved.
//
// Program's name: acovrp
// Purpose: performs three types of neighborhood search: sequence inversion, insertion, and exchange
//
// email: sunxq1991@gmail.com
//
// *********************************/

#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <math.h>
#include <assert.h>
#include<time.h>

#include "neighbourSearch.h"
#include "utilities.h"

using namespace std;

NeighbourSearch::NeighbourSearch(Problem *instance)
{
    srandom((unsigned int)time(NULL));
    this->instance = instance;
}

NeighbourSearch::~NeighbourSearch()
{
    routes.clear();
}

/*
 *
 */
void NeighbourSearch::reset_ant(AntStruct *ant)
{
    long int *tour = ant->tour;
    Point *nodes = instance->nodeptr;
    long int beg;
    long int load;
    
    this->ant = ant;
    
    routes.clear();
    beg = 0;
    load = 0;
    for (int i = 1; i < ant->tour_size; i++) {
        if (tour[i] == 0) {
            Route route(beg, i, load);
            routes.push_back(route);
            
            beg = i;
            load = 0;
        } else {
            load += nodes[tour[i]].demand;
        }
    }
}


Move *NeighbourSearch::search(AntStruct *ant)
{
    // 首先需要reset solution
    reset_ant(ant);
    
    long int *tour = ant->tour;
    long int tour_size = ant->tour_size;
    
    long int rnd = random() % 3;
    switch (rnd) {
        case 0:
            return exchange(tour, tour_size);
            break;
        case 1:
            return insertion(tour, tour_size);
            break;
        case 2:
            return inversion(tour, tour_size);
            break;
        default:
            return exchange(tour, tour_size);
            break;
    }
}

/*
 * function: randomly exchanging two nodes from two routes
 */
Move *NeighbourSearch::exchange(long int *tour, long int tour_size)
{
    long int n1 = 0, n2 = 0;   /* random node from route 1 and toure 2*/
    long int p_n1, p_n2, s_n1, s_n2;
    long int pos_n1 = 0, pos_n2 = 0;
    long int r1, r2;           /* idx of route 1 and route 2 */
    long int gain;
    long int **distance = instance->distance;
    Point *nodes = instance->nodeptr;
    long int load_r1, load_r2;
    bool valid = true;
    
    r1 = random_route();
    r2 = r1;
    while (r1 == r2) {
        r2 = random_route();
    }
    
    while (n1 == 0) {
        pos_n1 = random_pos_in_route(&routes[r1]);
        n1 = tour[pos_n1];
    }
    while (n2 == 0) {
        pos_n2 = random_pos_in_route(&routes[r2]);
        n2 = tour[pos_n2];
    }
    
    DEBUG(assert(pos_n1 > 0 && pos_n2 > 0);)
    DEBUG(assert(r1 != r2);)
    DEBUG(assert(n1 != 0 && n2 != 0);)
    
    p_n1 = tour[pos_n1-1];
    p_n2 = tour[pos_n2-1];
    s_n1 = tour[pos_n1+1];
    s_n2 = tour[pos_n2+1];
    
    gain = -(distance[p_n1][n1] + distance[n1][s_n1] + distance[p_n2][n2] + distance[n2][s_n2])
    +(distance[p_n1][n2] + distance[n2][s_n1] + distance[p_n2][n1] + distance[n1][s_n2]);
    
    load_r1 = routes[r1].load - nodes[n1].demand + nodes[n2].demand;
    load_r2 = routes[r2].load - nodes[n2].demand + nodes[n1].demand;
    
    if (load_r1 > instance->vehicle_capacity || load_r2 > instance->vehicle_capacity) {
        valid = false;
    }
    return new ExchangeMove(ant, valid, gain, pos_n1, pos_n2, load_r1, load_r2);
}

/*
 * function: randomly insert a node(n1) from route1 to route2 after pos_n2
 */
Move *NeighbourSearch::insertion(long int *tour, long int tour_size)
{
    long int n1 = 0, n2 = 0;   /* random node from route 1 and toure 2*/
    long int p_n1, s_n1, s_n2;
    long int pos_n1 = 0, pos_n2 = 0;
    long int r1, r2;           /* idx of route 1 and route 2 */
    long int gain;
    long int **distance = instance->distance;
    Point *nodes = instance->nodeptr;
    long int load_r1, load_r2;
    bool valid = true;
    
    r1 = random_route();
    r2 = r1;
    while (r1 == r2) {
        r2 = random_route();
    }
    
    while (n1 == 0) {
        pos_n1 = random_pos_in_route(&routes[r1]);
        n1 = tour[pos_n1];
    }

    pos_n2 = random_pos_in_route(&routes[r2]);
    n2 = tour[pos_n2];

    DEBUG(assert(r1 != r2);)
    DEBUG(assert(n1 != 0);)
    
    p_n1 = tour[pos_n1-1];
    s_n1 = tour[pos_n1+1];
    s_n2 = tour[pos_n2+1];
    
    gain = -(distance[p_n1][n1] + distance[n1][s_n1] + distance[n2][s_n2]) + (distance[n2][n1] + distance[n1][s_n2] + distance[p_n1][s_n1]);
    
    load_r1 = routes[r1].load - nodes[n1].demand;
    load_r2 = routes[r2].load + nodes[n1].demand;
    
    if (load_r2 > instance->vehicle_capacity) {
        valid = false;
    }
    return new InsertionMove(ant, valid, gain, pos_n1, pos_n2, load_r1, load_r2);
}

/*
 * chooses two nodes in a route randomly, and then inverts the substring between these two nodes(n1, n2)
 */
Move *NeighbourSearch::inversion(long int *tour, long int tour_size)
{
    Route *route = NULL;
    long int r_sz = 0;     /* route size */
    long int n1, n2;   /* random node from route 1 and toure 2*/
    long int pos_n1 = 0, pos_n2 = 0, tmp;
    long int s_n1, p_n2;
    long int gain;
    long int **distance = instance->distance;
    
    while (r_sz <= 5) {
        route = &routes[random_route()];
        r_sz = route->end -route->beg + 1;
    }
    
    pos_n1 = random_pos_in_route(route);
    
    pos_n2 = pos_n1;
    while (abs(pos_n2 - pos_n1) <= 1) {
        pos_n2 = random_pos_in_route(route);
    }
    
    if (pos_n1 > pos_n2) {
        tmp = pos_n1;
        pos_n1 = pos_n2;
        pos_n2 = tmp;
    }
    
    n1 = tour[pos_n1];
    n2 = tour[pos_n2];
    
    s_n1 = tour[pos_n1+1];
    p_n2 = tour[pos_n2-1];
    
    DEBUG(assert(n1 != n2);)
    DEBUG(assert(pos_n2 - pos_n1 > 1 && pos_n1 < route->end);)
    
    gain = -(distance[n1][s_n1] + distance[p_n2][n2]) + (distance[n1][p_n2] + distance[s_n1][n2]);
    
    return new InversionMove(ant, true, gain, pos_n1, pos_n2);
}


/*
 * get random idx from route
 * idx = [beg, end-1]
 */
long int NeighbourSearch::random_pos_in_route(Route *route)
{
    return (route->beg + random()%(route->end - route->beg));
}


long int NeighbourSearch::random_route(void)
{
    return random()%(routes.size());
}

