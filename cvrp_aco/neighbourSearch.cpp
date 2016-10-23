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
#include "io.h"

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
    double dist;
    
    this->ant = ant;
    
    routes.clear();
    beg = 0;
    load = 0;
    dist = 0;
    for (int i = 1; i < ant->tour_size; i++) {
        load += nodes[tour[i]].demand;
        dist += instance->distance[tour[i-1]][tour[i]];
        
        if (tour[i] == 0) {
            Route route(beg, i, load, dist);
            routes.push_back(route);
            
            beg = i;
            load = 0;
            dist = 0;
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
    long int n1, n2;          /* random node from route 1 and toure 2*/
    long int p_n1, p_n2, s_n1, s_n2;
    long int pos_n1 = 0, pos_n2 = 0;
    long int r1, r2;           /* idx of route 1 and route 2 */
    double gain;
    double **distance = instance->distance;
    Point *nodes = instance->nodeptr;
    long int load_r1, load_r2;
    double dist_r1, dist_r2;
    bool valid = true;
    
    r1 = random_route();
    r2 = random_route();
    if (r1 == r2) {
        return exchange_1(tour, tour_size);
    }
    
    n1 = 0;
    while (n1 == 0) {
        pos_n1 = random_pos_in_route(&routes[r1]);
        n1 = tour[pos_n1];
    }
    
    n2 = 0;
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
    
    dist_r1 = routes[r1].dist + (routes[r1].end - routes[r1].beg - 1) * instance->service_time
    - (distance[p_n1][n1] + distance[n1][s_n1]) + (distance[p_n1][n2] + distance[n2][s_n1]);
    
    dist_r2 = routes[r2].dist + (routes[r2].end - routes[r2].beg - 1) * instance->service_time
    - (distance[p_n2][n2] + distance[n2][s_n2]) + (distance[p_n2][n1] + distance[n1][s_n2]);
    
    load_r1 = routes[r1].load - nodes[n1].demand + nodes[n2].demand;
    load_r2 = routes[r2].load - nodes[n2].demand + nodes[n1].demand;
    
    if ((load_r1 > instance->vehicle_capacity || load_r2 > instance->vehicle_capacity)
        ||(dist_r1 > instance->max_distance || dist_r2 > instance->max_distance))
    {
        valid = false;
    }
    return new ExchangeMove(ant, valid, gain, pos_n1, pos_n2, load_r1, load_r2);
}

/*
 * function: randomly exchanging two non-zero nodes from one route
 */
Move *NeighbourSearch::exchange_1(long int *tour, long int tour_size)
{
    Route *route = NULL;
    long int n1, n2;
    long int p_n1, p_n2, s_n1, s_n2;
    long int pos_n1, pos_n2;
    double gain;
    double **distance = instance->distance;
    long int sz;
    double dist;
    bool valid = true;
    
    sz = 0;
    while (sz <= 3) {
        route = &routes[random_route()];
        sz = route->end - route->beg + 1;
    }
    
    // non-zero n1
    pos_n1 = route->beg;
    while (pos_n1 == route->beg) {
        pos_n1 = random_pos_in_route(route);
    }
    
    // non-zero n2
    pos_n2 = route->beg;
    while (pos_n2 == route->beg || pos_n1 == pos_n2) {
        pos_n2 = random_pos_in_route(route);
    }
    
    if (pos_n1 > pos_n2) {
        swap(&pos_n1, &pos_n2);
    }
    n1 = tour[pos_n1];
    n2 = tour[pos_n2];
    
    DEBUG(assert(pos_n1 > 0 && pos_n2 > 0 && pos_n1 < pos_n2);)
    DEBUG(assert(n1 != 0 && n2 != 0 && n1 != n2);)
    
    p_n1 = tour[pos_n1-1];
    p_n2 = tour[pos_n2-1];
    s_n1 = tour[pos_n1+1];
    s_n2 = tour[pos_n2+1];
    
    if (pos_n2 - pos_n1 == 1) {
        gain = -(distance[p_n1][n1] + distance[n2][s_n2]) +(distance[p_n1][n2] + distance[n1][s_n2]);
    } else {
        gain = -(distance[p_n1][n1] + distance[n1][s_n1] + distance[p_n2][n2] + distance[n2][s_n2])
        +(distance[p_n1][n2] + distance[n2][s_n1] + distance[p_n2][n1] + distance[n1][s_n2]);
    }
    
    dist = route->dist + (route->end - route->beg - 1) * instance->service_time + gain;
    if(dist > instance->max_distance) {
        valid = false;
    }
    
    return new ExchangeMove(ant, valid, gain, pos_n1, pos_n2, route->load, route->load);
}

/*
 * function: randomly insert a node(n1) from route1 to route2 after pos_n2
 * if pos_n1 < pos_n2, then insert node[pos_n1] after node[pos_n2]
 * if pos_n1 > pos_n2, then insert node[pos_n1] before node[pos_n1]
 */
Move *NeighbourSearch::insertion(long int *tour, long int tour_size)
{
    Route *route1, *route2;
    long int n1, n2;   /* random node from route 1 and toure 2*/
    long int p_n1, p_n2, s_n1, s_n2;
    long int pos_n1, pos_n2;
    long int r1 = 0, r2;           /* idx of route 1 and route 2 */
    double gain;
    double **distance = instance->distance;
    Point *nodes = instance->nodeptr;
    long int load_r1, load_r2;
    double dist_r2;
    bool valid = true;
    long int sz;
    
    sz = 0;
    while (sz <= 3) {
        r1 = random_route();
        sz = routes[r1].end - routes[r1].beg + 1;
    }
    
    r2 = random_route();
    if (r1 == r2) {
        return insertion_1(tour, tour_size);
    }
    
    // r1 != r2
    route1 = &routes[r1];
    route2 = &routes[r2];
    
    // non-zero n1 from route 1
    pos_n1 = route1->beg;
    while (pos_n1 == route1->beg) {
        pos_n1 = random_pos_in_route(route1);
    }

    // non-zero n2 from route 2
    pos_n2 = route2->beg;
    while (pos_n2 == route2->beg) {
        pos_n2 = random_pos_in_route(route2);
    }
    
    n1 = tour[pos_n1];
    n2 = tour[pos_n2];

    DEBUG(assert(r1 != r2);)
    DEBUG(assert(n1 != 0 && n2 != 0);)
    
    p_n1 = tour[pos_n1-1];
    p_n2 = tour[pos_n2-1];
    s_n1 = tour[pos_n1+1];
    s_n2 = tour[pos_n2+1];
    
    if (pos_n1 > pos_n2) {
        gain = -(distance[p_n1][n1] + distance[n1][s_n1] + distance[p_n2][n2])
        + (distance[n1][n2] + distance[p_n2][n1] + distance[p_n1][s_n1]);
        // r2多了一个元素
        dist_r2 = route2->dist + (route2->end - route2->beg) * instance->service_time
        - (distance[p_n2][n2]) + (distance[n1][n2] + distance[p_n2][n1]);
        
    } else {
        gain = -(distance[p_n1][n1] + distance[n1][s_n1] + distance[n2][s_n2])
        + (distance[n2][n1] + distance[n1][s_n2] + distance[p_n1][s_n1]);
        
        dist_r2 = route2->dist + (route2->end - route2->beg) * instance->service_time
        - (distance[n2][s_n2]) + (distance[n2][n1] + distance[n1][s_n2]);
        
    }
    
    load_r1 = route1->load - nodes[n1].demand;
    load_r2 = route2->load + nodes[n1].demand;
    
    if (load_r2 > instance->vehicle_capacity || dist_r2 > instance->max_distance) {
        valid = false;
    }
    return new InsertionMove(ant, valid, gain, pos_n1, pos_n2, load_r1, load_r2);
}

/*
 * insert node n1 after n2 in the same route
 * if pos_n1 < pos_n2, then insert node[pos_n1] after node[pos_n2]
 * if pos_n1 > pos_n2, then insert node[pos_n1] before node[pos_n1]
 */
Move *NeighbourSearch::insertion_1(long int *tour, long int tour_size)
{
    Route *route = NULL;
    long int n1, n2;
    long int p_n1, p_n2, s_n1, s_n2;
    long int pos_n1, pos_n2;
    double gain;
    double **distance = instance->distance;
    long int load_r1, load_r2;
    double dist;
    long int sz;
    bool valid = true;
    
    sz = 0;
    while (sz <= 3) {
        route = &routes[random_route()];
        sz = route->end - route->beg + 1;
    }
    
    // non-zero n1
    pos_n1 = route->beg;
    while (pos_n1 == route->beg) {
        pos_n1 = random_pos_in_route(route);
    }
    
    // non-zero n2
    pos_n2 = route->beg;
    while (pos_n2 == route->beg || pos_n1 == pos_n2) {
        pos_n2 = random_pos_in_route(route);
    }
    
    n1 = tour[pos_n1];
    n2 = tour[pos_n2];
    
    DEBUG(assert(n1 != 0 && n2 != 0);)
    
    p_n1 = tour[pos_n1-1];
    p_n2 = tour[pos_n2-1];
    s_n1 = tour[pos_n1+1];
    s_n2 = tour[pos_n2+1];
    
    if (pos_n1 > pos_n2) {
        gain = -(distance[p_n1][n1] + distance[n1][s_n1] + distance[p_n2][n2])
        + (distance[n1][n2] + distance[p_n2][n1] + distance[p_n1][s_n1]);
    } else {
        gain = -(distance[p_n1][n1] + distance[n1][s_n1] + distance[n2][s_n2])
        + (distance[n2][n1] + distance[n1][s_n2] + distance[p_n1][s_n1]);
    }
    
    dist = route->dist + (route->end - route->beg - 1) * instance->service_time + gain;
    if(dist > instance->max_distance) {
        valid = false;
    }
    
    return new InsertionMove(ant, valid, gain, pos_n1, pos_n2, load_r1, load_r2);
}

/*
 * chooses two nodes in a route randomly, 
 * and then inverts the substring between these two non-zero nodes[n1, n2],
 * pos_n1 and pos_n2 included.
 */
Move *NeighbourSearch::inversion(long int *tour, long int tour_size)
{
    Route *route = NULL;
    long int sz = 0;     /* route size */
    long int n1, n2;   /* random node from route 1 and toure 2*/
    long int pos_n1 = 0, pos_n2 = 0;
    long int p_n1, s_n2;
    double gain;
    double **distance = instance->distance;
    double dist;
    bool valid = true;
    
    sz = 0;
    while (sz <= 3) {
        route = &routes[random_route()];
        sz = route->end -route->beg + 1;
    }
    
    // non-zero n1
    pos_n1 = route->beg;
    while (pos_n1 == route->beg) {
        pos_n1 = random_pos_in_route(route);
    }
    
    // non-zero n2, n1 != n2
    pos_n2 = route->beg;
    while (pos_n2 == route->beg || pos_n1 == pos_n2) {
        pos_n2 = random_pos_in_route(route);
    }
    
    if (pos_n1 > pos_n2) {
        swap(&pos_n1, &pos_n2);
    }
    
    // useless
    if (pos_n2 == route->end - 1 && pos_n1 == route->beg + 1) {
        return NULL;
    }
    
    n1 = tour[pos_n1];
    n2 = tour[pos_n2];
    
    p_n1 = tour[pos_n1-1];
    s_n2 = tour[pos_n2+1];
    
    DEBUG(assert(n1 != n2);)
    DEBUG(assert(pos_n1 > 0 && pos_n2 > 0 && pos_n1 < pos_n2);)
    
    gain = -(distance[p_n1][n1] + distance[n2][s_n2]) + (distance[p_n1][n2] + distance[n1][s_n2]);
    
    dist = route->dist + (route->end - route->beg - 1) * instance->service_time + gain;
    if(dist > instance->max_distance) {
        valid = false;
    }
    
    return new InversionMove(ant, valid, gain, pos_n1, pos_n2);
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


