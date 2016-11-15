//
//  neighbourSearch.hpp
//  cvrp_aco
//
//  Created by 孙晓奇 on 2016/10/19.
//  Copyright © 2016年 xiaoqi.sxq. All rights reserved.
//

#ifndef neighbourSearch_h
#define neighbourSearch_h

#include <stdio.h>
#include <vector>
#include "problem.h"
#include "move.h"
#include "localSearch.h"

class NeighbourSearch {
public:
    NeighbourSearch(Problem *instance);
    ~NeighbourSearch();
    void reset_ant(AntStruct *ant);
    Move *search(AntStruct *ant);
    
private:
    Problem *instance;
    AntStruct *ant;
    vector<Route> routes;
    int rnd_seed;
    
    
    int random_pos_in_route(Route *route);
    int random_route();
    
    Move *exchange(int *tour, int tour_size);
    Move *exchange_1(int *tour, int tour_size);
    Move *insertion(int *tour, int tour_size);
    Move *insertion_1(int *tour, int tour_size);
    Move *inversion(int *tour, int tour_size);
};

#endif /* neighbourSearch_h */
