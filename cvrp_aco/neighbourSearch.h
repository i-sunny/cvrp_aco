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
    long int rnd_seed;
    
    
    long int random_pos_in_route(Route *route);
    long int random_route();
    
    Move *exchange(long int *tour, long int tour_size);
    Move *insertion(long int *tour, long int tour_size);
    Move *inversion(long int *tour, long int tour_size);
};

#endif /* neighbourSearch_h */
