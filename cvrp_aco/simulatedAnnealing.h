//
//  simulatedAnnealing.h
//  cvrp_aco
//
//  Created by 孙晓奇 on 2016/10/19.
//  Copyright © 2016年 xiaoqi.sxq. All rights reserved.
//

#ifndef simulatedAnnealing_h
#define simulatedAnnealing_h

#include <stdio.h>
#include <vector>
#include "problem.h"
#include "neighbourSearch.h"
#include "antColony.h"

#define TABU_LENGTH  3

struct Tabu {
    Tabu(Move move_, long int life_):move(move_), life(life_){}
    Move move;
    long int life;   /* tabu life */
};

class SimulatedAnnealing {
public:
    SimulatedAnnealing(Problem *instance, AntColony *ant_colony, double t0,
                       double alpha, long int epoch_length, long int terminal_ratio);
    ~SimulatedAnnealing();
    void run(void);
    bool step(void);
    
private:
    Problem *instance;
    AntStruct *best_ant;
    AntStruct *iter_ant;
    AntColony *ant_colony;
    NeighbourSearch *neighbour_search;
    LocalSearch *local_search;
    vector<Tabu> tabu_list;
    double alpha;
    double t0;
    double t;
    long int iteration;
    long int epoch_length;
    long int epoch_counter;
    long int terminal_ratio;
    
    // 用于统计
    long int test_cnt;
    long int improvement_cnt;
    long int accept_cnt;
    
    bool acceptable(Move *move);
    void accept(Move *move);
    void reject(Move *move);
    
    // tabu list
    bool is_tabu(Move *move);
    void update_tabu_list(Move *move);
    
};

#endif /* simulatedAnnealing_h */
