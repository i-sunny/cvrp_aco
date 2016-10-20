///*********************************
//
// Ant Colony Optimization algorithms RAS for CVRP
//
// Created by 孙晓奇 on 2016/10/19.
// Copyright © 2016年 xiaoqi.sxq. All rights reserved.
//
// Program's name: acovrp
// Purpose: Simulated Annealing algorithm
//
// email: sunxq1991@gmail.com
//
// *********************************/

#include <math.h>

#include "simulatedAnnealing.h"
#include "utilities.h"
#include "timer.h"
#include "io.h"

SimulatedAnnealing::SimulatedAnnealing(Problem *instance, double t0, double alpha, long int epoch_length, long int terminal_ratio)
{
    this->instance = instance;
    this->t = t0;
    this->t0 = t0;
    this->alpha = alpha;
    this->epoch_length = epoch_length;
    this->terminal_ratio = terminal_ratio;
    this->ant = instance->best_so_far_ant;
    best_length = ant->tour_length;
    
    ticks = 0;
    epoch_counter = 0;
    
    test_cnt = 0;
    improvement_cnt = 0;
    accept_cnt = 0;
    
    neighbour_search = new NeighbourSearch(instance);
    local_search = new LocalSearch(instance);
}

SimulatedAnnealing::~SimulatedAnnealing()
{
    delete neighbour_search;
}

void SimulatedAnnealing::run(void)
{
    printf("\n\n-----starting Simulated Annealing. pid: %d iter: %ld-----\n", instance->pid, instance->iteration);
    write_anneal_report(instance, NULL);
    
    while (t > (t0 / terminal_ratio)) {
        step();
    }
    
    printf("-----end Simulated Annealing. pid: %d iter: %ld-----\n\n", instance->pid, instance->iteration);
    
}


/*
 * 单步退火算法
 */
bool SimulatedAnnealing::step(void)
{
    bool accepted = false;
    bool is_valid;    /* mark if current move is valid */
    
    Move *move = neighbour_search->search(ant);
    
    if (move == NULL) {
        return false;
    }
    ticks++;
    
    is_valid = move->valid;
    if (!is_valid) {
        return false;
    }
    
    if(acceptable(move)) {
        accept(move);
        accepted = true;
        write_anneal_report(instance, move);
        
        DEBUG(check_solution(instance, ant->tour, ant->tour_size);)
    } else {
        accepted = false;
        reject(move);
    }
    
    epoch_counter++;
    if (epoch_counter >= epoch_length) {
        epoch_counter = 0;
        t = t * alpha;
        
        float ar = accept_cnt / (float) test_cnt;
        float ir = improvement_cnt / (float) test_cnt;
        printf("Time: %f, T: %f, ar: %f, ir: %f moves:%ld\n", elapsed_time(VIRTUAL), t, ar, ir, test_cnt);
        test_cnt = accept_cnt = improvement_cnt = 0;
    }
    
    delete move;
    return accepted;
}

/**
 * Accept this move
 */
void SimulatedAnnealing::accept(Move *move)
{
    // apply this neighbourhood move
    long int xxx = ant->tour_length;
    move->apply();
    
    local_search->do_local_search(ant);
    
    if (ant->tour_length < best_length) {
        best_length = ant->tour_length;
    }
}

/*
 * reject this move
 */
void SimulatedAnnealing::reject(Move *move)
{
    
}

bool SimulatedAnnealing::acceptable(Move *move)
{
    bool rv;
    test_cnt++;
    long int delta = move->gain;
    
    if (delta < 0) {
        rv = true;
        improvement_cnt++;
//        printf("Time: %f, T: %f, improvement: %ld\n", elapsed_time(VIRTUAL), t, delta);
    } else {
        rv = ran01(&instance->rnd_seed) < exp(-delta / t);
    }
    
    if (rv){
        accept_cnt++;
    }
    
    return rv;
}
