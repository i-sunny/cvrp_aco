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

SimulatedAnnealing::SimulatedAnnealing(Problem *instance, AntColony *ant_colony, double t0,
                                       double alpha, long int epoch_length, long int terminal_ratio)
{
    this->instance = instance;
    this->t = t0;
    this->t0 = t0;
    this->alpha = alpha;
    this->epoch_length = epoch_length;
    this->terminal_ratio = terminal_ratio;
    
    // copy best so far ant
    best_ant = new AntStruct();
    best_ant->tour = new long int[2*instance->num_node-1];
    AntColony::copy_solution_from_to(instance->best_so_far_ant, best_ant);
    
    iter_ant = new AntStruct();
    iter_ant->tour = new long int[2*instance->num_node-1];
    AntColony::copy_solution_from_to(best_ant, iter_ant);
    
    ticks = 0;
    epoch_counter = 0;
    
    test_cnt = 0;
    improvement_cnt = 0;
    accept_cnt = 0;
    
    neighbour_search = new NeighbourSearch(instance);
    local_search = new LocalSearch(instance);
    this->ant_colony = ant_colony;
}

SimulatedAnnealing::~SimulatedAnnealing()
{
    delete neighbour_search;
    delete local_search;
    
    delete iter_ant->tour;
    delete best_ant->tour;
    
    delete iter_ant;
    delete best_ant;
}

void SimulatedAnnealing::run(void)
{
    printf("\n\n-----starting Simulated Annealing. pid: %d iter: %ld-----\n", instance->pid, instance->iteration);
    write_anneal_report(instance, iter_ant, NULL);
    
    while (t > (t0 / terminal_ratio)) {
        step();
    }
    
    if (best_ant->tour_length < instance->best_so_far_ant->tour_length) {
        AntColony::copy_solution_from_to(best_ant, instance->best_so_far_ant);
        write_best_so_far_report(instance);
    }
    
    ant_colony->compute_total_information();
    
    printf("-----end Simulated Annealing. pid: %d iter: %ld-----\n\n", instance->pid, instance->iteration);
    
}


/*
 * 单步退火算法
 */
bool SimulatedAnnealing::step(void)
{
    bool accepted = false;
    bool is_valid;    /* mark if current move is valid */
    
    Move *move = neighbour_search->search(iter_ant);

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
        write_anneal_report(instance, iter_ant, move);
        
        DEBUG(check_solution(instance, iter_ant->tour, iter_ant->tour_size);)
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
    
    // update pheromone
//    if (move->gain < 0) {
//        ant_colony->global_update_pheromone_weighted(iter_ant, 0.01);
//    }
    
    delete move;
    return accepted;
}

/**
 * Accept this move
 */
void SimulatedAnnealing::accept(Move *move)
{
//    print_solution(instance, iter_ant->tour, iter_ant->tour_size);
    
    // apply this neighbourhood move
    move->apply();
    
    // seem to have worse performance
//    local_search->do_local_search(iter_ant);
    
    if (iter_ant->tour_length < best_ant->tour_length) {
        AntColony::copy_solution_from_to(iter_ant, best_ant);
        // update pheromone
        ant_colony->global_update_pheromone_weighted(iter_ant, 2 * ras_ranks);
    }
    
//    print_solution(instance, iter_ant->tour, iter_ant->tour_size);
}

/*
 * reject this move
 */
void SimulatedAnnealing::reject(Move *move)
{
    
}

bool SimulatedAnnealing::acceptable(Move *move)
{
    bool accepted = false;
    test_cnt++;
    long int delta = move->gain;
    
    if (delta < 0) {
        accepted = true;
        improvement_cnt++;
//        printf("Time: %f, T: %f, improvement: %ld\n", elapsed_time(VIRTUAL), t, delta);
    } else {
        accepted = ran01(&instance->rnd_seed) < exp(-delta / t);
    }
    
    if (accepted){
        accept_cnt++;
    }
    
    return accepted;
}
