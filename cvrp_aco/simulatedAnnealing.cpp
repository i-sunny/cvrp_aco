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
#include <assert.h>

#include "simulatedAnnealing.h"
#include "utilities.h"
#include "timer.h"
#include "io.h"

bool tabu_flag = true;

SimulatedAnnealing::SimulatedAnnealing(Problem *instance, AntColony *ant_colony, double t0,
                                       double alpha, int epoch_length, int terminal_ratio)
{
    this->instance = instance;
    this->t = t0;
    this->t0 = t0;
    this->alpha = alpha;
    this->epoch_length = epoch_length;
    this->terminal_ratio = terminal_ratio;
    
    // copy best so far ant
    best_ant = new AntStruct();
    best_ant->tour = new int[2*instance->num_node-1];
    AntColony::copy_solution_from_to(instance->best_so_far_ant, best_ant);
    
    iter_ant = new AntStruct();
    iter_ant->tour = new int[2*instance->num_node-1];
    AntColony::copy_solution_from_to(best_ant, iter_ant);
    
    iteration = 0;
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
    printf("\n----- Start SA. length: %f iter: %d time: %f-----\n", best_ant->tour_length, instance->iteration, elapsed_time( VIRTUAL));
    write_anneal_report(instance, iter_ant, NULL);
    
    tabu_list.clear();
    
    while (t > (t0 / terminal_ratio)) {
        step();
    }
    
    if (best_ant->tour_length - instance->best_so_far_ant->tour_length < -EPSILON) {
        AntColony::copy_solution_from_to(best_ant, instance->best_so_far_ant);
        if (instance->pid == 0) {
            instance->best_so_far_time = elapsed_time(VIRTUAL);
            write_best_so_far_report(instance);
        }
    }
    
    ant_colony->compute_total_information();
    
    printf("----- End SA. length: %f iter: %d time: %f-----\n", best_ant->tour_length, instance->iteration, elapsed_time( VIRTUAL));
    
}


/*
 * 单步退火算法
 */
bool SimulatedAnnealing::step(void)
{
    bool accepted = false;
    
    Move *move = neighbour_search->search(iter_ant);

    if (move != NULL) {
        iteration++;
        if (move->valid) {
            // valid move
            if(acceptable(move)) {
                accept(move);
                accepted = true;
                
                if (tabu_flag) {
                    update_tabu_list(move);
                }
//                write_anneal_report(instance, iter_ant, move);
                
                DEBUG(assert(check_solution(instance, iter_ant->tour, iter_ant->tour_size));)
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
                DEBUG(printf("Time: %f, T: %f, ar: %f, ir: %f moves:%ld\n", elapsed_time(VIRTUAL), t, ar, ir, test_cnt);)
                test_cnt = accept_cnt = improvement_cnt = 0;
            }
            
            // update pheromone
        //    if (move->gain < -EPSILON) {
        //        ant_colony->global_update_pheromone_weighted(iter_ant, 0.01);
        //    }
        }
    }
    delete move;
    return accepted;
}

bool SimulatedAnnealing::acceptable(Move *move)
{
    bool accepted = false;
    test_cnt++;
    double delta = move->gain;
    
    if (tabu_flag) {
        // this move is in tabu list
        if (is_tabu(move)) {
            return false;
        }
    }
    
    if (delta < -EPSILON) {
        accepted = true;
        improvement_cnt++;
        //        printf("Time: %f, T: %f, improvement: %ld\n", elapsed_time(VIRTUAL), t, delta);
    } else if (fabs(delta) < EPSILON) {
        accepted = true;
    }else {
        accepted = ran01(&instance->rnd_seed) < exp(-delta / t);
    }
    
    if (accepted){
        accept_cnt++;
    }
    
    return accepted;
}

/**
 * Accept this move
 */
void SimulatedAnnealing::accept(Move *move)
{
    // apply this neighbourhood move
    move->apply();
    
    // seem to have worse performance
//    local_search->do_local_search(iter_ant);
    
    if (iter_ant->tour_length - best_ant->tour_length < -EPSILON) {
        AntColony::copy_solution_from_to(iter_ant, best_ant);
        // update pheromone
        ant_colony->global_update_pheromone_weighted(iter_ant, 2 * ras_ranks);
        printf("[%d]SA better solution. length:%f, sa_iter:%d\n", instance->pid,iter_ant->tour_length, iteration);
    }
}

/*
 * reject this move
 */
void SimulatedAnnealing::reject(Move *move)
{
    
}


bool SimulatedAnnealing::is_tabu(Move *move)
{
    Move *tabu_move;
    int i;
    int pos_n1, pos_n2;
    
    pos_n1 = move->pos_n1;
    pos_n2 = move->pos_n2;
    if (move->type == INSERTION_MOVE) {
        pos_n1 = move->pos_n2;
        pos_n2 = move->pos_n1;
    }
    
    for (i = 0; i < tabu_list.size(); i++) {
        tabu_move = &tabu_list[i].move;
        if (tabu_move->type == move->type &&
            tabu_move->pos_n1 == pos_n1 &&
            tabu_move->pos_n2 == pos_n2 &&
            fabs(tabu_move->gain + move->gain) < EPSILON &&
            tabu_list[i].life > iteration)
        {
            return true;
        }
    }
    return false;
}

void SimulatedAnnealing::update_tabu_list(Move *move)
{
    Tabu tabu(*move, iteration + TABU_LENGTH);
    tabu_list.push_back(tabu);
    
    vector<Tabu>::iterator ibeg = tabu_list.begin();
    while (ibeg->life <= iteration) {
        tabu_list.erase(ibeg);
    }
}
