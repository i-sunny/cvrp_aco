//
//  move.hpp
//  cvrp_aco
//
//  Created by 孙晓奇 on 2016/10/19.
//  Copyright © 2016年 xiaoqi.sxq. All rights reserved.
//

#ifndef move_h
#define move_h

#include <stdio.h>
#include "problem.h"

class Move {
public:
    Move(AntStruct *ant_, bool valid_, long int gain_):ant(ant_), valid(valid_), gain(gain_){}
    virtual ~Move(){}
    virtual void apply() = 0;
    
    AntStruct *ant;
    bool valid;
    long int gain;
};

class ExchangeMove : public Move {
public:
    ExchangeMove(AntStruct *ant_, bool valid_, long int gain_, long int pos_n1_, long int pos_n2_,
                 long int load_r1_, long int load_r2_)
    :Move(ant_, valid_, gain_), pos_n1(pos_n1_), pos_n2(pos_n2_), load_r1(load_r1_), load_r2(load_r2_){}
    ~ExchangeMove(){}
    
    void apply();

    long int pos_n1;          /* 两个互换位置的node idx */
    long int pos_n2;
    long int load_r1;         /* new load of route 1 */
    long int load_r2;         /* new load of route 2 */
};


class InsertionMove : public Move {
public:
    InsertionMove(AntStruct *ant_, bool valid_, long int gain_, long int pos_n1_, long int pos_n2_,
                 long int load_r1_, long int load_r2_)
    :Move(ant_, valid_, gain_), pos_n1(pos_n1_), pos_n2(pos_n2_), load_r1(load_r1_), load_r2(load_r2_){}
    ~InsertionMove(){}
    
    void apply();

    long int pos_n1;
    long int pos_n2;
    long int load_r1;         /* new load of route 1 */
    long int load_r2;         /* new load of route 2 */
};

class InversionMove : public Move {
public:
    InversionMove(AntStruct *ant_, bool valid_, long int gain_, long int pos_n1_, long int pos_n2_)
    :Move(ant_, valid_, gain_), pos_n1(pos_n1_), pos_n2(pos_n2_){}
    ~InversionMove(){}
    
    void apply();

    long int pos_n1;
    long int pos_n2;
};

#endif /* move_h */
