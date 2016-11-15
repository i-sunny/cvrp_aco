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

enum MoveType{
    EXCHANGE_MOVE, INSERTION_MOVE, INVERSION_MOVE
};

class Move {
public:
    Move(AntStruct *ant_, bool valid_, double gain_, int pos_n1_, int pos_n2_, MoveType type_)
    :ant(ant_), valid(valid_), gain(gain_), pos_n1(pos_n1_), pos_n2(pos_n2_), type(type_){}
    virtual ~Move(){}
    virtual void apply(){}
    
    AntStruct *ant;
    bool valid;
    double   gain;
    int pos_n1;          /* 两个互换位置的node idx */
    int pos_n2;
    
    MoveType type;            /* move type */
};

class ExchangeMove : public Move {
public:
    ExchangeMove(AntStruct *ant_, bool valid_, double gain_, int pos_n1_, int pos_n2_,
                 int load_r1_, int load_r2_)
    :Move(ant_, valid_, gain_, pos_n1_, pos_n2_, EXCHANGE_MOVE), load_r1(load_r1_), load_r2(load_r2_){}
    ~ExchangeMove(){}
    
    void apply();
    
    int load_r1;         /* new load of route 1 */
    int load_r2;         /* new load of route 2 */
};


class InsertionMove : public Move {
public:
    InsertionMove(AntStruct *ant_, bool valid_, double gain_, int pos_n1_, int pos_n2_,
                 int load_r1_, int load_r2_)
    :Move(ant_, valid_, gain_, pos_n1_, pos_n2_, INSERTION_MOVE), load_r1(load_r1_), load_r2(load_r2_){}
    ~InsertionMove(){}
    
    void apply();

    int load_r1;         /* new load of route 1 */
    int load_r2;         /* new load of route 2 */
};

class InversionMove : public Move {
public:
    InversionMove(AntStruct *ant_, bool valid_, double gain_, int pos_n1_, int pos_n2_)
    :Move(ant_, valid_, gain_, pos_n1_, pos_n2_, INVERSION_MOVE){}
    ~InversionMove(){}
    
    void apply();
};

#endif /* move_h */
