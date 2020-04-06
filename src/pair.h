#ifndef PAIR_H
#define PAIR_H
#include "stdarg.h"

using namespace std;
typedef struct Pair
{
    bool include;
    string name;
    Pair(bool inc,string n){
        include = inc;
        name = n;
    }
} Pair;
#endif



