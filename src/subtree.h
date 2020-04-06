#ifndef SUBTREE_H
#define SUBTREE_H
#include "stdarg.h"
#include "pair.h"

using namespace std;
typedef struct Subtree
{
    Pair* root;
    vector<Pair*> tips;
    Subtree(Pair* r){
        root = r;
    }
    Subtree(Pair* r,vector<Pair*> t){
        root = r;
        tips = t;
    }
} Subtree;
#endif



