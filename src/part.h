#ifndef PART_H
#define PART_H
#include "stdarg.h"
#include "subtree.h"

using namespace std;
typedef struct Part
{
    string name;
    vector<Subtree*> subtrees;
    Part(string s){
        name =s;
    }
} Part;
#endif



