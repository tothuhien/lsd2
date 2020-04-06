#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <stack>
#include <list>
#include "stdarg.h"
#include <math.h>
#include "utils.h"

using namespace std;

bool without_constraint(Pr* pr,Node** nodes);


bool conditions(list<double>& ldLagrange,Pr* pr,Node** nodes);

void starting_point(Pr* pr,Node** nodes,list<int> & active_set);

list<double> computeLambda(list<int> active_set,Pr* pr,Node** nodes);

bool remove_ne_lambda(list<double> & lambda,list<int> & active_set,int& as);


bool without_constraint_active_set(Pr* pr,Node** nodes);

bool conditionsQP(list<double>& ldLagrange,Pr* pr,Node** nodes);

//bool temporalConstraintsConsistent(Pr* pr,Node** nodes);

bool starting_pointQP(Pr* pr,Node** nodes,list<int>& active_set);

bool with_constraint(Pr* pr,Node** &nodes,list<int> active_set,list<double>& lambda);

bool with_constraint_active_set(Pr* pr,Node** &nodes);

void calculateMultiplier(Pr* pr,Node** nodes);

bool without_constraint_multirates(Pr* pr,Node** nodes,bool reassign);

bool with_constraint_multirates(Pr* pr,Node** nodes,bool reassign);

