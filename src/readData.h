#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <stack>
#include <list>
#include "stdarg.h"
#include <math.h>
#include "utils.h"
#include "date.h"

using namespace std;

Node** tree2data(FILE * tree,Pr* pr,int & s);
void readDateFile(Pr* pr,Node** & nodes,bool& constraintConsistent);
void readPartitionFile(Pr* pr);
int tree2dataS(FILE *,Pr*,Node**);
//int extrait_outgroup(string,string,list<string>&,int);
void extrait_outgroup(Pr* pr,list<string> &outgroups);
