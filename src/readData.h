#include "stdarg.h"
#include "utils.h"
#include "date.h"

using namespace std;

Node** tree2data(istream& tree,Pr* pr,int & s);
void readDateFile(istream &dateFile, Pr* pr,Node** & nodes,bool& constraintConsistent);
void readPartitionFile(Pr* pr);
int tree2dataS(FILE *,Pr*,Node**);
//int extrait_outgroup(string,string,list<string>&,int);
void extrait_outgroup(InputOutputStream *io, Pr* pr,list<string> &outgroups);
