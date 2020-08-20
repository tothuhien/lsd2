#include "stdarg.h"
#include "utils.h"
#include "date.h"
#include "lsd.h"

using namespace std;

Node** tree2data(istream& tree,Pr* pr,int & s);
void readInputDate(InputOutputStream* io, Pr* pr,Node** &nodes,bool& constraintConsistent);
void readPartitionFile(istream &partFile, Pr* pr);
int tree2dataS(FILE *,Pr*,Node**);
void extrait_outgroup(InputOutputStream *io, Pr* pr, bool useBootstrapTree);
