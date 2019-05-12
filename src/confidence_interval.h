#include <random>
#include "dating.h"
#include "estimate_root.h"
#include "utils.h"

using namespace std;


int collapseTree(Pr* pr,Node** nodes,int* &tab,Pr* prReduced,Node** nodesReduced);

void collapse(int i,int j,Pr* pr,Node** nodes,Node** nodes_new,int &cc,int* &tab);

void collapseTreeReOrder(int cc,int* &tab,Pr* pr,Node** nodes,Pr* prReduced,Node** nodesReduced);

void computeIC(double br,Pr* pr,Node** nodes,double* &T_left,double* &T_right,double* &H_left,double* &H_right,double* &HD_left,double* &HD_right,double &rho_left,double& rho_right,double* &other_rhos_left,double* &other_rhos_right);


void output(double br,int y, Pr* pr,Node** nodes,FILE* f,FILE* tree1,FILE* tree2);
