#include <random>
#include "dating.h"
#include "estimate_root.h"
#include "utils.h"

using namespace std;

void computeIC(double br,Pr* pr,Node** nodes,double* &T_left,double* &T_right,double* &H_left,double* &H_right,double* &HD_left,double* &HD_right,double &rho_left,double& rho_right,double* &other_rhos_left,double* &other_rhos_right);


void output(double br,int y, Pr* pr,Node** nodes,ostream& f,ostream& tree1,ostream& tree2,ostream& tree3);
