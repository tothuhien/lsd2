#include "utils.h"
#include "dating.h"


using namespace std;

bool without_constraint_lambda(double br,Pr* &par,Node** &nodes,list<int> active_set, list<double> & ld);

bool without_constraint_active_set_lambda(double br,Pr* &pr,Node** &nodes,int whichStartingPoint);

bool with_constraint_lambda(double br,Pr* &pr,Node** &nodes,list<int> active_set,list<double> & ld);


bool with_constraint_active_set_lambda(double br,Pr* &pr,Node** &nodes,int whichStartingPoint);

int estimate_root_without_constraint_local_rooted(Pr* &pr,Node** &nodes);

int estimate_root_without_constraint_rooted(Pr* &pr,Node** &nodes);

int estimate_root_with_constraint_local_rooted(Pr* &pr,Node** &nodes);

int estimate_root_with_constraint_fast_rooted(Pr* &pr,Node** &nodes);

int estimate_root_with_constraint_rooted(Pr* &pr,Node** &nodes);

void calculateMultiplier_lambda(int r,int p_r,double br,Pr* &pr,Node** &nodes,bool* nan);

bool without_constraint_active_set_lambda_multirates(bool all, double br,Pr* &pr,Node** &nodes,bool reassign);

bool with_constraint_active_set_lambda_multirates(bool all,double br,Pr* &pr,Node** &nodes,bool reassign);

bool without_constraint_active_set_lambda_multirates_exclude_outlier(double br,Pr* &pr,Node** &nodes,bool reassign);

bool with_constraint_active_set_lambda_multirates_exclude_outlier(double br,Pr* &pr,Node** &nodes,bool reassign);

double regression_lambda(double br,double &lambda,Pr* pr, Node** nodes);

void estimate_root_rtt(Pr* pr, Node** & nodes);

/* Root the tree by root-to-tip regression, minimising the residual sum of squares over all root
 positions (TreeTime-like). If forcedBranch > 0, only that branch is searched. The tree is rerooted
 in place (nodes is replaced) and the index of the branch holding the root is returned. */
int estimate_root_rtt_residual(Pr* pr, Node** &nodes, int forcedBranch, bool verbose);

void imposeMinBlen(ostream& result,Pr* pr, Node** nodes,double minB, bool verbose);
