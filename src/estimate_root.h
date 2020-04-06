#include "utils.h"
#include "dating.h"


using namespace std;

bool without_constraint_lambda(double br,Pr* &par,Node** &nodes,list<int> active_set, list<double> & ld);

bool starting_point_without_constraint_lambda(double br,Pr* &pr,Node** &nodes,list<int> &active_set);

bool without_constraint_active_set_lambda(double br,Pr* &pr,Node** &nodes);

bool with_constraint_lambda(double br,Pr* &pr,Node** &nodes,list<int> active_set,list<double> & ld);


bool with_constraint_active_set_lambda(double br,Pr* &pr,Node** &nodes);

int estimate_root_without_constraint_local_rooted(Pr* &pr,Node** &nodes);

int estimate_root_without_constraint_rooted(Pr* &pr,Node** &nodes);

int estimate_root_with_constraint_local_rooted(Pr* &pr,Node** &nodes);

int estimate_root_with_constraint_fast_rooted(Pr* &pr,Node** &nodes);

int estimate_root_with_constraint_rooted(Pr* &pr,Node** &nodes);

void calculateMultiplier_lambda(int r,int p_r,double br,Pr* &pr,Node** &nodes,bool* nan);

bool without_constraint_active_set_lambda_multirates(double br,Pr* &pr,Node** &nodes,bool reassign);

bool with_constraint_active_set_lambda_multirates(double br,Pr* &pr,Node** &nodes,bool reassign);

bool without_constraint_active_set_lambda_multirates_exclude_outlier(double br,Pr* &pr,Node** &nodes,bool reassign);

bool with_constraint_active_set_lambda_multirates_exclude_outlier(double br,Pr* &pr,Node** &nodes,bool reassign);

double regression_lambda(double br,double &lambda,Pr* pr, Node** nodes);

void estimate_root_rtt(Pr* pr, Node** & nodes);
