
#include <stdio.h>
#include "utils.h"

using namespace std;

bool calculateOutliers(Pr* & pr,Node** & nodes);
void estimate_root_rtt(Pr* pr, Node** & nodes);
void estimate_root_local_rtt(Pr* pr, Node** & nodes);
double regression_lambda(double br,double &lambda,Pr* pr, Node** nodes);
double* residus_rtt(vector<double> paths,vector<double> dates,double slope,double intercept);
void outlier(Pr* pr,Node** nodes,double k);
bool remove_one_tip(Pr* pr,Node** nodes,int t,int* &tab);
void shift_node_id(Pr* &pr,Node** &nodes,int* &tab);
bool remove_outlier_tips(Pr* &pr,Node** &nodes);
void getOulier(double* sortedArray,double& mi,double& ma,int size,double k);
void calculateRtt(Pr* pr,Node** nodes,vector<double> &paths,vector<double> &dates);
void calculateRtt_lambda(double br, Pr* pr,Node** nodes,vector<double> & paths,vector<double> &paths_lambda,vector<double> &dates);
//void calculateRtt_lambda(double br,Pr* pr,Node** nodes,double* & rtt, double* & rtt_lambda);
void regression(Pr* pr,vector<double> paths,vector<double> dates,double & slope,double & intercept);
