
#include <stdio.h>
#include "utils.h"

using namespace std;

bool calculateOutliers(Pr* & pr,Node** & nodes);
bool remove_one_tip(Pr* pr,Node** nodes,int t,int* &tab);
void shift_node_id(Pr* &pr,Node** &nodes,int* &tab);
bool remove_outlier_nodes(Pr* &pr,Node** &nodes);
vector<int> sampleNoRepeat(vector<int> array,int ex, int size);
vector<double> residus_lsd(Pr* pr,Node** nodes,double& mean_res, double& var_res);
double median_rate(Pr* pr,Node** nodes, vector<double> dates, vector<int>* samples,bool addInternalDates);
void calculateRoot2DatedNode(Pr* pr,Node** nodes,vector<double> &paths);
void calculateRoot2DatedNode(Pr* pr,Node** nodes,vector<double> &paths,vector<double> &dates);
vector<int> outliers_rooted(Pr* pr,Node** nodes,vector<int>* samples, vector<double> dates,bool addInternalDates);
bool outliers_unrooted(Pr* &pr,Node** &nodes);
