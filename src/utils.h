
#ifndef UTILS_H_
#define UTILS_H_

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <stack>
#include <list>
#include <cstdlib>
#include <cmath>
#include "stdarg.h"
#include <algorithm>
#include "node.h"
#include "pr.h"
#include "pair.h"
#include "subtree.h"
#include "part.h"
#include "lsd.h"
#include <string.h>

using namespace std;
using namespace lsd;

/**
 input/output files for stand-alone program
 */
class InputOutputFile : public InputOutputStream {
public:
    /**
     constructor
     @param opt input option
     */
    InputOutputFile(Pr *opt);
    
    /** destructor */
    virtual ~InputOutputFile();

    /**
     set the content of the tree stream
     @param str a tree string
     */
    virtual void setTree(string str);
    
    /**
     set the content of the bootstrap tree stream
     @param str a tree string
     */
    virtual void setBootstrapTree(string str);

protected:
    /** true tree is a file, false other (e.g. stringstream) */
    bool treeIsFile;
    bool bootstrapTreeIsFile;
};


string readWord(istream& f,string fn);

string readWord(string line,int& pos);

char readChar(istream& f,string fn);

double readDouble(string line,int& pos);

double readdouble(istream& f,string fn);

string realToYearMonth(double y);

string realToYearMonthDay(double y);

double monthToReal(int m);

double monthDayToReal(int m,int d);

double readDate(istream& f,string fn,Pr* pr,double& month, double& day);

double readDate1(istream& f,string fn,char c,Pr* pr,double& month, double& day);

bool readDateFromString(const char* st,double& f);

void readWholeDate(istream &dateFile,Pr* pr,int& type,double& v1,double& v2, double& m1,double& m2,double& d1,double& d2,int& dateFormat);

vector<double> read_double_from_line(string line);

int readInt(istream& f,string msg);

int getLineNumber(istream &f);

string readSupport(istream& f,string fn);

void concat(list<int> & l1,list<int> l2);

void concatPos(list<int> l1,list<int> &l2);

void concat(stack<int> & l1,list<int> l2);

void concat(list<int> & l1,stack<int> l2);

int getPosition(Node**,string s,int n,int m);

list<int> intersect(list<int> l1,list<int> l2);

bool isAncestor(Node**,int i,int j);

int mrca(Node**,int i,int j);

int index(list<int> & L, int e);

int index(string s,string* & L,int n);

bool contain(int s,vector<int> l);

bool contain(int s,list<int> l);

bool contain(string s,list<string> l);

//string readLabel(istream& f,FILE *w);

string readLabel(char ch,istream& f,int& a);


char readBracket(istream& f,string fn);

char readCommaBracket(istream& f,string fn,string& s);

char read2P(istream& f,string fn);

void computeSuc(int* & Pre,int* & Suc1,int* & Suc2,int size,int n);

void computeVariance(Pr* pr,Node** nodes);

void computeNewVariance(Pr* pr,Node** nodes);

void myExit( string msg, ... );


void myErrorMsg( string msg, ... );

bool isReal( const char* str );

bool isInteger( const char* str );

void sort(int* & tab,int size);

void sort(double* & tab,int size);
 
int index(int* & tab,int value,int size);

int mrca(Node**,list<int> taxa);

int mrca(Node**,vector<int> taxa);

bool markLeaf(Node* no);

bool leaf(Node* no);

bool tc(Node* no);

bool limit(Node* no);

bool tc_limit(Node* no);

bool lower(Node* no);

bool upper(Node* no);

void activeMarkLeaf(Node* no);

void desactiveMarkLeaf(Node* no);

void activeTC(Node* no);

void activeUpper(Node* no);

void activeLower(Node* no);

void desactive(Node* no);

void desactiveTC(Node* no);

void desactiveLimit(Node* no);

list<int> down_leaf(int i,Pr* pr,Node** nodes);

stack<int>* computeFeuilles(list<int> ls,Pr* pr,Node** nodes);

list<int> suc(int i,int j,Pr* pr,Node** nodes,int* & Pre,list<int> &sucL,list<int> &sucI);

bool reroot_rootedtree(double& br,int r,int s1, int s2,Pr* pr,Node** nodes,Node** &nodes_new);

bool reroot_rootedtree(double& br,int r,Pr* pr,Node** nodes);

int reroot_rootedtree(int r,Pr* pr,Node** nodes,Node** &nodes_new);

bool reroot_rootedtree(double& br,int r,int s1, int s2,Pr* pr,Node** nodes,Node** &nodes_new,int* & P_ref,int* & tab);

void unrooted2rooted(Pr* &pr,Node** nodes);

void rooted2unrooted(Pr* &pr,Node** nodes);

Node** unrooted2rootedS(Pr* &pr,Node** nodes,int s);

void computeObjective(Pr* pr,Node** nodes);

void computeObjectiveMultiRates(Pr* pr,Node** nodes);

void computeObjectiveMultiRates(Pr* pr,Node** nodes,double* B, double* V);

void computeObjectiveEstimateRoot(int r,int p_r,double br,Pr* pr,Node** nodes);

string newick(int i,int terminate,Pr* pr,Node** nodes,int& nbTips);

string nexus(int i,Pr* pr,Node** nodes);

string nexusDate(int i,Pr* pr,Node** nodes);

string nexusIC(int i,Pr* pr,Node** nodes,double* D_min,double* D_max,double* H_min,double* H_max);

string nexusICDate(int i,Pr* pr,Node** nodes,double* D_min,double* D_max,double* H_min,double* H_max);


double* variance(int w,int m,double* B,int c,int s);

double variance(Pr* pr,double b);

double* newvariance(int m,double rho,int* P,double* T,int c,int s);

list<string> getOutgroup(istream &o, string fn);

bool initConstraintReRooted(Pr* pr,Node** nodes);

bool initConstraint(Pr* pr,Node** nodes);

void initialize_status(Pr* &pr,Node** &nodes);

list<int> getActiveSet(Pr* pr,Node** nodes);

Node** cloneLeaves(Pr* pr,Node** nodes,int f);

void cloneInternalNodes(Pr* pr,Node** nodes,Node** & nodes_new,int f);

void computeSuc_polytomy(Pr* pr,Node** nodes);

list<int> pos_polytomy(int i,Pr* pr,Node** nodes);

list<int> postorder_polytomy(Pr* pr,Node** nodes);

vector<int> pre_polytomy(int i,Pr* pr,Node** nodes);

vector<int> preorder_polytomy(Pr* pr,Node** nodes);

vector<int> pre_polytomy_withTips(int i,Pr* pr,Node** nodes);

vector<int> preorder_polytomy_withTips(Pr* pr,Node** nodes);

stack<int>* computeFeuilles_polytomy(list<int> ls,Pr* pr,Node** nodes);

list<int> down_polytomy(int,int,int* &,list<int>* &,list<int>* &,bool* &,int* &,double* &);

void reduceTree_polytomy(Pr* pr,Node** nodes,int* &Pre,list<int>* & Suc,double* & C,list<int>* &internal);

list<int> suc_polytomy(int i,int j,Pr* pr,Node** nodes,int* & Pre,double* & C,list<int> &suc);

void computeVarianceEstimateRoot(Pr* pr,Node** nodes,double br);

void computeNewVarianceEstimateRoot(Pr* pr,Node** nodes);

void checkRooted(Pr* & opt);

int getInternalNodeId(Pr* pr,Node** nodes,string s);

int firstCharacter(string s,int p);

int lastCharacter(string s,int p);

int assignRecursive(int r,Node** nodes,int g);

int assignRateGroupToSubTree(Subtree* subtree,Pr* pr,Node** nodes,int g);

void assignRateGroupToTree(Pr* pr,Node** nodes);

void calculateOulier(double* sortedArray,double& mi,double& ma,int size,double k);

double* calculateRtt(Pr* pr,Node** nodes);

double* sortTab(double* tab,int size);

int maxDate(int month);

void calculate_tree_height(Pr* pr,Node** & nodes);

void splitExternalBranches(Pr* pr,Node** nodes);

void splitLongBranches(Pr* pr,Node** nodes,double th);

bool isIn(int i,vector<int> v);

vector<int> intersect(vector<int> v1, vector<int> v2);

bool checkAllConstraintConsistent(Pr* pr,Node** nodes);

double median(vector<double> array);
    
double median_branch_lengths(Pr* pr,Node** nodes);

int collapseTree(Pr* pr,Node** nodes,Node** nodes_new,int* &tab, double toCollapse,bool& useSupport);
    
void collapse(int i,int j,Pr* pr,Node** nodes,Node** nodes_new,int &cc,int* &tab, double toCollapse, bool useSupport, double* support);

void collapseTreeReOrder(Pr* pr,Node** nodes,Pr* prReduced,Node** nodesReduced,int* &tab);

void collapseUnInformativeBranches(Pr* &pr,Node** &nodes,bool verbose);

void adjustNodeDateToYMD(Node*& node,int m1,int d1,int m2,int d2);

void adjustNodeDateToYM(Node*& node,int m1,int d1,int m2,int d2);

void adjustDateToYMD(Date*& date,int m1,int d1,int m2,int d2);

void adjustDateToYM(Date*& date,int m1,int d1,int m2,int d2);

bool checkTopology(Pr* pr,Node** nodes1, Node** nodes2);

double* rtt(Pr* pr,Node** nodes);

void starting_pointLower(Pr* pr,Node** nodes,list<int> & active_set);

void starting_pointUpper(Pr* pr,Node** nodes,list<int> & active_set);
#endif
