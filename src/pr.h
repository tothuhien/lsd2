#ifndef PR_H
#define PR_H
#include "stdarg.h"
#include <vector>
#include "date.h"
#include "pair.h"
#include "part.h"
#include <cmath>
#include "subtree.h"

using namespace std;

typedef struct Pr
{
    string  inFile;         //Nom du fichier d'arbres
    string  inDateFile;     //Nom du fichier contenant les dates
    string partitionFile;   //Nom du fichier contenant la partition des arretes pour different taux
    string bootstraps_file; //Name of the file containing bootstrap trees
    string  outFile;        //Nom du fichier de resultats.
    string treeFile2;       //Nom du fichier d'arbres sorties Nexus avec des longueurs de branches mesures par des temps ecoules
    string treeFile3;       //Nom du fichier d'arbres sorties Newick
    //bool relative;         //=true if all the leaves have the same date, the program estimate the relative dates
    double mrca;
    string MRCA;
    double leaves;
    string LEAVES;
    int    seqLength;      //Longueur des sequences dans l'alignement
    int    nbData;         //Nombre de cas a  traiter
    int    nbBootstrap;     //Number of bootstrap trees
    string fnOutgroup;
    string rate;           //le fichier contient les taux en entree
    string estimate_root;    //Method to estimate root
    bool rooted;
    bool constraint;       //Impose the constraints or not
    int variance;         //Use the variances or not
    bool ci;         //Compute confidence interval or not
    double c;                //either -b is specified or not
    double b;                 //var = b+c/s;
    double q; //standard deviation of lognormal distribution to simulate relaxed branch lengths for calculating confidence intervals
    double nullblen;
    double support;
    bool verbose;
    double rho_min;
    int nbINodes;
    int nbBranches;
    double minblen;
    double minblenL;
    double rho;
    double round_time;
    int inDateFormat;//3 for year-month; 2 for year-month-day; 1 for real from 9-9999; 0 for everything else
    int outDateFormat;//3 for year-month; 2 for year-month-day; 1 for real
    double rhoLower;
    double rhoUpper;
    bool haveLower;
    bool haveUpper;
    bool haveUnique;
    vector<double> multiplierRate;
    vector<bool> givenRate;
    double objective;
    int nbSampling;
    bool removeOutgroup;
    int m;
    double e;
    bool splitExternal;
    vector<Part* > ratePartition;
    vector<Date*> internalConstraints;
    vector<int> outlier;
    vector<string> warningMessage;
    vector<string> resultMessage;
    Pr(int n,int m){
        nbINodes=n;
        nbBranches=m;
    }
    void init(){
        internalConstraints.clear();
        rooted = false;
        if (fnOutgroup != "" ) rooted = true;
        warningMessage.clear();
        resultMessage.clear();
        outlier.clear();
    }
    void copy(Pr* pr){
        inFile = pr->inFile;
        inDateFile = pr->inDateFile;
        partitionFile = pr->partitionFile;
        bootstraps_file = pr->bootstraps_file;
        fnOutgroup = pr->fnOutgroup;
        outFile = pr->outFile;
        treeFile2 = pr->treeFile2;
        outlier = pr->outlier;
        mrca=pr->mrca;
        MRCA=pr->MRCA;
        leaves=pr->leaves;
        LEAVES=pr->LEAVES;
        seqLength=pr->seqLength;
        nbData=pr->nbData;
        nbBootstrap = pr->nbBootstrap;
        rate=pr->rate;
        estimate_root = pr->estimate_root;
        rooted=pr->rooted;
        constraint=pr->constraint;
        variance=pr->variance;
        ci=pr->ci;
        q=pr->q;
        c=pr->c;
        b=pr->b;
        e=pr->e;
        m=pr->m;
        round_time=pr->round_time;
        support=pr->support;
        minblen=pr->minblen;
        minblenL=pr->minblenL;
        nullblen=pr->nullblen;
        verbose=pr->verbose;
        splitExternal=pr->splitExternal;
        inDateFormat=pr->inDateFormat;
        outDateFormat=pr->outDateFormat;
        rho=pr->rho;
        rho_min=pr->rho_min;
        givenRate=pr->givenRate;
        nbSampling=pr->nbSampling;
        ratePartition = pr->ratePartition;
        internalConstraints = pr->internalConstraints;
        multiplierRate =  vector<double>();
        if (pr->ratePartition.size()>0) {
            for (int i=0;i<pr->multiplierRate.size();i++){
                multiplierRate.push_back(pr->multiplierRate[i]);
            }
        }
        rhoLower = pr->rhoLower;
        rhoUpper = pr->rhoUpper;
        haveLower = pr->haveLower;
        haveUpper = pr->haveUpper;
        haveUnique = pr->haveUnique;
        removeOutgroup = pr->removeOutgroup;
        m = pr->m;
        warningMessage = pr->warningMessage;
        resultMessage = pr->resultMessage;
    }
    Pr()
    {
        inFile = "";
        inDateFile = "";
        partitionFile = "";
        bootstraps_file = "";
        outFile = "";
        treeFile2 = "";
        fnOutgroup = "";
        seqLength = 1000;
        nbData = 1;
        nbBootstrap = 0;
        rate = "";
        //relative = false;
        mrca=0;
        leaves=1;
        MRCA="";
        LEAVES="";
        estimate_root = "";
        constraint = true;
        variance = 1;
        minblen = 0;
        minblenL = 0;
        nullblen = 0;
        support = -1;
        c = -1;
        b = -1;
        q = 0.2;
        rho_min = 1e-10;
        ci = false;
        nbSampling=100;
        round_time=-1;
        inDateFormat=0;
        outDateFormat=0;
        rooted=false;
        removeOutgroup=false;
        splitExternal=false;
        ratePartition = vector<Part* >();
        internalConstraints = vector<Date* >();
        outlier = vector<int>();
        e = -1;
        m = 10;
        verbose = false;
        nbINodes = 0;
        nbBranches = 0;
        rho = 0;
        givenRate = vector<bool>();
        givenRate.push_back(false);
        multiplierRate = vector<double>();
        rhoLower = 0;
        rhoUpper = 0;
        haveLower = false;
        haveUpper = false;
        haveUnique = false;
        objective = 0;
        warningMessage = vector<string>();
        resultMessage = vector<string>();
    }
}Pr;
#endif


