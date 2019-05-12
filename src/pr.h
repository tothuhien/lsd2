#ifndef PR_H
#define PR_H
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <vector>
#include "date.h"
#include "pair.h"
#include "part.h"
#include "subtree.h"

using namespace std;

typedef struct Pr
{
    string  inFile;         //Nom du fichier d'arbres
    string  inDateFile;     //Nom du fichier contenant les dates
    string partitionFile;   //Nom du fichier contenant la partition des arretes pour different taux
    string  outFile;        //Nom du fichier de resultats.
    string treeFile1;       //Nom du fichier d'abres sorties Nexus
    string treeFile2;       //Nom du fichier d'arbres sorties Nexus avec des longueurs de branches mesures par des temps ecoules
    bool relative;         //=true if all the leaves have the same date, the program estimate the relative dates
    double mrca;
    double leaves;
    int    seqLength;      //Longueur des sequences dans l'alignement
    int    nbData;         //Nombre de cas a  traiter (dans le cas de bootstrap)
    string fnOutgroup;
    string rate;           //le fichier contient les taux en entree
    string estimate_root;    //Method to estimate root
    bool rooted;
    bool constraint;       //Impose the constraints or not
    int variance;         //Use the variances or not
    bool ci;         //Compute confidence interval or not
    int  c;                //var = b+c/s;
    bool verbose;
    double rho_min;
    int nbINodes;
    int nbBranches;
    double rho;
    vector<double> multiplierRate;
    vector<bool> givenRate;
    double objective;
    int nbSampling;
    bool keepOutgroup;
    double k;
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
        warningMessage.clear();
        resultMessage.clear();
    }
    void copy(Pr* pr){
        inFile = pr->inFile;
        inDateFile = pr->inDateFile;
        partitionFile = pr->partitionFile;
        outFile = pr->outFile;
        treeFile1 = pr->treeFile1;
        treeFile2 = pr->treeFile2;
        relative=pr->relative;
        mrca=pr->mrca;
        leaves=pr->leaves;
        seqLength=pr->seqLength; 
        nbData=pr->nbData;
        rate=pr->rate;
        estimate_root = pr->estimate_root;
        rooted=pr->rooted;
        constraint=pr->constraint;  
        variance=pr->variance;  
        ci=pr->ci;     
        c=pr->c;
        verbose=pr->verbose;
        rho=pr->rho;
        rho_min=pr->rho_min;
        givenRate=pr->givenRate;
        nbSampling=pr->nbSampling;
        ratePartition = pr->ratePartition;
        multiplierRate =  vector<double>();
        if (pr->ratePartition.size()>0) {
            for (int i=0;i<pr->multiplierRate.size();i++){
                multiplierRate.push_back(pr->multiplierRate[i]);
            }
        }
        keepOutgroup = pr->keepOutgroup;
        k = pr->k;
        warningMessage = pr->warningMessage;
        resultMessage = pr->resultMessage;
    }
    Pr()
    {
        inFile = "";
        inDateFile = "";
        partitionFile = "";
        outFile = "";
        treeFile1 = "";
        treeFile2 = "";
        fnOutgroup = "";
        seqLength = 1000;
        nbData = 1;
        rate = "";
        relative = false;
        mrca=0;
        leaves=1;
        estimate_root = "";
        constraint = false;
        variance = 0;
        c = 10;
        rho_min = 1e-10;
        ci = false;
        nbSampling=100;
        rooted=true;
        keepOutgroup=false;
        ratePartition = vector<Part* >();
        internalConstraints = vector<Date* >();
        outlier = vector<int>();
        k = -1;
        verbose = false;
        nbINodes = 0;
        nbBranches = 0;
        rho = 0;
        givenRate = vector<bool>();
        givenRate.push_back(false);
        multiplierRate = vector<double>();
        objective = 0;
        warningMessage = vector<string>();
        resultMessage = vector<string>();
    }
}Pr;
#endif


