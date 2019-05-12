#ifndef DATE_H
#define DATE_H
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

using namespace std;

typedef struct Date
{
	char type;
	double lower;
	double upper;
    double date;
	int id;
    vector<int> mrca;
    Date(){
        type='n';
        id=-1;
    }
    Date(char t,double v1,double v2,int k){
        id=k;
        if (t=='l'){
            type=t;
            lower=v1;
        }
        else if (t=='u'){
            type=t;
            upper=v1;
        }
        else if (t=='b'){
            type=t;
            lower=v1;
            upper=v2;
        }
        else if (t=='p'){
            type=t;
            date=v1;
        }
        else{
            cout<<"unrecognized temporal constraint type"<<endl;
        }
    }
    Date(char t,double v1,double v2,vector<int> mr){
        mrca.clear();
        for (vector<int>::iterator iter=mr.begin();iter!=mr.end();iter++){
            mrca.push_back(*iter);
        }
        id=-1;
        if (t=='l'){
            type=t;
            lower=v1;
        }
        else if (t=='u'){
            type=t;
            upper=v1;
        }
        else if (t=='b'){
            type=t;
            lower=v1;
            upper=v2;
        }
        else if (t=='p'){
            type=t;
            date=v1;
        }
        else{
            cout<<"unrecognized temporal constraint type"<<endl;
        }
    }
} Date;
#endif


