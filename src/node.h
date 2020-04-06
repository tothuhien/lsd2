#ifndef NODE_H
#define NODE_H
#include <vector>
#include "date.h"

using namespace std;

class Node
{
public:
    int P;
    double B;
    string L;
    double V;
    double H;
    double HD;
    vector<int> suc;
    char type;
    double lower;
    double upper;
    double D;
    double minblen;
    int rateGroup;
    int status;
    Node(){
        type='n';
        status=0;
        rateGroup=0;
        minblen = 0;
    }
    bool addConstraint(char t,double v){
        if (type=='n'){
            type=t;
            if (t=='l') lower=v;
            if (t=='u') upper=v;
            if (t=='p') {
                D=v;
                status=8;
            }
            return true;
        }
        else if (type=='l'){
            if (t=='l'){
                lower=max(v,lower);
                return true;
            }
            if (t=='u' && lower<=v){
                if (lower<v) {
                    type='b';
                    upper=v;
                }
                else{
                    type='p';
                    D=lower;
                }
                return true;
            }
            if (t=='p' && lower<=v){
                type='p';
                D=v;
                status=8;
                return true;
            }
        }
        else if (type=='u'){
            if (t=='l' && v<=upper){
                if (v<upper){
                    type='b';
                    upper=v;
                }
                else{
                    type='p';
                    D=upper;
                }
                return true;
            }
            if (t=='u'){
                upper=min(v,upper);
                return true;
            }
            if (t=='p' && v<=upper){
                type='p';
                D=v;
                status=8;
                return true;
            }
        }
        else if (type=='b'){
            if (t=='l' && v<=upper){
                type='b';
                lower=max(v,lower);
                return true;
            }
            if (t=='u' && lower<=v){
                type='b';
                upper=min(v,upper);
                return true;
            }
            if (t=='p' && v<=upper && v>=lower){
                type='p';
                D=v;
                status=8;
                return true;
            }
        }
        else if ((type=='p') && ((t=='l' && v<=D) || (t=='u' && D<=v) || (t=='p' && v==D))){
            return true;
        }
        return false;
    }
    
    bool addConstraint(char t,double l,double u){
        if (t=='b'){
            if (type=='n'){
                type=t;
                lower=l;
                upper=u;
                return true;
            }
            else if (type=='l' && lower<=u){
                if (lower<u){
                    type='b';
                    lower=max(l,lower);
                    upper=u;
                }
                else{
                    type='p';
                    D=u;
                    status=8;
                }
                return true;
            }
            else if (type=='u' && l<=upper){
                if (l<upper){
                    type='b';
                    lower=l;
                    upper=min(u,upper);
                }
                else {
                    type='p';
                    D=l;
                    status=8;
                }
                return true;
            }
            else if (type=='b' && lower<=u && l<=upper){
                if (lower<u && l<upper){
                    lower=max(lower,l);
                    upper=min(upper,u);
                }
                else if (lower==u){
                    type='p';
                    D=lower;
                    status=8;
                }
                else {
                    type='p';
                    D=upper;
                    status=8;
                }
                return true;
            }
            else if ((type=='p') && (l<=D && D<=u)){
                return true;
            }
        }
        return false;
    }
    bool addConstraint(Node* no){
        if (no->type=='b'){
            return addConstraint('b',no->lower,no->upper);
        }
        else if (no->type=='p'){
            return addConstraint('p',no->D);
        }
        else if (no->type=='l'){
            return addConstraint('l',no->lower);
        }
        else if (no->type=='u'){
            return addConstraint('u',no->upper);
        }
        return false;

    }
    bool addConstraint(Date* newdate){
        if (newdate->type=='b'){
            return addConstraint('b',newdate->lower,newdate->upper);
        }
        else if (newdate->type=='p'){
            return addConstraint('p',newdate->date);
        }
        else if (newdate->type=='l'){
            return addConstraint('l',newdate->lower);
        }
        else if (newdate->type=='u'){
            return addConstraint('u',newdate->upper);
        }
        return false;
    }
    void newPConstraint(char t,double v){
        if (t=='p') {
            type=t;
            D=v;
            status=8;
        }
    }
    void removeConstraint(){
        type='n';
        status=0;
    }
    bool addConstraintOfP(Node* no){
        if (no->type=='p'){
            return addConstraint('u',no->D);
        }
        else if (no->type=='b' || no->type=='u'){
            return addConstraint('u',no->upper);
        }
        return true;
    }
    bool addConstraintOfS(Node* no){
        if (no->type=='p'){
            return addConstraint('l',no->D);
        }
        else if (no->type=='b' || no->type=='l'){
            return addConstraint('l',no->lower);
        }
        return true;
    }
};

#endif

