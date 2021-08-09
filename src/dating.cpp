//    LSD - Least Square Dating for etimating substitution rate and divergence dates
//    Copyright (C) <2015> <Thu-Hien To>

//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.

//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.

//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#include "dating.h"


bool without_constraint(Pr* pr,Node** nodes){
    //This function implements LD algorithm using only precise date constrains of the tips (interval date constrains of tips are also ignored)
    //The estimated substitution rate and dates of all internal nodes are returned in variable paramaters.
    //It returns the value of the objective function
    list<int> pos = postorder_polytomy(pr,nodes);
    double *W= new double[pr->nbINodes];//nodes[i]->D=W[i].T[a(i)]+C[i]+X[i]/rho
    double *C = new double[pr->nbINodes];
    double *X =new double[pr->nbINodes];
    for (list<int>::iterator iter=pos.begin();iter!=pos.end();iter++){
        int i =  *iter;
        if (leaf(nodes[i])) {
            W[i]=0;
            C[i]=nodes[i]->D;
            X[i]=0;
        }
        else{
            if (i==0){
                double coefs=0;
                double xtemp=0;
                double ctemp=0;
                for (vector<int>::iterator iter = nodes[i]->suc.begin();iter!=nodes[i]->suc.end();iter++){
                    if (leaf(nodes[*iter])) {
                        coefs+=1/nodes[*iter]->V;
                        ctemp+=nodes[*iter]->D/nodes[*iter]->V;
                        xtemp-=nodes[*iter]->B/nodes[*iter]->V;
                    }
                    else if (*iter<pr->nbINodes){
                        coefs+=(1-W[*iter])/nodes[*iter]->V;
                        ctemp+=C[*iter]/nodes[*iter]->V;
                        xtemp+=(X[*iter]-nodes[*iter]->B)/nodes[*iter]->V;
                    }
                }
                if (coefs==0) {
                    return false;
                }
                C[i]=ctemp/coefs;
                X[i]=xtemp/coefs;
            }
            else{
                double coefs=1/nodes[i]->V;
                double wtemp=1/nodes[i]->V;
                double ctemp=0;
                double xtemp=nodes[i]->B/nodes[i]->V;
                for (vector<int>::iterator iter = nodes[i]->suc.begin();iter!=nodes[i]->suc.end();iter++){
                    if (leaf(nodes[*iter])) {
                        coefs+=1/nodes[*iter]->V;
                        ctemp+=nodes[*iter]->D/nodes[*iter]->V;
                        xtemp-=nodes[*iter]->B/nodes[*iter]->V;
                    }
                    else if (*iter<pr->nbINodes){
                        coefs+=(1-W[*iter])/nodes[*iter]->V;
                        ctemp+=C[*iter]/nodes[*iter]->V;
                        xtemp+=(X[*iter]-nodes[*iter]->B)/nodes[*iter]->V;
                    }
                }
                W[i]=wtemp/coefs;
                C[i]=ctemp/coefs;
                X[i]=xtemp/coefs;
            }
        }
    }
    vector<int> pre = preorder_polytomy(pr,nodes);
    for (vector<int>::iterator iter = pre.begin();iter!=pre.end();iter++){
        int i = *iter;
        if (i!=0 && !leaf(nodes[i])){
            C[i]=W[i]*C[nodes[i]->P]+C[i];
            X[i]=W[i]*X[nodes[i]->P]+X[i];
        }
    }
    if (!pr->givenRate[0]){
        double *F = new double[pr->nbBranches+1];
        double *G = new double[pr->nbBranches+1];
        for (int i=1;i<=pr->nbBranches;i++){
            if (i<pr->nbINodes || leaf(nodes[i])){
                if (leaf(nodes[nodes[i]->P]) && leaf(nodes[i])){
                    F[i]=nodes[nodes[i]->P]->D-nodes[i]->D;
                    G[i]=0;
                }
                else if (leaf(nodes[nodes[i]->P])){
                    F[i]=nodes[nodes[i]->P]->D-C[i];
                    G[i]=-X[i];
                }
                else if (leaf(nodes[i])){
                    F[i]=C[nodes[i]->P]-nodes[i]->D;
                    G[i]=X[nodes[i]->P];
                }
                else {
                    F[i]=C[nodes[i]->P]-C[i];
                    G[i]=X[nodes[i]->P]-X[i];
                }
            }
        }
        double a = 0;
        double b = 0;
        for (int i=1;i<=pr->nbBranches;i++){
            if ( i<pr->nbINodes || leaf(nodes[i])){
                a += F[i]*F[i]/nodes[i]->V;
                b += 2*(nodes[i]->B+G[i])*F[i]/nodes[i]->V;
            }
        }
        if (abs(a) < 1e-10) return false;
        pr->rho = -b/(2*a);
        if (pr->rho<pr->rho_min) pr->rho=pr->rho_min;
        delete[] F;
        delete[] G;
    }
    for (int i =0;i<pr->nbINodes;i++) {
        if (!leaf(nodes[i])) {
            nodes[i]->D = C[i]+X[i]/pr->rho;
        }
        
    }
    for (int i=pr->nbINodes; i<=pr->nbBranches; i++) {
        if (!leaf(nodes[i])) {
            nodes[i]->D=nodes[nodes[i]->P]->D+nodes[i]->B/pr->rho;
        }
    }
    computeObjective(pr,nodes);
    delete[] W;
    delete[] C;
    delete[] X;
    return true;
}

bool conditions(list<double>& ldLagrange,Pr* pr,Node** nodes){
    //the stop condition of active set method
    for (list<double>::iterator iter = ldLagrange.begin();iter!=ldLagrange.end();iter++){
        if (*iter < 0)  {
            return false;
        }
    }
    for (int i=0;i<=pr->nbBranches;i++){
        if (((nodes[i]->type=='l' || nodes[i]->type=='b') && nodes[i]->D<nodes[i]->lower) || ((nodes[i]->type=='u' || nodes[i]->type=='b') && nodes[i]->D>nodes[i]->upper)){
            return false;
        }
        
    }
    return true;
}


list<double> computeLambda(list<int> active_set,Pr* pr,Node** nodes){
    //this method computes the optimized solution of the Lagrange function.
    //it returns the list of Lagrange multipliers
    int *as=new int[pr->nbBranches+1];
    for (int i=0;i<=pr->nbBranches;i++) as[i]=-1;
    int count = 0;
    for (list<int>::iterator iter = active_set.begin();iter!=active_set.end();iter++){
        int a = *iter;
        as[-a]=count;
        count++;
    }
    double* lambda = new double[count];
    list<double> ld;
    for (list<int>::iterator it = active_set.begin();it!=active_set.end();it++){
        int i=-(*it);
        lambda[as[i]]=0;
        if (lower(nodes[i])){
            if (i!=0) lambda[as[i]]-=2*pr->rho*(nodes[i]->B+pr->rho*nodes[nodes[i]->P]->D-pr->rho*nodes[i]->D)/nodes[i]->V;
            if (i<pr->nbINodes) {
                for (vector<int>::iterator iter=nodes[i]->suc.begin(); iter!=nodes[i]->suc.end(); iter++) {
                    int s=*iter;
                    lambda[as[i]] += 2*pr->rho*(nodes[s]->B+pr->rho*nodes[i]->D-pr->rho*nodes[s]->D)/nodes[s]->V;
                }
            }
            if (abs(lambda[as[i]])<1e-10) lambda[as[i]]=0;
            ld.push_back(lambda[as[i]]);
        }
        else if (upper(nodes[i])){
            if (i!=0) {
                lambda[as[i]]=2*pr->rho*(nodes[i]->B+pr->rho*nodes[nodes[i]->P]->D-pr->rho*nodes[i]->D)/nodes[i]->V;
            }
            if (i<pr->nbINodes) {
                for (vector<int>::iterator iter=nodes[i]->suc.begin(); iter!=nodes[i]->suc.end(); iter++) {
                    int s=*iter;
                    lambda[as[i]] -= 2*pr->rho*(nodes[s]->B+pr->rho*nodes[i]->D-pr->rho*nodes[s]->D)/nodes[s]->V;
                }
            }
            if (abs(lambda[as[i]])<1e-10) lambda[as[i]]=0;
            ld.push_back(lambda[as[i]]);
        }
    }
    /*  double sr=0;
     for (int i=0; i<=pr->nbBranches; i++) {
     if (nodes[i]->type!='p') {
     double s=0;
     if (i>0) {
     s=-2*pr->rho*(nodes[i]->B-pr->rho*(nodes[i]->D-nodes[nodes[i]->P]->D))/nodes[i]->V;
     }
     if (i<pr->nbINodes) {
     for (vector<int>::iterator iter=nodes[i]->suc.begin(); iter!=nodes[i]->suc.end(); iter++) {
     int s1=*iter;
     s += 2*pr->rho*(nodes[s1]->B+pr->rho*nodes[i]->D-pr->rho*nodes[s1]->D)/nodes[s1]->V;
     }
     }
     if (lower(nodes[i])){
     s-=lambda[as[i]];
     }
     if (upper(nodes[i])){
     s+=lambda[as[i]];
     }
     if (abs(s)>1e-10) {
     cout<<"TEST PROBLEM "<<i<<" "<<s<<" "<<nodes[i]->P<<endl;
     }
     
     }
     if (i>0) sr+=(nodes[i]->D-nodes[nodes[i]->P]->D)*(nodes[i]->B-pr->rho*(nodes[i]->D-nodes[nodes[i]->P]->D))/nodes[i]->V;
     }
     if (abs(sr)>1e-6) {
     cout<<"TEST PROBLEM rho "<<sr<<endl;
     }*/
    delete[] lambda;
    delete[] as;
    return ld;
}

bool remove_ne_lambda(list<double> & lambda,list<int> & active_set,int& as){
    //remove the active set corresponding to the most negative lagrange multiplier
    double ne_lambda=0;
    double min = 0;
    list<int>::iterator iteras=active_set.begin();
    for (list<double>::iterator iter = lambda.begin();iter!=lambda.end();iter++){
        ne_lambda = *iter;
        if (ne_lambda<min){
            min = ne_lambda;
            as=*iteras;
        }
        iteras++;
    }
    if (min<0) {
        return true;
    }
    else return false;
}


bool without_constraint_active_set(Pr* pr,Node** nodes,int whichStartingPoint){
    //this methods implements the LD algorithm (active set method)
    initialize_status(pr,nodes);
    list<int> active_set;
    if (whichStartingPoint==-1) starting_pointLower(pr,nodes,active_set);
    if (whichStartingPoint==1) starting_pointUpper(pr,nodes,active_set);
    //starting_point(pr,nodes,active_set);
    double* D_old = new double[pr->nbBranches+1];
    for (int i=0; i<=pr->nbBranches; i++) {
        D_old[i]=nodes[i]->D;
    }
    bool val = without_constraint(pr,nodes);
    /*if (!val) {
        cerr<<"Error: There's not enough signal in the input temporal constraints to have unique solution."<<endl;
        exit(EXIT_FAILURE);
    }*/
    list<double> lambda = computeLambda(active_set,pr,nodes);
    int nb_iter=0;
    double alpha;
    double* dir = new double[pr->nbBranches+1];
    while (val && !conditions(lambda,pr,nodes)){
        for (int i=0;i<=pr->nbBranches;i++)  {dir[i]=nodes[i]->D-D_old[i];}
        alpha=1;
        int as=0;
        double a;
        for (int i=0;i<=pr->nbBranches;i++){
            if (nodes[i]->status==0 && (nodes[i]->type=='l' || nodes[i]->type=='u' || nodes[i]->type=='b')){
                if (dir[i]<0 && (nodes[i]->type=='l' || nodes[i]->type=='b') ){
                    a = (nodes[i]->lower-D_old[i])/dir[i];
                    if (a<alpha){
                        alpha = a;
                        as = i+pr->nbBranches+1;
                    }
                }
                if (dir[i]>0 && (nodes[i]->type=='u' || nodes[i]->type=='b')){
                    a = (nodes[i]->upper-D_old[i])/dir[i];
                    if (a<alpha){
                        alpha = a;
                        as = i+2*pr->nbBranches+2;
                    }
                }
            }
        }
        int asrm;
        if (remove_ne_lambda(lambda,active_set,asrm)){
            active_set.remove(asrm);
            desactive(nodes[-asrm]);
        }
        for (int i=0;i<=pr->nbBranches;i++) D_old[i]=D_old[i]+alpha*dir[i];
        if (as!=0) {
            if (as>2*pr->nbBranches+1) {
                active_set.push_back(-(as-2*pr->nbBranches-2));
                activeUpper(nodes[as-2*pr->nbBranches-2]);
                nodes[as-2*pr->nbBranches-2]->D=nodes[as-2*pr->nbBranches-2]->upper;
            }
            else {
                active_set.push_back(-(as-pr->nbBranches-1));
                activeLower(nodes[as-pr->nbBranches-1]);
                nodes[as-pr->nbBranches-1]->D=nodes[as-pr->nbBranches-1]->lower;
            }
        }
        val = without_constraint(pr,nodes);
        lambda=computeLambda(active_set,pr,nodes);
        nb_iter++;
    }
    /*if (nb_iter>1000){
        for (int i=0;i<=pr->nbBranches;i++) nodes[i]->D=D_old[i];
    }*/
    computeObjective(pr,nodes);
    delete[] D_old;
    delete[] dir;
    return val;
}



/////////////////////////////////////////////////////////////////////////////////////////////////
//////////////This part implements QPD algorithm (active set methods)////////////////////////////

bool conditionsQP(list<double>& ldLagrange,Pr* pr,Node** nodes){
    //the stop condition of active set method
    for (list<double>::iterator iter = ldLagrange.begin();iter!=ldLagrange.end();iter++){
        if (*iter < 0)  {
            return false;
        }
    }
    for (int i=0;i<=pr->nbBranches;i++){
        if (i>0 && (nodes[i]->D - nodes[nodes[i]->P]->D - nodes[i]->minblen) < -1e-10) {
            return false;
        }
        if (((nodes[i]->type=='l' || nodes[i]->type=='b') && nodes[i]->D<nodes[i]->lower) || ((nodes[i]->type=='u' || nodes[i]->type=='b') && nodes[i]->D>nodes[i]->upper)){
            return false;
        }
        
    }
    return true;
}


bool starting_pointQP(Pr* pr,Node** nodes,list<int> &active_set,int whichStartingPoint){
    //compute a starting feasible point, which is the solution from without constraint and then collapse all the branches that violate the constraints.
    bool val = without_constraint_active_set(pr,nodes,whichStartingPoint);
    if (!val) return val;
    double* lowerX = new double[pr->nbBranches+1];
    bool* bl = new bool[pr->nbBranches+1];
    for (int i=0;i<=pr->nbBranches;i++){
        if (nodes[i]->type=='l' || nodes[i]->type=='b'){
            bl[i]=true;
            lowerX[i]=nodes[i]->lower;
        }
        else if (nodes[i]->type=='p'){
            bl[i]=true;
            lowerX[i]=nodes[i]->D;
        }
        else bl[i]=false;
    }
    vector<int> pre = preorder_polytomy(pr,nodes); //top-bottom order
    for (vector<int>::iterator iter = pre.begin();iter!=pre.end();iter++){
        int i=*iter;
        if (bl[i]){
            for (vector<int>::iterator iter=nodes[i]->suc.begin(); iter!=nodes[i]->suc.end(); iter++) {
                int s=*iter;
                if (!bl[s] || (nodes[s]->type == 'l' && (nodes[s]->lower - lowerX[i] < nodes[s]->minblen))){
                    lowerX[s] = lowerX[i] + nodes[s]->minblen;
                    nodes[s]->lower = lowerX[i] + nodes[s]->minblen;
                    bl[s]=true;
                }
                else if (nodes[s]->type == 'u' && nodes[s]->upper - lowerX[i] < nodes[s]->minblen){
                    return false;
                }
                else if (nodes[s]->type == 'b'){
                    if (nodes[s]->upper - lowerX[i] >= nodes[s]->minblen){
                        if (nodes[s]->lower - lowerX[i] < nodes[s]->minblen){
                            lowerX[s] = lowerX[i] + nodes[s]->minblen;
                            nodes[s]->lower = lowerX[i] + nodes[s]->minblen;
                            bl[s]=true;
                        }
                    }
                    else {
                        return false;
                    }
                }
                else if ((nodes[s]->type == 'p') && nodes[s]->D - lowerX[i] < nodes[s]->minblen){
                    return false;
                }
                if (nodes[s]->D < lowerX[s]){
                    nodes[s]->D = lowerX[s];
                }
            }
        }
    }
    for (int i=pr->nbINodes; i<=pr->nbBranches; i++) {
        if (lower(nodes[i]) || upper(nodes[i])) {
            active_set.push_back(-i);
        }
        else if (bl[i] && nodes[i]->type != 'p' && nodes[i]->D < lowerX[i]){
            nodes[i]->D = lowerX[i];
            if ((nodes[i]->type=='l' || nodes[i]->type=='b')) {
                activeLower(nodes[i]);
                active_set.push_back(-i);
            }
        }
    }
    list<int> pos = postorder_polytomy(pr,nodes);
    for (list<int>::iterator iter = pos.begin();iter!=pos.end();iter++){
        int i = *iter;
        if (lower(nodes[i]) || upper(nodes[i])) {
            active_set.push_back(-i);
        }
        else if (bl[i] && nodes[i]->type != 'p' && nodes[i]->D < lowerX[i]){
            nodes[i]->D = lowerX[i];
            if ((nodes[i]->type=='l' || nodes[i]->type == 'b')) {
                activeLower(nodes[i]);
                active_set.push_back(-i);
            }
        }
        bool conflictU = ((nodes[i]->type == 'u' || nodes[i]->type == 'b') && nodes[i]->D > nodes[i]->upper);
        bool conflictTC = false;
        int minI = i;
        double minS = nodes[i]->D;
        for (vector<int>::iterator iter=nodes[i]->suc.begin(); iter!=nodes[i]->suc.end(); iter++) {
            int s=*iter;
            if (nodes[s]->D - nodes[i]->D < nodes[s]->minblen){
                conflictTC = true;
            }
            if (nodes[s]->D - nodes[s]->minblen < minS) {
                minS = nodes[s]->D - nodes[s]->minblen;
                minI = s;
            }
        }
        if (conflictU || conflictTC){
            if (((nodes[i]->type == 'u' || nodes[i]->type == 'b') && minS < nodes[i]->upper) || (nodes[i]->type != 'u' && nodes[i]->type != 'b')){
                if (lower(nodes[i]) || upper(nodes[i])){
                    desactiveLimit(nodes[i]);
                    active_set.remove(-i);
                }
                nodes[i]->D = minS;
                if (!tc(nodes[minI])){
                    activeTC(nodes[minI]);
                    active_set.push_back(minI);
                }
            }
            else{
                if (lower(nodes[i]) || upper(nodes[i])){
                    desactiveLimit(nodes[i]);
                    active_set.remove(-i);
                }
                nodes[i]->D = nodes[i]->upper;
                activeUpper(nodes[i]);
                active_set.push_back(-i);
            }
        }
    }
    delete[] bl;
    delete[] lowerX;
    computeObjective(pr,nodes);
    return true;
}


bool with_constraint(Pr* pr,Node** &nodes,list<int> active_set,list<double>& ld){
    //this method computes the optimized solution of the Lagrange function.
    //it returns the list of Lagrange multipliers
    list<int>* Suc = new list<int>[pr->nbINodes];
    int* Pre=new int[pr->nbBranches+1];
    for (int i=0;i<=pr->nbBranches;i++){
        Pre[i]=-1;
        if (leaf(nodes[i])) {
            activeMarkLeaf(nodes[i]);
        }
        else desactiveMarkLeaf(nodes[i]);
    }
    list<int> ls;
    for (list<int>::iterator iter = active_set.begin();iter!=active_set.end();iter++){
        int i = *iter;
        if (i>0) {//time constraint
            if (leaf(nodes[i])) {
                ls.push_back(i);
            }
        }
        else if (!tc(nodes[-i]) && leaf(nodes[-i])) ls.push_back(-i);
    }
    stack<int>* feuilles = computeFeuilles_polytomy(ls,pr,nodes);
    list<int> top;
    for (int i=0;i<pr->nbINodes;i++){
        bool bl=false;
        for (vector<int>::iterator iter=nodes[i]->suc.begin(); iter!=nodes[i]->suc.end(); iter++) {
            if (tc(nodes[*iter])) {
                bl=true;
            }
        }
        if (!tc(nodes[i]) && (!markLeaf(nodes[i]) || nodes[i]->type=='p')  && bl){
            top.push_back(i);
        }
    }
    list<int>* internal = new list<int>[top.size()];
    double *add = new double[pr->nbBranches+1];
    for (int i=0;i<=pr->nbBranches;i++) add[i] = 0;
    reduceTree_polytomy(pr,nodes,Pre,Suc,add,internal);
    list<int> pos = postorder_polytomy(pr,nodes);
    vector<int> pre = preorder_polytomy(pr,nodes);
    double *W= new double[pr->nbINodes];//nodes[i]->D=W[i].T[a(i)]+C[i]+X[i]/rho
    double *C = new double[pr->nbINodes];
    double *X =new double[pr->nbINodes];
    for (int i=0;i<pr->nbINodes;i++){
        if (markLeaf(nodes[i])) C[i]=nodes[i]->D;
        else C[i] = 0;
        W[i]=0;X[i]=0;
    }
    for (list<int>::iterator it = pos.begin();it!=pos.end();it++){
        int i = *it;
        if (nodes[i]->status==0){
            int p = Pre[i];
            if (p==-1){
                double coefs=0;
                double xtemp=0;
                double ctemp=0;
                for (list<int>::iterator iter = Suc[i].begin();iter!=Suc[i].end();iter++){
                    if (markLeaf(nodes[*iter])) {
                        coefs+=1/nodes[*iter]->V;
                        ctemp+=(nodes[*iter]->D-add[*iter])/nodes[*iter]->V;
                        xtemp-=nodes[*iter]->B/nodes[*iter]->V;
                    }
                    else{
                        coefs+=(1-W[*iter])/nodes[*iter]->V;
                        ctemp+=(C[*iter]-add[*iter])/nodes[*iter]->V;
                        xtemp+=(X[*iter]-nodes[*iter]->B)/nodes[*iter]->V;
                    }
                }
                if (coefs==0) {
                    return false;
                }
                C[i]=ctemp/coefs;
                X[i]=xtemp/coefs;
            }
            else{
                double coefs=1/nodes[i]->V;
                double wtemp=1/nodes[i]->V;
                double ctemp=add[i]/nodes[i]->V;
                double xtemp=nodes[i]->B/nodes[i]->V;
                for (list<int>::iterator iter = Suc[i].begin();iter!=Suc[i].end();iter++){
                    if (markLeaf(nodes[*iter])) {
                        coefs+=1/nodes[*iter]->V;
                        ctemp+=(nodes[*iter]->D-add[*iter])/nodes[*iter]->V;
                        xtemp-=nodes[*iter]->B/nodes[*iter]->V;
                    }
                    else{
                        coefs+=(1-W[*iter])/nodes[*iter]->V;
                        ctemp+=(C[*iter]-add[*iter])/nodes[*iter]->V;
                        xtemp+=(X[*iter]-nodes[*iter]->B)/nodes[*iter]->V;
                    }
                }
                W[i]=wtemp/coefs;
                C[i]=ctemp/coefs;
                X[i]=xtemp/coefs;
            }
        }
    }
    for (vector<int>::iterator iter = pre.begin();iter!=pre.end();iter++){
        int i = *iter;
        if (Pre[i]!=-1 && nodes[i]->status==0){
            C[i]=W[i]*C[Pre[i]]+C[i];
            X[i]=W[i]*X[Pre[i]]+X[i];
        }
        for (vector<int>::iterator iter=nodes[i]->suc.begin(); iter!=nodes[i]->suc.end(); iter++) {
            int s=*iter;
            if (tc(nodes[s]) && !markLeaf(nodes[s]) && s<pr->nbINodes) {
                C[s]=C[i] + nodes[s]->minblen;
                X[s]=X[i];
            }
        }
    }
    if (!pr->givenRate[0]){//if the rate is not given
        double *F=new double[pr->nbBranches+1];
        double *G=new double[pr->nbBranches+1];//nodes[nodes[i]->P]->D-nodes[i]->D=F[i]+G[i]/rho;
        for (int i=0;i<=pr->nbBranches;i++){
            F[i]=0;G[i]=0;
        }
        for (int i=1;i<=pr->nbBranches;i++){
            if (i<pr->nbINodes || markLeaf(nodes[i])){
                if (markLeaf(nodes[nodes[i]->P]) && markLeaf(nodes[i])){
                    F[i]=nodes[nodes[i]->P]->D-nodes[i]->D;
                    G[i]=0;
                }
                else if (markLeaf(nodes[nodes[i]->P])){
                    F[i]=nodes[nodes[i]->P]->D-C[i];
                    G[i]=-X[i];
                }
                else if (markLeaf(nodes[i])){
                    F[i]=C[nodes[i]->P]-nodes[i]->D;
                    G[i]=X[nodes[i]->P];
                }
                else {
                    F[i]=C[nodes[i]->P]-C[i];
                    G[i]=X[nodes[i]->P]-X[i];
                }
            }
        }
        double a = 0;
        double b = 0;//phi = a*rho^2+b*rho+c
        for (int i=1;i<=pr->nbBranches;i++){
            if (i<pr->nbINodes || markLeaf(nodes[i])){
                a += F[i]*F[i]/nodes[i]->V;
                b += 2*(nodes[i]->B+G[i])*F[i]/nodes[i]->V;
            }
        }
        pr->rho = -b/(2*a);
        if (pr->rho < pr->rho_min) {
            pr->rho=pr->rho_min;
        }
        delete[] F;
        delete[] G;
    }
    for (vector<int>::iterator iter = pre.begin();iter!=pre.end();iter++){
        int i = *iter;
        nodes[i]->D = C[i]+X[i]/pr->rho;
    }
    for (int i=pr->nbINodes; i<=pr->nbBranches; i++) {
        if (!markLeaf(nodes[i])) {
            if (tc(nodes[i])) {
                nodes[i]->D = nodes[nodes[i]->P]->D + nodes[i]->minblen;
            }
            else nodes[i]->D = nodes[nodes[i]->P]->D + nodes[i]->B/pr->rho;
        }
    }
    computeObjective(pr,nodes);
    delete[] W;
    delete[] C;
    delete[] X;
    delete[] Pre;
    delete[] Suc;
    delete[] add;
    int *as=new int[2*pr->nbBranches+1];
    for (int i=0;i<=2*pr->nbBranches;i++) as[i]=-1;
    int count = 0;
    for (list<int>::iterator iter = active_set.begin();iter!=active_set.end();iter++){
        int a = *iter;
        if (a>0) {
            as[a]=count;
        }
        else{
            as[-a+pr->nbBranches+1]=count;
        }
        count++;
    }
    bool *bl=new bool[pr->nbBranches+1];
    for (int i=0;i<=pr->nbBranches;i++) bl[i]=false;
    double* lambda = new double[count];
    for (int y=0;y<ls.size();y++){
        int i=0;
        while (!feuilles[y].empty()){
            i=feuilles[y].top();
            feuilles[y].pop();
            if (tc(nodes[i])){//compute lambda of active tc
                int p=nodes[i]->P;
                lambda[as[i]] = -2*pr->rho * (nodes[i]->B - pr->rho*nodes[i]->D + pr->rho*nodes[p]->D) / nodes[i]->V;
                if (i>=pr->nbINodes && !leaf(nodes[i])) bl[i]=true;
                if (i<pr->nbINodes && !leaf(nodes[i])){
                    if (!limit(nodes[i])){
                        bl[i]=true;
                        for (vector<int>::iterator iter=nodes[i]->suc.begin(); iter!=nodes[i]->suc.end(); iter++) {
                            int s=*iter;
                            if (!tc(nodes[s]) || bl[s]) {
                                lambda[as[i]] += 2*pr->rho * (nodes[s]->B - pr->rho*nodes[s]->D + pr->rho*nodes[i]->D) / nodes[s]->V;
                                if (tc(nodes[s])) {
                                    lambda[as[i]] += lambda[as[s]];
                                }
                            }
                            else bl[i]=false;
                        }
                    }
                }
                if (!bl[i]){//top-down order
                    lambda[as[i]]=0;
                    for (vector<int>::iterator iter=nodes[p]->suc.begin(); iter!=nodes[p]->suc.end(); iter++) {
                        int j=*iter;
                        lambda[as[i]] -= 2*pr->rho * (nodes[j]->B - pr->rho*nodes[j]->D + pr->rho*nodes[p]->D) / nodes[j]->V;
                        if (j!=i && tc(nodes[j])) {
                            lambda[as[i]] -= lambda[as[j]];
                        }
                    }
                    if (nodes[p]->P!=-1){
                        lambda[as[i]] += 2*pr->rho * (nodes[p]->B - pr->rho*nodes[p]->D + pr->rho*nodes[nodes[p]->P]->D) / nodes[p]->V;
                    }
                    if (tc(nodes[p])) {
                        lambda[as[i]] += lambda[as[p]];
                    }
                    bl[i]=true;
                }
            }
        }
        if (limit(nodes[i])){//compute lambda of active limit constraints
            double lambdaPC=0;
            for (vector<int>::iterator iter=nodes[i]->suc.begin(); iter!=nodes[i]->suc.end(); iter++) {
                int s=*iter;
                lambdaPC+=2*pr->rho*(nodes[s]->B-pr->rho*nodes[s]->D+pr->rho*nodes[i]->D)/nodes[s]->V;
                if (tc(nodes[s])){
                    lambdaPC+=lambda[as[s]];
                }
            }
            if (i!=0) {
                lambdaPC-=2*pr->rho*(nodes[i]->B-pr->rho*nodes[i]->D+pr->rho*nodes[nodes[i]->P]->D)/nodes[i]->V;
            }
            if (tc(nodes[i])){
                lambdaPC-=lambda[as[i]];
            }
            if (lower(nodes[i])){
                lambda[as[i+pr->nbBranches+1]]=lambdaPC;
            }
            else if (upper(nodes[i])){
                lambda[as[i+pr->nbBranches+1]]=-lambdaPC;
            }
        }
    }
    for (int y=0;y<top.size();y++){
        for (list<int>::iterator iter = internal[y].begin();iter!=internal[y].end();iter++){
            int i= *iter;
            int p = nodes[i]->P;
            lambda[as[i]] = -2*pr->rho * (nodes[i]->B - pr->rho*nodes[i]->D + pr->rho*nodes[p]->D) / nodes[i]->V;
            if (i<pr->nbINodes) {
                for (vector<int>::iterator iter=nodes[i]->suc.begin(); iter!=nodes[i]->suc.end(); iter++) {
                    int s=*iter;
                    lambda[as[i]] += 2*pr->rho*(nodes[s]->B - pr->rho*nodes[s]->D + pr->rho*nodes[i]->D) / nodes[s]->V;
                    if (tc(nodes[s])) {
                        lambda[as[i]] += lambda[as[s]];
                    }
                }
                bl[i]=true;
            }
        }
    }
/*
    double sr=0;
    for (int i=0; i<=pr->nbBranches; i++) {
        if (nodes[i]->type!='p') {
            double s=0;
            if (i>0) {
                s = -2*pr->rho * (nodes[i]->B - pr->rho*(nodes[i]->D - nodes[nodes[i]->P]->D)) / nodes[i]->V;
            }
            if (tc(nodes[i])) {
                s -= lambda[as[i]];
            }
            if (lower(nodes[i])){
                s -= lambda[as[i+pr->nbBranches+1]];
            }
            if (upper(nodes[i])){
                s += lambda[as[i+pr->nbBranches+1]];
            }
            if (i<pr->nbINodes) {
                for (vector<int>::iterator iter=nodes[i]->suc.begin(); iter!=nodes[i]->suc.end(); iter++) {
                    int s1=*iter;
                    s += 2*pr->rho * (nodes[s1]->B - pr->rho*(nodes[s1]->D - nodes[i]->D)) / nodes[s1]->V;
                    if (tc(nodes[s1])) {
                        s+=lambda[as[s1]];
                    }
                }
            }
            if (abs(s)>1e-10) {
                cout<<"TEST  PROBLEM "<<i<<" "<<s<<" "<<nodes[i]->P<<endl;
            }
        }
        if (i>0) sr += (nodes[i]->D-nodes[nodes[i]->P]->D) * (nodes[i]->B-pr->rho*(nodes[i]->D-nodes[nodes[i]->P]->D)) /nodes[i]->V;
    }
    if (abs(sr)>1e-6) {
        cout<<"TEST  PROBLEM rho "<<sr<<endl;
    }*/
    for (int i=0;i<count;i++){
        if (abs(lambda[i])<(1e-10))
            lambda[i]=0;
        ld.push_back(lambda[i]);
    }
    delete[] as;
    delete[] feuilles;
    delete[] bl;
    delete[] lambda;
    delete[] internal;
    return true;
}

bool with_constraint_active_set(Pr* pr,Node** &nodes,int whichStartingPoint){
    //this methods implements the QPD algorithm (active set method)
    list<int> active_set;
    bool consistent=starting_pointQP(pr,nodes,active_set,whichStartingPoint);
    if (consistent){
        list<double> lambda;
        double* dir = new double[pr->nbBranches+1];
        double* D_old = new double[pr->nbBranches+1];
        
        for (int i=0; i<=pr->nbBranches; i++) {
            D_old[i]=nodes[i]->D;
        }
        bool val = with_constraint(pr,nodes,active_set,lambda);
        double alpha;
        int nb_iter = 0;
        while (val && !conditionsQP(lambda,pr,nodes)){
            //for (list<double>::iterator iter = lambda.begin(); iter != lambda.end(); iter++) cout<<*iter<<" "; cout<<endl;
            for (int i=0;i<=pr->nbBranches;i++)  {
                dir[i] = nodes[i]->D - D_old[i];
            }
            alpha=1;
            int as=2*pr->nbBranches+2;
            double a;
            for (int i=0; i<=pr->nbBranches; i++) {
                int p = nodes[i]->P;
                if (p!=-1){
                    if (dir[p] > dir[i] && !tc(nodes[i])){
                        a = (D_old[i]-D_old[p]-nodes[i]->minblen)/(dir[p]-dir[i]);
                        if (a<alpha){
                            alpha = a;
                            as = i;
                        }
                    }
                }
                if (dir[i]<0 && (nodes[i]->type=='l' || nodes[i]->type=='b') && !lower(nodes[i])){
                    a = (nodes[i]->lower-D_old[i])/dir[i];
                    if (a<alpha){
                        alpha = a;
                        as = -i;
                    }
                }
                if (dir[i]>0 && (nodes[i]->type=='u' || nodes[i]->type=='b') && !upper(nodes[i])){
                    a = (nodes[i]->upper-D_old[i])/dir[i];
                    if (a<alpha){
                        alpha = a;
                        as = i+pr->nbBranches+1;
                    }
                }
            }
            int asrm;
            if (remove_ne_lambda(lambda,active_set,asrm)){
                active_set.remove(asrm);//cout<<"remove "<<asrm<<endl;
                if (asrm>0) {
                    desactiveTC(nodes[asrm]);
                }
                else {
                    desactiveLimit(nodes[-asrm]);
                }
            }
            for (int i=0;i<=pr->nbBranches;i++) D_old[i]=D_old[i]+alpha*dir[i];
            if (as<2*pr->nbBranches+2) {
                if (as>pr->nbBranches) {
                    active_set.push_back(-(as-pr->nbBranches-1));//cout<<"add "<<-(as-pr->nbBranches-1)<<endl;
                    nodes[(as-pr->nbBranches-1)]->D=nodes[(as-pr->nbBranches-1)]->upper;
                    activeUpper(nodes[(as-pr->nbBranches-1)]);
                }
                else if (as>0) {
                    active_set.push_back(as);//cout<<"add "<<as<<endl;
                    activeTC(nodes[as]);
                }
                else{
                    active_set.push_back(as);//cout<<"add "<<as<<endl;
                    nodes[-as]->D=nodes[-as]->lower;
                    activeLower(nodes[-as]);
                }
            }
            lambda.clear();
            val=with_constraint(pr,nodes,active_set,lambda);
            nb_iter++;
        }
        /*if (nb_iter>1000){
            for (int i=0;i<=pr->nbBranches;i++) nodes[i]->D=D_old[i];
        }*/
        computeObjective(pr,nodes);
        delete[] D_old;
        delete[] dir;
        return val;
    }
    else return false;
}


// Multi rates


void calculateMultiplier(Pr* pr,Node** nodes){
    int nbG = pr->ratePartition.size()+1;
    double* A =  new double[nbG];
    double* B =  new double[nbG];
    double* C =  new double[nbG];
    for (int i=1;i<nbG;i++){
        A[i] = 0;
        B[i] = 0;
        C[i] = 0;
    }
    for (int i=1; i<=pr->nbBranches; i++) {
        int g = nodes[i]->rateGroup;
        A[g] += pr->rho*pr->rho*(nodes[i]->D-nodes[nodes[i]->P]->D)*(nodes[i]->D-nodes[nodes[i]->P]->D)/nodes[i]->V;
        B[g] += -2*nodes[i]->B*pr->rho*(nodes[i]->D-nodes[nodes[i]->P]->D)/nodes[i]->V;
        C[g] += nodes[i]->B*nodes[i]->B/nodes[i]->V;
    }
    for (int i=1; i<nbG; i++) {
        if (!pr->givenRate[i]){
            pr->multiplierRate[i] = -B[i]/2/A[i];
            if (pr->multiplierRate[i]*pr->rho < pr->rho_min){
                pr->multiplierRate[i] = pr->rho_min/pr->rho;
            }
            if (A[i]==0) {
                pr->multiplierRate[i] = -1;
            }
        }
    }
}

bool without_constraint_active_set(bool all,Pr* pr,Node** nodes){
    bool val = false;
    if (pr->haveUnique) {
        val = without_constraint_active_set(pr,nodes,0);
    }
    if (!val){
        bool valL = false;
        bool valU = false;
        if (pr->haveLower){
            valL = without_constraint_active_set(pr,nodes,-1);
            if (valL){
                pr->rhoLower = pr->rho;
                val = true;
            }
        }
        if (!pr->haveLower || (all && pr->haveUpper)){
            valU = without_constraint_active_set(pr,nodes,1);
            if (valU){
                pr->rhoUpper = pr->rho;
                val = true;
            }
        }
    }
    return val;
}

bool without_constraint_multirates(Pr* pr,Node** nodes,bool reassign){
    double* B = new double[pr->nbBranches+1];
    double* V = new double[pr->nbBranches+1];
    for (int i=1; i<=pr->nbBranches; i++) {
        B[i] = nodes[i]->B;
        V[i] = nodes[i]->V;
    }
    if (pr->ratePartition.size()>0){
        if (reassign) assignRateGroupToTree(pr,nodes);
        for (int r=1; r<=pr->nbBranches; r++) {
            double m = pr->multiplierRate[nodes[r]->rateGroup];
            nodes[r]->B = B[r]/m;
            nodes[r]->V = V[r]/m/m;
        }
    }
    bool val = without_constraint_active_set(true,pr,nodes);
    if (!val) return false;
    if (pr->ratePartition.size()>0) {
        if (pr->verbose){
            printf("ROUND 0 , objective function %.15e , rate %.15f , other_rates ",pr->objective,pr->rho);
            for (int r=1; r<=pr->ratePartition.size(); r++) {
                printf(" %.15f",pr->rho*pr->multiplierRate[r]);
            }
            cout<<endl;
        }
        double old_phi = 0;
        double old_rho = 0;
        double* old_multi = new double[pr->ratePartition.size()+1];
        old_multi[0] = 1;
        int i= 1;
        bool cont = false;
        do {
            old_phi = pr->objective;
            old_rho = pr->rho;
            for (int r=1; r<=pr->ratePartition.size(); r++) {
                old_multi[r] = pr->multiplierRate[r];
            }
            for (int r=1; r<=pr->nbBranches; r++) {
                nodes[r]->B = B[r];
                nodes[r]->V = V[r];
            }
            calculateMultiplier(pr,nodes);
            for (int r=1; r<=pr->nbBranches; r++) {
                double m = pr->multiplierRate[nodes[r]->rateGroup];
                nodes[r]->B = B[r]/m;
                nodes[r]->V = V[r]/m/m;
            }
            val = without_constraint_active_set(false,pr,nodes);
            if (pr->verbose){
                for (int i=1;i<pr->multiplierRate.size();i++) cout<<pr->multiplierRate[i]<<" ";
                cout<<pr->rho<<endl;
            }
            cont = val && (abs((old_rho-pr->rho)/pr->rho) >= 1e-5);
            for (int r=1; r<=pr->ratePartition.size(); r++) {
                cont = cont || (pr->multiplierRate[r]>0 && abs((old_multi[r]*old_rho-pr->multiplierRate[r]*pr->rho)/pr->multiplierRate[r]/pr->rho)>=1e-5);
            }
            if (pr->verbose){
                printf("ROUND %d , objective function %.15e , rate %.15f , other_rates ",i,pr->objective,pr->rho);
                for (int r=1; r<=pr->ratePartition.size(); r++) {
                    printf(" %.15f",pr->rho*pr->multiplierRate[r]);
                }
                printf(", diff %.15f",abs((old_rho-pr->rho)/pr->rho));
                cout<<endl;
            }
            i++;
        } while (cont);
        if (!pr->haveUnique && (!pr->haveLower || pr->haveUpper)){
            val = without_constraint_active_set(pr,nodes,1);
            if (val){
                pr->rhoUpper = pr->rho;
            }
        }
        for (int i=1; i<=pr->nbBranches; i++) {
            nodes[i]->B = B[i];
            nodes[i]->V = V[i];
        }
    }
    delete[] B;
    delete[] V;
    return val;
}

bool with_constraint_active_set(bool all,Pr* pr,Node** nodes){
    bool val = false;
    if (pr->haveUnique) {
        val = with_constraint_active_set(pr,nodes,0);
    }
    if (!val){
        bool valL = false;
        bool valU = false;
        if (pr->haveLower){
            valL = with_constraint_active_set(pr,nodes,-1);
            if (valL){
                pr->rhoLower = pr->rho;
                val = true;
            }
        }
        if (!pr->haveLower || (all && pr->haveUpper)){
            valU = with_constraint_active_set(pr,nodes,1);
            if (valU){
                pr->rhoUpper = pr->rho;
                val = true;
            }
        }
    }
    return val;
}

bool with_constraint_multirates(Pr* pr,Node** nodes,bool reassign){
    double* B = new double[pr->nbBranches+1];
    double* V = new double[pr->nbBranches+1];
    for (int i=1; i<=pr->nbBranches; i++) {
        B[i] = nodes[i]->B;
        V[i] = nodes[i]->V;
    }
    if (pr->ratePartition.size()>0){
        if (reassign) assignRateGroupToTree(pr,nodes);
        for (int r=1; r<=pr->nbBranches; r++) {
            double m = pr->multiplierRate[nodes[r]->rateGroup];
            nodes[r]->B = B[r]/m;
            nodes[r]->V = V[r]/m/m;
        }
    }
    bool val = with_constraint_active_set(true,pr,nodes);
    if (!val) return false;
    if (pr->ratePartition.size()>0) {
        if (pr->verbose){
            printf("ROUND 0 , objective function %.15e , rate %.15f , other_rates ",pr->objective,pr->rho);
            for (int r=1; r<=pr->ratePartition.size(); r++) {
                printf(" %.15f",pr->rho*pr->multiplierRate[r]);
            }
            cout<<endl;
        }
        double old_phi = 0;
        double old_rho = 0;
        double* old_multi = new double[pr->ratePartition.size()+1];
        old_multi[0]=1;
        int i= 1;
        bool bl = false;
        bool cont = false;
        do {
            old_phi = pr->objective;
            old_rho = pr->rho;
            for (int r=1; r<=pr->ratePartition.size(); r++) {
                old_multi[r] = pr->multiplierRate[r];
            }
            for (int r=1; r<=pr->nbBranches; r++) {
                nodes[r]->B = B[r];
                nodes[r]->V = V[r];
            }
            calculateMultiplier(pr,nodes);
            for (int r=1; r<=pr->nbBranches; r++) {
                double m = pr->multiplierRate[nodes[r]->rateGroup];
                nodes[r]->B = B[r]/m;
                nodes[r]->V = V[r]/m/m;
            }
            if (!bl) {
                val = with_constraint_active_set(false,pr,nodes);
                bl = abs((old_rho-pr->rho)/pr->rho)<=1e-4;
                for (int r=1; r<=pr->ratePartition.size(); r++) {
                    bl = bl && (pr->multiplierRate[r]<0 || abs((old_multi[r]*old_rho-pr->multiplierRate[r]*pr->rho)/pr->multiplierRate[r]/pr->rho)<=1e-4);
                }
            }
            else {
                val = with_constraint_active_set(false,pr,nodes);
            }
            if (pr->verbose){
                for (int i=1;i<pr->multiplierRate.size();i++) cout<<pr->multiplierRate[i]<<" ";
                cout<<pr->rho<<endl;
            }
            cont = val && (abs((old_rho-pr->rho)/pr->rho) >= 1e-5);
            for (int r=1; r<=pr->ratePartition.size(); r++) {
                cont = cont || (pr->multiplierRate[r]>0 && abs((old_multi[r]*old_rho-pr->multiplierRate[r]*pr->rho)/pr->multiplierRate[r]/pr->rho)>=1e-5);
            }
            if (pr->verbose){
                printf("ROUND %d , objective function %.15e , rate %.15f , other_rates ",i,pr->objective,pr->rho);
                for (int r=1; r<=pr->ratePartition.size(); r++) {
                    printf(" %.15f",pr->rho*pr->multiplierRate[r]);
                }
                printf(", diff %.15f",abs((old_rho-pr->rho)/pr->rho));
                cout<<endl;
            }
            i++;
        } while (cont);
        if (!pr->haveUnique && (!pr->haveLower || pr->haveUpper)){
            val = with_constraint_active_set(pr,nodes,1);
            if (val){
                pr->rhoUpper = pr->rho;
            }
        }
        for (int r=1; r<=pr->nbBranches; r++) {
            nodes[r]->B = B[r];
            nodes[r]->V = V[r];
        }
    }
    delete[] B;
    delete[] V;
    return val;
}

