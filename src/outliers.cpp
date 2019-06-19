#include "outliers.h"

bool calculateOutliers(Pr* & pr,Node** & nodes){
    pr->outlier.clear();
    if (pr->partitionFile!="") {
        std::ostringstream oss;
        oss<<" - Rate partition can not be included in estimating outliers.\n";
        pr->warningMessage.push_back(oss.str());
    }
    if (pr->estimate_root==""){
        cout<<"Calculating the outlier tips ..."<<endl;
        outlier(pr,nodes,pr->k);
        if (pr->outlier.size()>0){
            std::ostringstream oss;
            oss<<"- There are "<<pr->outlier.size()<<" tips that are considered as outliers and were excluded from the analysis: ";
            for (int i=0;i<pr->outlier.size();i++){
                oss<<" "<<(nodes[pr->outlier[i]]->L).c_str();
            }
            oss<<"\n";
            pr->resultMessage.push_back(oss.str());
            bool bl = remove_outlier_tips(pr,nodes);
            if (!bl) {
                cout<<"Removing outliers will make the root lost. Check the root position and the dates of the outlier tips: ";
                for (int i=0;i<pr->outlier.size();i++){
                    cout<<(nodes[pr->outlier[i]]->L).c_str()<<" ";
                }
                cout<<". You can also try to not remove outliers or increase value of option -e to exclude some outliers."<<endl;
                return false;
            }
        }
        else{
            std::ostringstream oss;
            oss<<"- There is not any outlier tip.\n";
            pr->resultMessage.push_back(oss.str());
        }
    }
    else {
        cout<<"Calculating the outlier tips ..."<<endl;
        if (pr->estimate_root=="l"){
            estimate_root_local_rtt(pr,nodes);
        } else if (pr->estimate_root=="a" || pr->estimate_root=="as"){
            estimate_root_rtt(pr,nodes);
        }
        outlier(pr,nodes,pr->k);
        if (pr->outlier.size()>0){
            std::ostringstream oss;
            oss<<"- There are "<<pr->outlier.size()<<" tips that are considered as outliers and were excluded from the analysis: ";
            for (int i=0;i<pr->outlier.size();i++){
                oss<<" "<<(nodes[pr->outlier[i]]->L).c_str();
            }
            oss<<"\n";
            pr->resultMessage.push_back(oss.str());
            bool bl = remove_outlier_tips(pr,nodes);//cout<<pr->nbBranches<<" "<<pr->leaves<<" "<<pr->mrca<<" "<<pr->relative<<endl;
            if (!bl) {
                cout<<"Removing outliers will make the root lost. Check the root position and the dates of the outlier tips: ";
                for (int i=0;i<pr->outlier.size();i++){
                    cout<<(nodes[pr->outlier[i]]->L).c_str()<<" ";
                }
                cout<<". You can also try to not remove outliers or increase value of option -e to exclude some outliers."<<endl;
                return false;
            }
        }
        else{
            std::ostringstream oss;
            oss<<"- There is not any outlier tip.\n";
            pr->resultMessage.push_back(oss.str());
        }
    }
    return true;
}

double regression_lambda(double br,double &lambda,Pr* pr, Node** nodes){
    //double regression_lambda(int n,double* paths,double* paths_lambda,double* dates,double & lambda){
    /*double* dates = new double[pr->nbBranches-pr->nbINodes+1];
    for (int i=0;i<=pr->nbBranches-pr->nbINodes;i++) {
        dates[i] = nodes[i+pr->nbINodes]->D;
    }*/
    //double* paths = new double[pr->nbBranches-pr->nbINodes+1];
    //double* paths_lambda = new double[pr->nbBranches-pr->nbINodes+1];
    vector<double> paths;
    vector<double> dates;
    vector<double> paths_lambda;
    double slope = 0;
    double slope_lambda = 0;
    double intercept = 0;
    double intercept_lambda = 0;
    double det = 0;
    //for (int i=pr->nbINodes;i<=pr->nbBranches;i++) cout<<nodes[i]->D<<" ";
    //cout<<br<<endl;
    if (br == 0){
        vector<int> s = nodes[0]->suc;
        int s1 = s[0];
        int s2 = s[1];
        nodes[s1]->B = br/2;
        nodes[s2]->B = br/2;
        calculateRtt(pr,nodes,paths,dates);
        regression(pr,nodes,paths,dates,slope,intercept);
        double* res = residus_rtt(paths,dates,slope,intercept);
        double resSquare = 0;
        for (int i=0; i<= pr->nbBranches-pr->nbINodes; i++){
            resSquare += res[i]*res[i];
        }
        return resSquare;
    }
    calculateRtt_lambda(br,pr,nodes,paths,paths_lambda,dates);
    double mean_dates = 0;
    double mean_paths = 0;
    double mean_paths_lambda = 0;
    int n = dates.size();
    for (int i=0;i<n;i++){
        mean_dates += dates[i];
        mean_paths += paths[i];
        mean_paths_lambda += paths_lambda[i];
    }
    mean_dates /= n;
    mean_paths /= n;
    mean_paths_lambda /= n;
    
    if (pr->givenRate[0] && nodes[0]->type=='p'){
        slope = pr->rho;
        slope_lambda = 0;
        intercept = -slope*nodes[0]->D;
        intercept_lambda = 0;
    }
    else if (pr->relative){
        slope = (mean_paths)/(pr->leaves - pr->mrca);
        slope_lambda = (mean_paths_lambda)/(pr->leaves - pr->mrca);
        intercept = -slope*pr->mrca;
        intercept_lambda = -slope_lambda*pr->mrca;
    }
    else if (pr->givenRate[0]){
        slope = pr->rho;
        slope_lambda = 0;
        intercept =  mean_paths - slope*mean_dates;
        intercept_lambda = mean_paths_lambda;
    }
    else if (nodes[0]->type=='p'){
        intercept = nodes[0]->D;
        intercept_lambda = 0;
        for (int i=0;i<n;i++){
            det += dates[i]*dates[i];
            slope += dates[i]*(paths[i] - intercept);
            slope_lambda += dates[i]*(paths_lambda[i]);
        }
        slope = slope / det;
        slope_lambda = slope_lambda / det;
    }
    else{
        for (int i=0;i<n;i++){
            slope += (dates[i] - mean_dates)*(paths[i] - mean_paths);
            slope_lambda += (dates[i] - mean_dates)*(paths_lambda[i] - mean_paths_lambda);
            det += (dates[i] - mean_dates)*(dates[i] - mean_dates);
        }
        slope /= det;
        slope_lambda /= det;
        intercept =  mean_paths - slope*mean_dates;
        intercept_lambda = mean_paths_lambda - slope_lambda*mean_dates;
    }

    double A = 0;
    double B = 0;
    double C = 0;
    //square_residus = A*lambda^2 + B*lambda + C
    for (int i=0;i<n;i++){
        A += (paths_lambda[i] - slope_lambda*dates[i] - intercept_lambda)*(paths_lambda[i] - slope_lambda*dates[i] - intercept_lambda);
        B += 2*(paths_lambda[i] - slope_lambda*dates[i] - intercept_lambda)*(paths[i] - slope*dates[i] - intercept);
        C += (paths[i] - slope*dates[i] - intercept)*(paths[i] - slope*dates[i] - intercept);
    }
    lambda = -B/2/A;
    if (lambda <0) lambda = 0;
    if (lambda >1) lambda = 1;
    return (A*lambda*lambda + B*lambda + C);
}

void estimate_root_rtt(Pr* pr, Node** & nodes){
    //only apply for non-flexible given dates and non internal node date
    double phi1;
    Node** nodes_new = cloneLeaves(pr,nodes,0);
    int y=1;
    double l=0;
    double br=0;
    int r=0;
    double L=0;
    double BR=0;
    vector<int> s=nodes[0]->suc;
    int s1=s[0];
    int s2=s[1];
    vector<double> originalD;
    for (int i=0;i<=pr->nbBranches;i++) originalD.push_back(nodes[i]->D);
    bool bl=reroot_rootedtree(br,y,s1,s2,pr,nodes,nodes_new);
    for (int i=pr->nbINodes;i<=pr->nbBranches;i++){
        if (nodes_new[i]->type=='p') nodes_new[i]->D = originalD[i];
    }
    phi1 = regression_lambda(br,l,pr,nodes_new);
    r=y;
    y++;
    double phi;
    while (y<=pr->nbBranches){
        for (int i=pr->nbINodes; i<=pr->nbBranches; i++) {
            nodes_new[i]->status=nodes[i]->status;
        }
        bl=reroot_rootedtree(br,y,s1,s2,pr,nodes,nodes_new);
        for (int i=pr->nbINodes;i<=pr->nbBranches;i++){
            if (nodes_new[i]->type=='p') nodes_new[i]->D = originalD[i];
        }
        phi = regression_lambda(br,l,pr,nodes_new);
        if (phi1>phi || r==0){
            phi1=phi;
            r=y;
            L=l;
            BR=br;
        }
        y++;
    }
    nodes_new = cloneLeaves(pr,nodes,0);
    for (int i=pr->nbINodes; i<=pr->nbBranches; i++) {
        nodes_new[i]->status=nodes[i]->status;
    }
    reroot_rootedtree(BR,r,s1,s2,pr,nodes,nodes_new);
    nodes = nodes_new;
    s = nodes[0]->suc;
    s1 = s[0];
    s2 = s[1];
    nodes[s1]->B = L*BR;
    nodes[s2]->B = BR-L*BR;
}

void estimate_root_local_rtt(Pr* pr, Node** & nodes){
    //only apply for non-flexible given dates and non internal node date
    double phi1,phi;
    Node** nodes_new = cloneLeaves(pr,nodes,0);
    double* cv = new double[pr->nbBranches+1];
    for (int i=0;i<=pr->nbBranches;i++) cv[i]=0;
    double l=0;
    double br=0;
    int r=0;
    double L=0;
    double BR=0;
    vector<int> s=nodes[0]->suc;
    int s1=s[0];
    int s2=s[1];
    vector<double> originalD;
    for (int i=0;i<=pr->nbBranches;i++) originalD.push_back(nodes[i]->D);
    bool bl=reroot_rootedtree(br,s1,s1,s2,pr,nodes,nodes_new);
    for (int i=pr->nbINodes;i<=pr->nbBranches;i++){
        if (nodes_new[i]->type=='p') nodes_new[i]->D = originalD[i];
    }
    phi1 = regression_lambda(br,l,pr,nodes_new);
    cv[s1]=phi1;
    cv[s2]=cv[s1];
    phi1=cv[s1];
    r=s1;
    list<int> next;
    if (s1<pr->nbINodes){
        for (vector<int>::iterator iter=nodes[s1]->suc.begin(); iter!=nodes[s1]->suc.end(); iter++) {
            next.push_back(*iter);
        }
    }
    if (s2<pr->nbINodes){
        for (vector<int>::iterator iter=nodes[s2]->suc.begin(); iter!=nodes[s2]->suc.end(); iter++) {
            next.push_back(*iter);
        }
    }
    while (!next.empty()){
        int i = next.back();
        for (int i=pr->nbINodes; i<=pr->nbBranches; i++) {
            nodes_new[i]->status=nodes[i]->status;
        }
        bl=reroot_rootedtree(br,i,s1,s2,pr,nodes,nodes_new);
        for (int i=pr->nbINodes;i<=pr->nbBranches;i++){
            if (nodes_new[i]->type=='p') nodes_new[i]->D = originalD[i];
        }
        phi = regression_lambda(br,l,pr,nodes_new);
        cv[i]=phi;
        if (cv[i]<cv[nodes[i]->P] || r==0){
            if (i<pr->nbINodes){
                for (vector<int>::iterator iter=nodes[i]->suc.begin(); iter!=nodes[i]->suc.end(); iter++) {
                    next.push_back(*iter);
                }
            }
            if (cv[i]<phi1 || r==0){
                phi1=cv[i];
                r=i;
                L=l;
                BR=br;
            }
        }
        next.remove(i);
    }
    nodes_new = cloneLeaves(pr,nodes,0);
    for (int i=pr->nbINodes; i<=pr->nbBranches; i++) {
        nodes_new[i]->status=nodes[i]->status;
    }
    reroot_rootedtree(BR,r,s1,s2,pr,nodes,nodes_new);
    nodes = nodes_new;
    s = nodes[0]->suc;
    s1 = s[0];
    s2 = s[1];
    nodes[s1]->B = L*BR;
    nodes[s2]->B = BR-L*BR;
}

void regression(Pr* pr,Node** & nodes,vector<double> paths,vector<double> dates,double & slope,double & intercept){
    int n = paths.size();
    double mean_dates = 0;
    double mean_paths = 0;
    double det = 0;
    for (int i=0;i<n;i++){
        mean_dates += dates[i];
        mean_paths += paths[i];
    }
    mean_dates /= n;
    mean_paths /= n;
    slope = 0;
    if (pr->givenRate[0] && nodes[0]->type=='p'){
        slope = pr->rho;
        intercept = -slope*nodes[0]->D;
    }
    else if (pr->relative){
        slope = (mean_paths)/(pr->leaves - pr->mrca);
        intercept = -slope*pr->mrca;
    }
    else if (pr->givenRate[0]){
        slope = pr->rho;
        intercept =  mean_paths - slope*mean_dates;
    }
    else if (nodes[0]->type=='p'){
        intercept = nodes[0]->D;
        for (int i=0;i<n;i++){
            det += dates[i]*dates[i];
            slope += dates[i]*(paths[i] - intercept);
        }
        slope = slope / det;
    }
    else{
        for (int i=0;i<n;i++){
            slope += (paths[i] - mean_paths)*(dates[i] - mean_dates);
            det += (dates[i] - mean_dates)*(dates[i] - mean_dates);
        }
        slope = slope / det;
        intercept =  mean_paths - slope*mean_dates;
    }
}

void calculateRtt(Pr* pr,Node** nodes,vector<double> &paths,vector<double> &dates){
    for (int i=pr->nbINodes;i<=pr->nbBranches;i++){
        if (nodes[i]->type == 'p'){
            int k=i;
            double rtt = 0;
            while (k!=0){
                rtt += nodes[k]->B;
                k = nodes[k]->P;
            }
            paths.push_back(rtt);
            dates.push_back(nodes[i]->D);
        }
    }
}

void calculateRtt_lambda(double br,Pr* pr,Node** nodes,vector<double> & paths, vector<double> & paths_lambda,vector<double> & dates){
    //double* rtt = new double[pr->nbBranches-pr->nbINodes+1];
    //double* rtt_lambda = new double[pr->nbBranches-pr->nbINodes+1];
    vector<int> s = nodes[0]->suc;
    int s1 = s[0];
    int s2 = s[1];
    for (int i=pr->nbINodes;i<=pr->nbBranches;i++){
        if (nodes[i]->type == 'p'){
            int k=i;
            double rtt = 0;
            double rtt_lambda = 0;
            while (k != s1 && k!= s2){
                rtt += nodes[k]->B;
                k = nodes[k]->P;
            }
            if (k == s1){
                rtt_lambda = br;
            }
            if (k == s2){
                rtt += br;
                rtt_lambda = -br;
            }
            paths.push_back(rtt);
            paths_lambda.push_back(rtt_lambda);
            dates.push_back(nodes[i]->D);
        }
    }
}


double* residus_rtt(vector<double> paths,vector<double> dates,double slope,double intercept){
    int n = paths.size();
    double* res = new double[n];
    for (int i=0;i<n;i++){
        res[i] = (intercept + slope*dates[i]) - paths[i];
    }
    return res;
}

void outlier(Pr* pr,Node** nodes,double k){
    vector<double> paths, dates;
    string nodes_imprecise_date="";
    for (int i = pr->nbINodes; i <= pr->nbBranches; i++){
        if (nodes[i]->type != 'p'){
            nodes_imprecise_date = nodes_imprecise_date+" "+nodes[i]->L;
        }
    }
    if (nodes_imprecise_date.size()>0){
        std::ostringstream oss;
        oss<<" - The tip(s)"+nodes_imprecise_date+" do not have precise date value, so will not be included in estimating outliers.\n";
        pr->warningMessage.push_back(oss.str());
    }
    calculateRtt(pr,nodes,paths,dates);
    double slope,intercept;
    regression(pr,nodes,paths,dates,slope,intercept);
    double* res = residus_rtt(paths,dates,slope,intercept);
    double mi=0;
    double ma=0;
    double* sortedTab = sortTab(res,pr->nbBranches - pr->nbINodes + 1);
    getOulier(sortedTab,mi,ma,pr->nbBranches - pr->nbINodes + 1,k);
    delete[] sortedTab;
    pr->outlier.clear();
    for (int i=0;i<pr->nbBranches - pr->nbINodes + 1;i++){
        if (res[i]<mi || res[i]>ma){
            pr->outlier.push_back(i+pr->nbINodes);
        }
    }
    delete[] res;
    //check if outliers form a clan under the root
}


bool remove_one_tip(Pr* pr,Node** nodes,int t,int* &tab){//return false if remove outlier tip would remove the root
    int p = nodes[t]->P;
    if (p==0){
        return false;
    }
    tab[t] = -1;
    vector<int> suct;
    for (int i=0;i<nodes[p]->suc.size();i++){
        int ps = nodes[p]->suc[i];
        if (ps!=t) suct.push_back(ps);
    }
    if (suct.size()==1){
        int sibling_t = suct[0];
        tab[p] = -1;
        if (p>0){
            int pp = nodes[p]->P;
            nodes[sibling_t]->B = nodes[sibling_t]->B + nodes[p]->B;
            nodes[sibling_t]->P = pp;
            vector<int> sucpp;
            sucpp.push_back(sibling_t);
            for (int i=0;i<nodes[pp]->suc.size();i++){
                int ps = nodes[pp]->suc[i];
                if (ps!=p) sucpp.push_back(ps);
            }
            nodes[pp]->suc = sucpp;
            if (nodes[p]->type!='n'){
                nodes[sibling_t]->addConstraintOfP(nodes[p]);
                nodes[pp]->addConstraintOfS(nodes[p]);
            }
        }
        else{
            nodes[sibling_t]->P = -1;
        }
        
    }
    return true;
}

void shift_node_id(Pr* &pr,Node** &nodes,int* &tab){//t is a tip
    int shift = 0;
    int inodes_shift=0;
    for (int i=0; i<= pr->nbBranches;i++){
        if (tab[i]!=-1){
            tab[i] = (i-shift);
        }
        else {
            shift++;
            if (i<pr->nbINodes) inodes_shift++;
        }
    }
    Pr* prReduced = new Pr(pr->nbINodes-inodes_shift,pr->nbBranches-shift);
    prReduced->copy(pr);
    prReduced->outlier.clear();
    
    Node** nodesReduced = new Node*[prReduced->nbBranches+1];
    nodesReduced[0]=new Node();
    nodesReduced[0]->P=-1;
    nodesReduced[0]->type=nodes[0]->type;
    nodesReduced[0]->lower=nodes[0]->lower;
    nodesReduced[0]->upper=nodes[0]->upper;
    nodesReduced[0]->D=nodes[0]->D;
    for (int i=1; i<=pr->nbBranches; i++) {
        if (tab[i]!=-1) {
            nodesReduced[tab[i]]=new Node();
            nodesReduced[tab[i]]->P=tab[nodes[i]->P];
            nodesReduced[tab[i]]->B=nodes[i]->B;
            nodesReduced[tab[i]]->type=nodes[i]->type;
            nodesReduced[tab[i]]->lower=nodes[i]->lower;
            nodesReduced[tab[i]]->upper=nodes[i]->upper;
            nodesReduced[tab[i]]->D=nodes[i]->D;
            nodesReduced[tab[i]]->L=nodes[i]->L;
        }
    }
    if (pr->ratePartition.size()>0) {
        for (int i=0;i<=pr->nbBranches;i++){
            if (tab[i]!=-1) {
                nodesReduced[tab[i]]->rateGroup = nodes[i]->rateGroup;
            }
        }
    }
    for (int i=0;i<pr->internalConstraints.size();i++){
        Date* d = pr->internalConstraints[i];
        d->id = tab[d->id];
        prReduced->internalConstraints.push_back(d);
    }
    computeSuc_polytomy(prReduced, nodesReduced);
    computeVariance(prReduced,nodesReduced);
    nodes = nodesReduced;
    pr = prReduced;
}

bool remove_outlier_tips(Pr* &pr,Node** &nodes){
    int* tab = new int[pr->nbBranches+1];
    for (int i=0;i<=pr->nbBranches;i++) tab[i] = i;
    bool canRemove = true;
    int i=0;
    while (i<pr->outlier.size() && canRemove){
        canRemove = canRemove && remove_one_tip(pr,nodes,pr->outlier[i],tab);
        i++;
    }
    if (canRemove){
        shift_node_id(pr,nodes,tab);
        return true;
    }
    else{
        return false;
    }
}

/*double slope2(double* paths,double* dates,int i,int j){
 return (paths[j]-paths[i])/(dates[j]-dates[i]);
 }
 
 double intercept2(double* paths,double* dates,int i,int j){
 return (dates[j]*paths[i]-dates[i]*paths[j])/(dates[j]-dates[i]);
 }
 */

void getOulier(double* sortedArray,double& mi,double& ma,int size,double k){
    int q1 = round(size/4);
    int q3 = round(3*size/4);
    double iqr = sortedArray[q3-1]-sortedArray[q1-1];
    mi = sortedArray[q1-1]-k*iqr;
    ma = sortedArray[q3-1]+k*iqr;
}


