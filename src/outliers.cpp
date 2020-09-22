#include "outliers.h"
#include "dating.h"
#include "estimate_root.h"

bool calculateOutliers(Pr* & pr,Node** & nodes,double & median_rate,bool verbose){
    pr->outlier.clear();
    if (verbose){
        if (pr->partitionFile!="" || pr->splitExternal) {
            std::ostringstream oss;
            oss<<"- Multiple rates can not be included when estimating outliers.\n";
            pr->warningMessage.push_back(oss.str());
            cout<<" Multiple rates can not be included when estimating outliers"<<endl;
        }
        cout<<"Calculating the outlier nodes with Zscore threshold "<<pr->e<<" (setable via option -e)..."<<endl;
    }
    if (pr->estimate_root=="" || pr->estimate_root=="k"){
        bool givenRate = pr->givenRate[0];
        vector<double> dates_min;
        vector<double> dates_max;
        vector<int> index;
        bool both = false;
        vector<int> ignoredNodes;
        for (int i=pr->nbINodes;i<=pr->nbBranches;i++) {
            if (nodes[i]->type == 'p'){
                dates_min.push_back(nodes[i]->D);
                dates_max.push_back(nodes[i]->D);
            }
            if (nodes[i]->type == 'b'){
                both = true;
                dates_min.push_back(nodes[i]->lower);
                dates_max.push_back(nodes[i]->upper);
            }
            if (nodes[i]->type == 'u' || nodes[i]->type == 'l'){
                ignoredNodes.push_back(i);
            }
        }
        for (int i=0; i< pr->nbINodes;i++) {
            if (nodes[i]->type == 'p'){
                dates_min.push_back(nodes[i]->D);
                dates_max.push_back(nodes[i]->D);
            }
            if (nodes[i]->type == 'b') {
                both = true;
                dates_min.push_back(nodes[i]->lower);
                dates_max.push_back(nodes[i]->upper);
            }
            if (nodes[i]->type == 'u' || nodes[i]->type == 'l'){
                ignoredNodes.push_back(i);
            }
        }
        if (verbose && ignoredNodes.size()>0){
            cout<<"Ignore the nodes that only have upper/lower values in estimating outliers: ";
            std::ostringstream oss;
            oss<<" - The following nodes only have upper/lower values, so were ignored in estimating outliers: ";
            for (int i = 0; i<ignoredNodes.size();i++){
                cout<<nodes[ignoredNodes[i]]->L<<" ";
                oss<<nodes[ignoredNodes[i]]->L<<" ";
            }
            cout<<"\n";
            oss<<"\n";
            pr->warningMessage.push_back(oss.str());
        }
        for (int i=0;i<dates_min.size();i++){
            index.push_back(i);
        }
        int m = pr->m;
        if ((m +1) > index.size()){
            m = index.size() - 1;
            if (m<2){
                cout<<"Ignore estimating outliers: the temporal constraints provided are not enough, or conflict."<<endl;
                std::ostringstream oss;
                oss<<"- Ignore estimating outliers: the temporal constraints provided are not enough, or conflict.\n";
                pr->warningMessage.push_back(oss.str());
                return false;
            }
        }
        vector<int>* samples = new vector<int>[index.size()];
        for (int i = 0; i< index.size();i++){
            samples[i] = sampleNoRepeat(index, i , m);
        }
        if (pr->estimate_root==""){
            bool consistent;
            pr->outlier = outliers_rooted(pr,nodes,samples,dates_min,dates_max,false,both,median_rate,consistent);
            pr->givenRate[0] = givenRate;
            delete[] samples;
            if (!consistent){
                return false;
            }
        }
        if (pr->estimate_root=="k"){
            double br=0;
            int s1,s2;
            vector<int>::iterator iter=nodes[0]->suc.begin();
            s1=(*iter);
            iter++;
            s2=(*iter);
            br=nodes[s1]->B+nodes[s2]->B;
            bool consistent1, consistent2;
            double median_rate1, median_rate2;
            nodes[s1]->B = 0;
            nodes[s2]->B = br;
            vector<int> outliers1 = outliers_rooted(pr,nodes,samples,dates_min,dates_max,false,both,median_rate1,consistent1);
            pr->givenRate[0] = givenRate;
            nodes[s1]->B = br;
            nodes[s2]->B = 0;
            vector<int> outliers2 = outliers_rooted(pr,nodes,samples,dates_min,dates_max,false,both,median_rate2,consistent2);
            pr->givenRate[0] = givenRate;
            delete[] samples;
            if (consistent1 && consistent2){
                pr->outlier = intersect(outliers1,outliers2);
                median_rate = (median_rate1+median_rate2)/2;
            } else if (consistent1){
                pr->outlier = outliers1;
                median_rate = median_rate1;
            } else if (consistent2){
                pr->outlier = outliers2;
                median_rate = median_rate2;
            } else {
                cout<<"Ignore estimating outliers: the temporal constraints provided are not enough, or conflict."<<endl;
                std::ostringstream oss;
                oss<<"- Ignore estimating outliers: the temporal constraints provided are not enough, or conflict.\n";
                pr->warningMessage.push_back(oss.str());
                return false;
            }
        }
        if (pr->outlier.size()>0){
            std::ostringstream oss;
            oss<<"- The input dates associated with the following "<<pr->outlier.size()<<" nodes are considered as\n outliers, so those nodes were removed from the analysis:\n";
            for (int i=0;i<pr->outlier.size();i++){
                //if (pr->outlier[i] >= pr->nbINodes){
                    oss<<" "<<(nodes[pr->outlier[i]]->L).c_str();
                //} else {
                //    oss<<" "<<pr->internalConstraints[pr->outlier[i]]->label;
                //}
            }
            oss<<"\n";
            pr->resultMessage.push_back(oss.str());
            bool bl = remove_outlier_nodes(pr,nodes);
            if (!bl) {
                std::ostringstream oss;
                oss<<"- Removing outliers make the initial root lost. If you don't want that, you can\n try to reroot the tree or increase the Zscore threshold in option -e to exclude\n some outliers.\n"<<endl;
                pr->warningMessage.push_back(oss.str());
            }
            
        }
        else{
            std::ostringstream oss;
            oss<<"- There is not any outlier date.\n";
            pr->resultMessage.push_back(oss.str());
        }
    }
    else {
        bool bl = outliers_unrooted(pr,nodes,median_rate);
        if (bl){
            if (pr->outlier.size()>0){
                std::ostringstream oss;
                oss<<"- The input dates associated with the following "<<pr->outlier.size()<<" nodes are considered as\n outliers, so the nodes were excluded from the analysis:\n";
                for (int i=0;i<pr->outlier.size();i++){
                    //if (pr->outlier[i] >= pr->nbINodes){
                        oss<<" "<<(nodes[pr->outlier[i]]->L).c_str();
                    //} else {
                    //    oss<<" "<<pr->internalConstraints[pr->outlier[i]]->label;
                    //}
                }
                oss<<"\n";
                pr->resultMessage.push_back(oss.str());
                remove_outlier_nodes(pr,nodes);
            }
            else{
                std::ostringstream oss;
                oss<<"- There is not any outlier date.\n";
                pr->resultMessage.push_back(oss.str());
            }
        } else {
            cout<<"Ignore estimating outliers: the temporal constraints are not enough or conflict."<<endl;
            std::ostringstream oss;
            oss<<"- Ignore estimating outliers: the temporal constraints are not enough or conflict.\n";
            pr->warningMessage.push_back(oss.str());
            return false;
        }
    }
    return true;
}


bool remove_one_tip(Pr* pr,Node** nodes,int t,int* &tab){//return false if remove outlier tip would remove the root
    tab[t] = -1;
    int p = nodes[t]->P;
    vector<int> suct;
    bool bl = true;
    for (int i=0;i<nodes[p]->suc.size();i++){
        int ps = nodes[p]->suc[i];
        if (ps!=t) {
            suct.push_back(ps);
        }
    }
    if (suct.size()==1){
        int sibling_t = suct[0];
        if (p==0 && (nodes[sibling_t]->suc.size()==2 || pr->estimate_root=="" || pr->estimate_root=="k")) {
            tab[sibling_t] = -1;
            nodes[0]->L = nodes[sibling_t]->L;
            nodes[0]->V = nodes[sibling_t]->V;
            nodes[0]->type = nodes[sibling_t]->type;
            nodes[0]->lower = nodes[sibling_t]->lower;
            nodes[0]->upper = nodes[sibling_t]->upper;
            nodes[0]->D = nodes[sibling_t]->D;
            nodes[0]->rateGroup = nodes[sibling_t]->rateGroup;
            nodes[0]->status = nodes[sibling_t]->status;
            nodes[0]->suc.clear();
            for (int i=0;i<nodes[sibling_t]->suc.size();i++){
                nodes[nodes[sibling_t]->suc[i]]->P = 0;
                nodes[0]->suc.push_back(nodes[sibling_t]->suc[i]);
            }
            bl = false;
        } else if (p==0){
            int s = nodes[sibling_t]->suc[0];
            double br=nodes[s]->B;
            nodes[s]->B=br/2;
            nodes[sibling_t]->B=br/2;
            nodes[s]->P=0;
            nodes[sibling_t]->P=0;
            nodes[sibling_t]->suc.erase(nodes[sibling_t]->suc.begin());
            nodes[0]->suc.clear();
            nodes[0]->suc.push_back(sibling_t);
            nodes[0]->suc.push_back(s);
        }
        else {
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
    } else{
        nodes[p]->suc = suct;
    }
    return bl;
}


void shift_node_id(Pr* &pr,Node** &nodes,int* &tab){
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
    
    Node** nodesReduced = new Node*[prReduced->nbBranches+1];
    nodesReduced[0]=new Node();
    nodesReduced[0]->P=-1;
    nodesReduced[0]->type=nodes[0]->type;
    nodesReduced[0]->lower=nodes[0]->lower;
    nodesReduced[0]->upper=nodes[0]->upper;
    nodesReduced[0]->D=nodes[0]->D;
    nodesReduced[0]->L=nodes[0]->L;
    for (int i=1; i<=pr->nbBranches; i++) {
        if (tab[i]!=-1) {
            nodesReduced[tab[i]]=new Node();
            nodesReduced[tab[i]]->P=tab[nodes[i]->P];
            nodesReduced[tab[i]]->B=nodes[i]->B;
            nodesReduced[tab[i]]->type=nodes[i]->type;
            nodesReduced[tab[i]]->status=nodes[i]->status;
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
        if (tab[d->id] != -1){
            int m = d->mrca.size();
            if (m>0 && pr->estimate_root!="" && pr->estimate_root!="k"){
                vector<int> ignoredNodes;
                for (int j=0;j<d->mrca.size();j++){
                    if (tab[d->mrca[j]]==-1){
                        d->mrca.erase(d->mrca.begin()+j);
                        ignoredNodes.push_back(d->mrca[j]);
                    } else {
                        d->mrca[j] = tab[d->mrca[j]];
                    }
                }
                if (m>=2 && d->mrca.size()<2){
                    ostringstream ss;
                    ss<<"- The input date constraint of the node "<<d->label<<" was ignored due to outliers removal.\n";
                    prReduced->warningMessage.push_back(ss.str());
                } else if (ignoredNodes.size()>0){
                    ostringstream ss;
                    ss<<"- The tips ";
                    for (int j=0;j<ignoredNodes.size();j++){
                        ss<<nodes[ignoredNodes[j]]->L<<" ";
                    }
                    ss<<"are outliers and were removed, that might effect the input date constraint of the node "<<d->label<<"\n";
                    prReduced->warningMessage.push_back(ss.str());
                    nodesReduced[tab[d->id]]->type = 'n';
                    nodesReduced[tab[d->id]]->status = 0;
                } else {
                    d->id = tab[d->id];
                    prReduced->internalConstraints.push_back(d);
                }
            } else {
                d->id = tab[d->id];
                d->mrca.clear();
                prReduced->internalConstraints.push_back(d);
            }
        }
    }
    computeSuc_polytomy(prReduced, nodesReduced);
    computeVariance(prReduced,nodesReduced);
    nodes = nodesReduced;
    pr = prReduced;
    initConstraint(pr, nodes);
}

bool remove_outlier_nodes(Pr* &pr,Node** &nodes){
    int* tab = new int[pr->nbBranches+1];
    for (int i=0;i<=pr->nbBranches;i++) tab[i] = i;
    bool keepRoot = true;
    int i=0;
    //int shift = 0;
    while (i<pr->outlier.size()){
        //if (pr->outlier[i] >= pr->nbINodes){
            keepRoot = remove_one_tip(pr,nodes,pr->outlier[i],tab) && keepRoot;
        //} else {//ignored the input date of the node
        //    pr->internalConstraints.erase(pr->internalConstraints.begin()+pr->outlier[i] - shift);
        //    shift++;
        //}
        i++;
    }
    shift_node_id(pr,nodes,tab);
    delete[] tab;
    return keepRoot;
}

vector<int> sampleNoRepeat(vector<int> array,int ex, int size){
    swap(array[ex],array[array.size()-1]);
    vector<int> result;
    for (int i=0; i< size;i++){
        int a = std::rand() % (array.size()-1-i);
        result.push_back(array[a]);
        swap(array[a],array[array.size()-2-i]);
    }
    return result;
}

vector<double> residus_lsd(Pr* pr,Node** nodes,double& mean_res, double& var_res){
    vector<double> res;
    for (int i=1;i<=pr->nbBranches;i++){
        double r = (nodes[i]->B - pr->rho*(nodes[i]->D - nodes[nodes[i]->P]->D))/sqrt(nodes[i]->V);
        res.push_back(r);
        mean_res += r;
    }
    mean_res /= pr->nbBranches;
    for (int i=0;i<pr->nbBranches;i++){
        var_res += (res[i] - mean_res)*(res[i] - mean_res);
    }
    var_res /= (pr->nbBranches-1);
    return res;
}

vector<double> residus_lsd_rtt(Pr* pr,Node** nodes,double& mean_res, double& var_res){
    vector<double> res;
    for (int i=1;i<=pr->nbBranches;i++){
        double r = (nodes[i]->B - pr->rho*(nodes[i]->D - nodes[nodes[i]->P]->D))/sqrt(nodes[i]->V);
        res.push_back(r);
    }
    return res;
}

void calculateRoot2DatedNode(Pr* pr,Node** nodes,vector<double> &paths){
    for (int i=pr->nbINodes;i<=pr->nbBranches;i++){
        if (nodes[i]->type == 'p' || nodes[i]->type == 'b'){
            int k=i;
            double rtt = 0;
            while (k!=0){
                rtt += nodes[k]->B;
                k = nodes[k]->P;
            }
            paths.push_back(rtt);
        }
    }
    for (int i=0;i<pr->nbINodes;i++){
        if (nodes[i]->type == 'p' || nodes[i]->type == 'b'){
            int k=i;
            double rtt = 0;
            while (k!=0){
                rtt += nodes[k]->B;
                k = nodes[k]->P;
            }
            paths.push_back(rtt);
        }
    }
}

void calculateRoot2DatedNode(Pr* pr,Node** nodes,vector<double> &paths,vector<double> &dates){
    for (int i=pr->nbINodes;i<=pr->nbBranches;i++){
        if (nodes[i]->type == 'p'){
            int k=i;
            double rtt = 0;
            while (k!=0){
                rtt += nodes[k]->B;
                k = nodes[k]->P;
            }
            paths.push_back(rtt);
        }
    }
    vector<int> visited;
    for (int i=0;i<pr->internalConstraints.size();i++){
        Date* no = pr->internalConstraints[i];
        if ((no->type=='p') && !isIn(no->id,visited)){
            int k=no->id;
            double rtt = 0;
            while (k!=0){
                rtt += nodes[k]->B;
                k = nodes[k]->P;
            }
            paths.push_back(rtt);
            dates.push_back(nodes[no->id]->D);
            visited.push_back(no->id);
        }
    }
}

void calculateRoot2DatedNode(Pr* pr,Node** nodes,vector<double> &paths,vector<double> &dates_min,vector<double> &dates_max){
    for (int i=pr->nbINodes;i<=pr->nbBranches;i++){
        if (nodes[i]->type == 'p' || nodes[i]->type == 'b'){
            int k=i;
            double rtt = 0;
            while (k!=0){
                rtt += nodes[k]->B;
                k = nodes[k]->P;
            }
            paths.push_back(rtt);
        }
    }
    vector<int> visited;
    for (int i=0;i<pr->internalConstraints.size();i++){
        Date* no = pr->internalConstraints[i];
        if ((no->type=='p' || no->type=='b') && !isIn(no->id,visited)){
            int k=no->id;
            double rtt = 0;
            while (k!=0){
                rtt += nodes[k]->B;
                k = nodes[k]->P;
            }
            paths.push_back(rtt);
            if (no->type == 'p'){
                dates_min.push_back(nodes[no->id]->D);
                dates_max.push_back(nodes[no->id]->D);
            } else {
                dates_min.push_back(nodes[no->id]->lower);
                dates_max.push_back(nodes[no->id]->upper);
            }
            visited.push_back(no->id);
        }
    }
}

bool median_rate(Pr* pr,Node** nodes, vector<double> dates, vector<int>* samples,bool addInternalDates,double& rate){
    if (pr->givenRate[0]){
        rate = pr->rho;
        return true;
    }
    else{
        vector<double> paths;
        if (addInternalDates) calculateRoot2DatedNode(pr,nodes,paths,dates);
        else calculateRoot2DatedNode(pr,nodes,paths);
        vector<double> rates;
        for (int i = 0; i< paths.size();i++){
            vector<double> slopes;
            for (vector<int>::iterator iter = samples[i].begin(); iter != samples[i].end(); iter++){
                if ((*iter < paths.size()) && dates[*iter] != dates[i]){
                    double slope = (paths[*iter] - paths[i])/(dates[*iter] - dates[i]);
                    slopes.push_back(slope);
                }
            }
            if (!slopes.empty()){
                double m = median(slopes);
                if (m>0) rates.push_back(m);
            }
        }
        if (rates.size()==0) return false;
        else {
            rate = median(rates);
            return true;
        }
    }
}


bool median_rate(Pr* pr,Node** nodes, vector<double> dates_min,vector<double> dates_max, vector<int>* samples,bool addInternalDates,double& rate1, double& rate2){
    if (pr->givenRate[0]){
        rate1 = pr->rho;
        rate2 = pr->rho;
        return true;
    }
    else{
        vector<double> paths;
        if (addInternalDates) calculateRoot2DatedNode(pr,nodes,paths,dates_min,dates_max);
        else calculateRoot2DatedNode(pr,nodes,paths);
        vector<double> rates_min;
        vector<double> rates_max;
        for (int i = 0; i< paths.size();i++){
            vector<double> slopeMin;
            vector<double> slopeMax;
            for (vector<int>::iterator iter = samples[i].begin(); iter != samples[i].end(); iter++){
                if ((*iter < paths.size()) && dates_min[*iter] != dates_min[i]){
                    double slope_min = (paths[*iter] - paths[i])/(dates_min[*iter] - dates_min[i]);
                    slopeMin.push_back(slope_min);
                }
                if ((*iter < paths.size()) && dates_max[*iter] != dates_max[i]){
                    double slope_max = (paths[*iter] - paths[i])/(dates_max[*iter] - dates_max[i]);
                    slopeMax.push_back(slope_max);
                }
            }
            if (!slopeMin.empty()){
                double m = median(slopeMin);
                if (m>0) rates_min.push_back(m);
            }
            if (!slopeMax.empty()){
                double m = median(slopeMax);
                if (m>0) rates_max.push_back(median(slopeMax));
            }
        }
        if (rates_min.size()>0 && rates_max.size()>0){
            rate1 = median(rates_min);
            rate2 = median(rates_max);
            return true;
        } else {
            if (rates_max.size()>0){
                rate2 = median(rates_max);
                return true;
            }
            else if (rates_min.size()>0){
                rate1 = median(rates_min);
                return true;
            }
            else return false;
        }
    }
}

vector<int> outliers_rooted(Pr* pr,Node** nodes,vector<int>* samples, vector<double> dates_min, vector<double> dates_max, bool addInternalDates, bool both,double & med_rate,bool& consistent){
    double rate_min, rate_max;
    if (!both) {
        consistent = median_rate(pr,nodes,dates_min, samples,addInternalDates, rate_min);
        med_rate = rate_min;
    } else {
        consistent = median_rate(pr,nodes,dates_min, dates_max,samples,addInternalDates, rate_min, rate_max);
        med_rate = (rate_min+rate_max)/2;
    }
    vector<int> outliers_min;
    if (!consistent) return outliers_min;
    pr->givenRate[0] = true;
    pr->rho = rate_min;
    without_constraint_multirates(pr,nodes,true); 
    double mean_res = 0;
    double var_res = 0;
    vector<double> res = residus_lsd(pr,nodes,mean_res,var_res);
    for (int i=0;i<res.size();i++){
        res[i] = (res[i]-mean_res)/sqrt(var_res);
    }
    for (int i=pr->nbINodes;i<=pr->nbBranches;i++){
        if (abs(res[i-1])>pr->e){
            if (nodes[i]->type == 'p' || nodes[i]->type == 'b'){
                outliers_min.push_back(i);
            }
        }
    }
    /*for (int k=0;k<pr->internalConstraints.size();k++){
        Date* no = pr->internalConstraints[k];
        int i = no->id;
        if ((!pr->relative || i!=0) && (nodes[i]->type == 'p' || nodes[i]->type == 'b')){
                bool bl = false;
                if (i>0){
                    bl = (abs(res[i-1])>pr->e);
                }
                for (int j=0;j<nodes[i]->suc.size();j++){
                    int s = nodes[i]->suc[j];
                    bl = bl || (abs(res[s-1])>pr->e);
                }
                if (bl){
                    outliers_min.push_back(k);
                }
        }
    }*/
    if (!both) return outliers_min;
    else {
        pr->rho = rate_max;
        without_constraint_multirates(pr,nodes,true);
        mean_res = 0;
        var_res = 0;
        res = residus_lsd(pr,nodes,mean_res,var_res);
        for (int i=0;i<res.size();i++) {
            res[i] = (res[i]-mean_res)/sqrt(var_res);
        }
        vector<int> outliers_max;
        for (int i=pr->nbINodes;i<=pr->nbBranches;i++){
            if (abs(res[i-1])>pr->e){
                if (nodes[i]->type == 'p' || nodes[i]->type == 'b'){
                    outliers_max.push_back(i);
                }
            }
        }
        /*for (int k=0;k<pr->internalConstraints.size();k++){
            Date* no = pr->internalConstraints[k];
            int i = no->id;
            if (nodes[i]->type == 'p' || nodes[i]->type == 'b'){
                bool bl = false;
                if (i>0){
                    bl = (abs(res[i-1])>pr->e);
                }
                for (int j=0;j<nodes[i]->suc.size();j++){
                    int s = nodes[i]->suc[j];
                    bl = bl || (abs(res[s-1])>pr->e);
                }
                if (bl){
                    outliers_max.push_back(k);
                }
            }
        }*/
        return intersect(outliers_min,outliers_max);
    }
}

bool outliers_unrooted(Pr* &pr,Node** &nodes,double& median_rate){
    Node** nodes_new = cloneLeaves(pr,nodes,0);
    double br=0;
    vector<int>::iterator iter=nodes[0]->suc.begin();
    int s1=(*iter);
    iter++;
    int s2=(*iter);
    vector<double> originalD_min;
    vector<double> originalD_max;
    int nbFixedNodes = 0;
    bool both = false;
    for (int i=pr->nbINodes;i<=pr->nbBranches;i++) {
        if (nodes[i]->type == 'p') {
            originalD_min.push_back(nodes[i]->D);
            originalD_max.push_back(nodes[i]->D);
            nbFixedNodes++;
        }
        if (nodes[i]->type == 'b') {
            originalD_min.push_back(nodes[i]->lower);
            originalD_max.push_back(nodes[i]->upper);
            nbFixedNodes++;
            both = true;
        }
    }
    bool bl = false;
    int y=1;
    while (!bl && y<=pr->nbBranches){
        bl=reroot_rootedtree(br,y,s1,s2,pr,nodes,nodes_new);
        y++;
    }
    bool consistent;
    if (bl){
        for (int i=0; i< pr->internalConstraints.size();i++){
            if (pr->internalConstraints[i]->type == 'p' || pr->internalConstraints[i]->type == 'b') nbFixedNodes++;
        }
        vector<int> index;
        for (int i=0;i<nbFixedNodes;i++){
            index.push_back(i);
        }
        int m = pr->m;
        if (m > (index.size()-1)){
            m = index.size()-1;
            if (m<2){
                cout<<"Ignore estimating outliers: the temporal constraints are not enough or conflict."<<endl;
                std::ostringstream oss;
                oss<<"- Ignore estimating outliers: the temporal constraints are not enough or conflict.\n";
                pr->warningMessage.push_back(oss.str());
                return false;
            }
        }
        vector<int>* samples = new vector<int>[index.size()];
        for (int i = 0; i< index.size();i++){
            samples[i] = sampleNoRepeat(index, i , m);
        }
        vector<int> outliersFreq;
        for (int i=0;i<=pr->nbBranches;i++){
            outliersFreq.push_back(0);
        }
        nodes_new[y]->B = 0;
        nodes_new[nodes[y]->P]->B = br;
        bool givenRate = pr->givenRate[0];
        vector<double> median_rates;
        double med_rate;
        vector<int> outs=outliers_rooted(pr,nodes_new,samples,originalD_min,originalD_max,true,both,med_rate,consistent);
        if (consistent){
            median_rates.push_back(med_rate);
            for (int i=0;i<outs.size();i++){
                outliersFreq[outs[i]] ++;
            }
        }
        pr->givenRate[0] = givenRate;
        y++;
        while (y<=pr->nbBranches){
            bl=reroot_rootedtree(br,y,s1,s2,pr,nodes,nodes_new);
            if (bl){
                nodes_new[y]->B = 0;
                nodes_new[nodes[y]->P]->B = br;
                vector<int> outs=outliers_rooted(pr,nodes_new,samples,originalD_min,originalD_max,true,both,med_rate,consistent);
                if (consistent){
                    median_rates.push_back(med_rate);
                    for (int i=0;i<outs.size();i++){
                        outliersFreq[outs[i]] ++;
                    }
                }
            }
            pr->givenRate[0] = givenRate;
            y++;
        }
        if (median_rates.empty()){
            return false;
        }
        median_rate = median(median_rates);
        delete[] samples;
        for (int i=0;i<=pr->nbBranches;i++) delete nodes_new[i];
        delete[] nodes_new;
        pr->outlier.clear();
        double mean_freq = 0;
        double var_freq = 0;
        for (int i=1;i<=pr->nbBranches;i++){
            mean_freq += outliersFreq[i];
        }
        mean_freq /= pr->nbBranches;
        for (int i=0;i<pr->nbBranches;i++){
            var_freq += (outliersFreq[i] - mean_freq)*(outliersFreq[i] - mean_freq);
        }
        var_freq /= (pr->nbBranches-1);
        
        for (int i=0;i<outliersFreq.size();i++){
            outliersFreq[i] = (outliersFreq[i] - mean_freq)/sqrt(var_freq);
            if ((outliersFreq[i])>pr->e){
                pr->outlier.push_back(i);
            }
        }
        return true;
    } else {
        return false;
    }
}



/*
bool calculateMedianRate(Pr* pr,Node** nodes,double& med_rate){
    double rate_min, rate_max;
    if (pr->estimate_root=="" || pr->estimate_root=="k"){
        bool givenRate = pr->givenRate[0];
        vector<double> dates_min;
        vector<double> dates_max;
        vector<int> index;
        bool both = false;
        for (int i=pr->nbINodes;i<=pr->nbBranches;i++) {
            if (nodes[i]->type == 'p'){
                dates_min.push_back(nodes[i]->D);
                dates_max.push_back(nodes[i]->D);
            }
            if (nodes[i]->type == 'b'){
                both = true;
                dates_min.push_back(nodes[i]->lower);
                dates_max.push_back(nodes[i]->upper);
            }
        }
        for (int i=0; i< pr->nbINodes;i++) {
            if (nodes[i]->type == 'p'){
                dates_min.push_back(nodes[i]->D);
                dates_max.push_back(nodes[i]->D);
            }
            if (nodes[i]->type == 'b') {
                both = true;
                dates_min.push_back(nodes[i]->lower);
                dates_max.push_back(nodes[i]->upper);
            }
        }
        for (int i=0;i<dates_min.size();i++){
            index.push_back(i);
        }
        int m = pr->m;
        if ((m +1) > index.size()){
            m = index.size() - 1;
            if (m<2) return false;
        }
        vector<int>* samples = new vector<int>[index.size()];
        for (int i = 0; i< index.size();i++){
            samples[i] = sampleNoRepeat(index, i , m);
        }
        bool consistent;
        if (pr->estimate_root==""){
            if (!both) {
                consistent = median_rate(pr,nodes,dates_min, samples,false, med_rate);
            } else {
                consistent = median_rate(pr,nodes,dates_min, dates_max,samples,false, rate_min, rate_max);
                med_rate = (rate_min+rate_max)/2;
            }
        }
        if (pr->estimate_root=="k"){
            double br=0;
            int s1,s2;
            vector<int>::iterator iter=nodes[0]->suc.begin();
            s1=(*iter);
            iter++;
            s2=(*iter);
            br=nodes[s1]->B+nodes[s2]->B;
            double median_rate1, median_rate2;
            nodes[s1]->B = 0;
            nodes[s2]->B = br;
            if (!both) {
                consistent = median_rate(pr,nodes,dates_min, samples,false, median_rate1);
            } else {
                consistent = median_rate(pr,nodes,dates_min, dates_max,samples,false, rate_min, rate_max);
                median_rate1 = (rate_min+rate_max)/2;
            }
            nodes[s1]->B = br;
            nodes[s2]->B = 0;
            if (!both) {
                consistent = consistent && median_rate(pr,nodes,dates_min, samples,false, median_rate2);
            } else {
                consistent = consistent && median_rate(pr,nodes,dates_min, dates_max,samples,false, rate_min, rate_max);
                median_rate2 = (rate_min+rate_max)/2;
            }
            med_rate = (median_rate1+median_rate2)/2;
        }
        pr->givenRate[0] = givenRate;
        delete[] samples;
        return consistent;
    }
    else {
        Node** nodes_new = cloneLeaves(pr,nodes,0);
        double br=0;
        vector<int>::iterator iter=nodes[0]->suc.begin();
        int s1=(*iter);
        iter++;
        int s2=(*iter);
        vector<double> originalD_min;
        vector<double> originalD_max;
        int nbFixedNodes = 0;
        bool both = false;
        for (int i=pr->nbINodes;i<=pr->nbBranches;i++) {
            if (nodes[i]->type == 'p') {
                originalD_min.push_back(nodes[i]->D);
                originalD_max.push_back(nodes[i]->D);
                nbFixedNodes++;
            }
            if (nodes[i]->type == 'b') {
                originalD_min.push_back(nodes[i]->lower);
                originalD_max.push_back(nodes[i]->upper);
                nbFixedNodes++;
                both = true;
            }
        }
        bool bl = false;
        int y=1;
        while (!bl && y<=pr->nbBranches){
            bl=reroot_rootedtree(br,y,s1,s2,pr,nodes,nodes_new);
            y++;
        }
        if (bl){
            for (int i=0; i< pr->internalConstraints.size();i++){
                if (pr->internalConstraints[i]->type == 'p' || pr->internalConstraints[i]->type == 'b') nbFixedNodes++;
            }
            vector<int> index;
            for (int i=0;i<nbFixedNodes;i++){
                index.push_back(i);
            }
            int m = pr->m;
            if (m > (index.size()-1)){
                m = index.size()-1;
                if (m<2) return false;
            }
            vector<int>* samples = new vector<int>[index.size()];
            for (int i = 0; i< index.size();i++){
                samples[i] = sampleNoRepeat(index, i , m);
            }
            nodes_new[y]->B = 0;
            nodes_new[nodes[y]->P]->B = br;
            bool givenRate = pr->givenRate[0];
            vector<double> median_rates;
            double mrate;
            bool consistent;
            if (!both) {
                 consistent = median_rate(pr,nodes_new,originalD_min,samples,true, mrate);
            } else {
                consistent = median_rate(pr,nodes_new,originalD_min,originalD_max, samples,true, rate_min, rate_max);
                mrate = (rate_min+rate_max)/2;
            }
            if (consistent){
                median_rates.push_back(mrate);
            }
            pr->givenRate[0] = givenRate;
            y++;
            while (y<=pr->nbBranches){
                bl=reroot_rootedtree(br,y,s1,s2,pr,nodes,nodes_new);
                if (bl){
                    nodes_new[y]->B = 0;
                    nodes_new[nodes[y]->P]->B = br;
                    if (!both) {
                        consistent = median_rate(pr,nodes_new,originalD_min,samples,true, mrate);
                    } else {
                        consistent = median_rate(pr,nodes_new,originalD_min,originalD_max, samples,true, rate_min, rate_max);
                        mrate = (rate_min+rate_max)/2;
                    }
                    if (consistent) median_rates.push_back(mrate);
                }
                pr->givenRate[0] = givenRate;
                y++;
            }
            delete[] samples;
            for (int i=0;i<=pr->nbBranches;i++) delete nodes_new[i];
            delete[] nodes_new;
            if (median_rates.empty()){
                return false;
            } else {
                med_rate = median(median_rates);
                return true;
            }
        } else {
            return false;
        }
    }
}


double rate_rtt(Pr* pr,Node** & nodes){
    vector<double> paths;
    vector<double> dates;
    calculateRoot2DatedNode(pr,nodes,paths,dates);
    double slope;
    double intercept;
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
        intercept = nodes[0]->D;
    }
    else if (pr->relative){
        slope = (mean_paths)/(pr->leaves - pr->mrca);
        intercept = pr->mrca;
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
    return slope;
}

void calculateRoot2DatedNode_lambda(double br,Pr* pr,Node** nodes,vector<double> & paths, vector<double> & paths_lambda,vector<double> & dates){
    vector<int> s = nodes[0]->suc;
    int s1 = s[0];
    int s2 = s[1];
    for (int i=pr->nbINodes;i<=pr->nbBranches;i++){
        if (nodes[i]->type == 'p'){
            int k=i;
            double rtt = 0;
            double rtt_lambda = 0;
            while (k!=s1 && k!=s2){
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
    vector<int> visited;
    for (int i=0;i<pr->internalConstraints.size();i++){
        Date* no = pr->internalConstraints[i];
        if ((no->type=='p') && !isIn(no->id,visited)){
            int k=no->id;
            double rtt = 0;
            double rtt_lambda = 0;
            while (k!=s1 && k!=s2){
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
            visited.push_back(no->id);
        }
    }
}


double regression_lambda(double br,Pr* pr, Node** nodes){
    double lambda;
    vector<double> paths;
    vector<double> dates;
    vector<double> paths_lambda;
    double slope = 0;
    double slope_lambda = 0;
    double intercept = 0;
    double intercept_lambda = 0;
    double det = 0;
    if (br == 0){
        return rate_rtt(pr,nodes);
    }
    calculateRoot2DatedNode_lambda(br,pr,nodes,paths,paths_lambda,dates);
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
    return (slope + lambda*slope_lambda);//or (slope + slope_lambda/lambda);
}

double median_rate_rtt(Pr* pr, Node** & nodes){
    //only apply for non-flexible given dates and non internal node date
    Node** nodes_new = cloneLeaves(pr,nodes,0);
    int y=1;
    double l=0;
    double br=0;
    vector<int> s=nodes[0]->suc;
    int s1=s[0];
    int s2=s[1];
    vector<double> originalD;
    for (int i=0;i<=pr->nbBranches;i++) originalD.push_back(nodes[i]->D);
    bool bl=reroot_rootedtree(br,y,s1,s2,pr,nodes,nodes_new);
    for (int i=pr->nbINodes;i<=pr->nbBranches;i++){
        if (nodes_new[i]->type=='p') nodes_new[i]->D = originalD[i];
    }
    vector<double> slopes;
    double slope = regression_lambda(br,pr,nodes_new);
    slopes.push_back(slope);
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
        slope = regression_lambda(br,pr,nodes_new);
        slopes.push_back(slope);
        y++;
    }
    return median(slopes);
}
*/
