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

#include "confidence_interval.h"


void computeIC(double br,Pr* pr,Node** nodes,double* &T_left,double* &T_right,double* &H_left,double* &H_right,double* &HD_left,double* &HD_right,double &rho_left,double& rho_right,double* &other_rhos_left,double* &other_rhos_right){
    Node** nodes_new = new Node*[pr->nbBranches+1];
    int* tab = new int[pr->nbBranches+1];
    bool useSupport = false;
    int nbC = collapseTree(pr, nodes, nodes_new,tab,0,useSupport);//nbC is the number of internal nodes reduced
    Node** nodesReduced = new Node*[nbC+pr->nbBranches-pr->nbINodes+1];
    Pr* prReduced = new Pr(nbC,nbC+pr->nbBranches-pr->nbINodes);
    prReduced->copy(pr);
    collapseTreeReOrder( pr, nodes_new, prReduced, nodesReduced,tab);
    initConstraint(prReduced, nodesReduced);
    for (int i=0;i<pr->nbBranches+1;i++){
        delete nodes_new[i];
    }
    delete[] nodes_new;
    computeSuc_polytomy(prReduced, nodesReduced);
    double** B_simul = new double*[pr->nbSampling];
    for (int i=0;i<pr->nbSampling;i++){
        B_simul[i]=new double[prReduced->nbBranches+1];
    }
    std::default_random_engine generator;
    
    double minB=nodesReduced[1]->B;
    for (int i=2;i<=prReduced->nbBranches;i++){
        if (nodesReduced[i]->B<minB && nodesReduced[i]->B>0) {
            minB=nodesReduced[i]->B;
        }
    }
    int seqLength_forIC = pr->seqLength;//min(pr->seqLength,1000);
    double m = -log(prReduced->q*prReduced->q+1)/2;
    double sigma = sqrt(log(prReduced->q*prReduced->q+1));
    for (int i=1;i<=prReduced->nbBranches;i++){
        std::poisson_distribution<int> distribution(nodesReduced[i]->B*seqLength_forIC);
        std::lognormal_distribution<double> distributionlogNormal(m,sigma);
        for (int j=0;j<pr->nbSampling;j++){
            double coef = (double)distributionlogNormal(generator);
            B_simul[j][i]=(double)distribution(generator)*coef/(seqLength_forIC);
        }
    }
    double maxD = nodesReduced[prReduced->nbINodes]->D; // the most recent tip date
    for (int i = prReduced->nbINodes+1; i <= prReduced->nbBranches; i++){
        if (nodesReduced[i]->D > maxD){
            maxD = nodesReduced[i]->D;
        }
    }
    double** T_simul = new double*[pr->nbSampling];
    double** H_simul = new double*[pr->nbSampling];
    double** HD_simul = new double*[pr->nbSampling];
    double* rho_simul = new double[pr->nbSampling];
    double* rho_lower = new double[pr->nbSampling];
    double* rho_upper = new double[pr->nbSampling];
    vector<double>* other_rhos_simul = new vector<double>[pr->nbSampling];
    std::poisson_distribution<int> distribution(br*pr->seqLength);
    double* tipDates = new double[prReduced->nbBranches - prReduced->nbINodes + 1];
    for (int i=0;i<prReduced->nbBranches - prReduced->nbINodes + 1;i++) {
        tipDates[i] = nodesReduced[i + prReduced->nbINodes]->D;
    }
    computeNewVariance(prReduced,nodesReduced);
    bool unique = true;
    for (int r=0;r<pr->nbSampling;r++){
        for (int j=0;j<=prReduced->nbBranches;j++){
            nodesReduced[j]->B=B_simul[r][j];
        }
        initialize_status(prReduced, nodesReduced);
        for (int g=1; g<=pr->ratePartition.size(); g++) {
            prReduced->multiplierRate[g] = pr->multiplierRate[g];
        }
        if (pr->ratePartition.size()>0){
            for (int i=0;i<=pr->ratePartition.size();i++) prReduced->multiplierRate[i] = 1;
        }
        if (pr->constraint) with_constraint_multirates(prReduced,nodesReduced,false);
        else without_constraint_multirates(prReduced,nodesReduced,false);
        rho_lower[r]  = prReduced->rhoLower;
        rho_upper[r] = prReduced->rhoUpper;
        if (prReduced->rhoLower != 0 && prReduced->rhoUpper!=0) unique  = false;
        T_simul[r]=new double[prReduced->nbBranches+1];
        for (int j=0;j<=prReduced->nbBranches;j++) T_simul[r][j]=nodesReduced[j]->D;
        rho_simul[r]=prReduced->rho;
        for (int g=1; g<=pr->ratePartition.size(); g++) {
            other_rhos_simul[r].push_back(prReduced->rho*prReduced->multiplierRate[g]);
        }
        calculate_tree_height(prReduced,nodesReduced);
        H_simul[r]=new double[prReduced->nbBranches+1];
        HD_simul[r]=new double[prReduced->nbBranches+1];
        for (int j=0;j<=prReduced->nbBranches;j++){
            H_simul[r][j]=nodesReduced[j]->H;
            HD_simul[r][j]=nodesReduced[j]->HD;
        }
        delete[] B_simul[r];
    }
    delete[] B_simul;
    delete prReduced;
    for (int i=0;i<nbC+pr->nbBranches-pr->nbINodes+1;i++){
        delete nodesReduced[i];
    }
    delete[] nodesReduced;
    if (unique){
        sort(rho_simul,pr->nbSampling);
        rho_left=rho_simul[int(0.025*pr->nbSampling)];
        rho_right=rho_simul[pr->nbSampling-int(0.025*pr->nbSampling)-1];
        if (pr->rho<rho_left) rho_left=pr->rho;
        if (pr->rho>rho_right) rho_right=pr->rho;
        double* T_sort = new double[pr->nbSampling];
        double* H_sort = new double[pr->nbSampling];
        double* HD_sort = new double[pr->nbSampling];
        vector<int> pre = preorder_polytomy_withTips(pr,nodes);
        for (vector<int>::iterator iter=pre.begin();iter!=pre.end();iter++){
            int i =  *iter;
            if (tab[i]!=-1) {
                for (int j=0;j<pr->nbSampling;j++) {
                    T_sort[j]=T_simul[j][tab[i]];
                    H_sort[j]=H_simul[j][tab[i]];
                    HD_sort[j]=HD_simul[j][tab[i]];
                }
                sort(T_sort,pr->nbSampling);
                sort(H_sort,pr->nbSampling);
                sort(HD_sort,pr->nbSampling);
                
                T_left[i]=T_sort[int(0.025*pr->nbSampling)];
                if (T_left[i]>nodes[i]->D) T_left[i]=nodes[i]->D;
                T_right[i]=T_sort[pr->nbSampling-int(0.025*pr->nbSampling)-1];
                if (T_right[i]<nodes[i]->D) T_right[i]=nodes[i]->D;
                
                H_left[i]=H_sort[int(0.025*pr->nbSampling)];
                if (H_left[i]>nodes[i]->H) H_left[i]=nodes[i]->H;
                H_right[i]=H_sort[pr->nbSampling-int(0.025*pr->nbSampling)-1];
                if (H_right[i]<nodes[i]->H) H_right[i]=nodes[i]->H;
                
                HD_left[i]=HD_sort[int(0.025*pr->nbSampling)];
                if (HD_left[i]>nodes[i]->HD) HD_left[i]=nodes[i]->HD;
                HD_right[i]=HD_sort[pr->nbSampling-int(0.025*pr->nbSampling)-1];
                if (HD_right[i]<nodes[i]->HD) HD_right[i]=nodes[i]->HD;
            } else {
                T_left[i] = T_left[nodes[i]->P];
                T_right[i] = T_right[nodes[i]->P];
                H_left[i] = H_left[nodes[i]->P];
                H_right[i] = H_right[nodes[i]->P];
                HD_left[i] = HD_left[nodes[i]->P];
                HD_right[i] = HD_right[nodes[i]->P];
            }
        }
        delete[] T_sort;
        delete[] H_sort;
        delete[] HD_sort;
        double* other_rhos_simul_sort = new double[pr->nbSampling];
        for (int g=1;g<=pr->ratePartition.size();g++){
            for (int r=0;r<pr->nbSampling;r++) {
                other_rhos_simul_sort[r]=other_rhos_simul[r][g-1];
            }
            sort(other_rhos_simul_sort,pr->nbSampling);
            other_rhos_left[g]=other_rhos_simul_sort[int(0.025*pr->nbSampling)];
            if (other_rhos_left[g]>pr->rho*pr->multiplierRate[g]) other_rhos_left[g]=pr->rho*pr->multiplierRate[g];
            other_rhos_right[g]=other_rhos_simul_sort[pr->nbSampling-int(0.025*pr->nbSampling)-1];
            if (other_rhos_right[g]<pr->rho*pr->multiplierRate[g]) other_rhos_right[g]=pr->rho*pr->multiplierRate[g];
        }
        delete[] other_rhos_simul_sort;
    } else {
        double* rho_simul_2 = new double[2*pr->nbSampling];
        for (int i=0;i<pr->nbSampling;i++){
            rho_simul_2[i] = rho_lower[i];
            rho_simul_2[2*pr->nbSampling-i-1] = rho_upper[i];
        }
        sort(rho_simul_2,2*pr->nbSampling);
        rho_left=rho_simul_2[int(0.025*pr->nbSampling*2)];
        rho_right=rho_simul_2[2*pr->nbSampling-int(0.025*pr->nbSampling*2)-1];
        if (pr->rho<rho_left) rho_left=pr->rho;
        if (pr->rho>rho_right) rho_right=pr->rho;
        double* T_sort = new double[2*pr->nbSampling];
        double* H_sort = new double[2*pr->nbSampling];
        double* HD_sort = new double[2*pr->nbSampling];
        vector<int> pre = preorder_polytomy_withTips(pr,nodes);
        for (vector<int>::iterator iter=pre.begin();iter!=pre.end();iter++){
            int i =  *iter;
            if (tab[i]!=-1) {
                for (int j=0;j<pr->nbSampling;j++) {
                    T_sort[j]=T_simul[j][tab[i]]*rho_simul[j] / rho_lower[j];
                    T_sort[2*pr->nbSampling-j-1]=T_simul[j][tab[i]]*rho_simul[j] / rho_upper[j];
                    H_sort[j]=H_simul[j][tab[i]]*rho_simul[j] / rho_lower[j];
                    H_sort[2*pr->nbSampling-j-1]=H_simul[j][tab[i]]*rho_simul[j] / rho_upper[j];
                    HD_sort[j]=HD_simul[j][tab[i]]*rho_simul[j] /  rho_lower[j];
                    HD_sort[2*pr->nbSampling-j-1]=HD_simul[j][tab[i]]*rho_simul[j] /  rho_upper[j];
                }
                sort(T_sort,2*pr->nbSampling);
                sort(H_sort,2*pr->nbSampling);
                sort(HD_sort,2*pr->nbSampling);
                
                T_left[i]=T_sort[int(0.025*pr->nbSampling*2)];
                if (T_left[i]>nodes[i]->D*pr->rho/pr->rhoUpper) T_left[i]=nodes[i]->D*pr->rho/pr->rhoUpper;
                T_right[i]=T_sort[2*pr->nbSampling-int(0.025*pr->nbSampling*2)-1];
                if (T_right[i]<nodes[i]->D*pr->rho/pr->rhoLower) T_right[i]=nodes[i]->D*pr->rho/pr->rhoLower;
                
                H_left[i]=H_sort[int(0.025*pr->nbSampling*2)];
                if (H_left[i]>nodes[i]->H*pr->rho/pr->rhoUpper) H_left[i]=nodes[i]->H*pr->rho/pr->rhoUpper;
                H_right[i]=H_sort[2*pr->nbSampling-int(0.025*pr->nbSampling*2)-1];
                if (H_right[i]<nodes[i]->H*pr->rho/pr->rhoLower) H_right[i]=nodes[i]->H*pr->rho/pr->rhoLower;
                
                HD_left[i]=HD_sort[int(0.025*pr->nbSampling*2)];
                if (HD_left[i]>nodes[i]->HD*pr->rho/pr->rhoUpper) HD_left[i]=nodes[i]->HD*pr->rho/pr->rhoUpper;
                HD_right[i]=HD_sort[2*pr->nbSampling-int(0.025*pr->nbSampling*2)-1];
                if (HD_right[i]<nodes[i]->HD*pr->rho/pr->rhoLower) HD_right[i]=nodes[i]->HD*pr->rho/pr->rhoLower;
                
            } else {
                T_left[i] = T_left[nodes[i]->P];
                T_right[i] = T_right[nodes[i]->P];
                H_left[i] = H_left[nodes[i]->P];
                H_right[i] = H_right[nodes[i]->P];
                HD_left[i] = HD_left[nodes[i]->P];
                HD_right[i] = HD_right[nodes[i]->P];
            }
        }
        delete[] T_sort;
        delete[] H_sort;
        delete[] HD_sort;
        delete[] rho_simul_2;
        double* other_rhos_simul_sort = new double[2*pr->nbSampling];
        for (int g=1;g<=pr->ratePartition.size();g++){
            for (int r=0;r<2*pr->nbSampling;r++) {
                other_rhos_simul_sort[r]=other_rhos_simul[r][g-1]*pr->rho/pr->rhoLower;
                other_rhos_simul_sort[2*pr->nbSampling-r-1]=other_rhos_simul[r][g-1]*pr->rho/pr->rhoUpper;
            }
            sort(other_rhos_simul_sort,2*pr->nbSampling);
            other_rhos_left[g]=other_rhos_simul_sort[int(0.025*pr->nbSampling*2)];
            if (other_rhos_left[g]>pr->rhoLower*pr->multiplierRate[g]) other_rhos_left[g]=pr->rhoLower*pr->multiplierRate[g];
            other_rhos_right[g]=other_rhos_simul_sort[2*pr->nbSampling-int(0.025*pr->nbSampling*2)-1];
            if (other_rhos_right[g]<pr->rhoUpper*pr->multiplierRate[g]) other_rhos_right[g]=pr->rhoUpper*pr->multiplierRate[g];
        }
        delete[] other_rhos_simul_sort;
    }
    delete[] tab;
    delete[] rho_lower;
    delete[] rho_upper;
    delete[] other_rhos_simul;
    for (int i=0;i<pr->nbSampling;i++){
        delete[] T_simul[i];
        delete[] H_simul[i];
        delete[] HD_simul[i];
    }
    delete[] T_simul;
    delete[] H_simul;
    delete[] HD_simul;
}

bool computeIC_bootstraps(InputOutputStream *io, Pr* pr,Node** nodes,double* &T_left,double* &T_right,double* &H_left,double* &H_right,double* &HD_left,double* &HD_right,double &rho_left,double& rho_right,double* &other_rhos_left,double* &other_rhos_right,int r, bool diffTopology){
    double** T_bootstrap = new double*[pr->nbBootstrap];
    double** H_bootstrap = new double*[pr->nbBootstrap];
    double** HD_bootstrap = new double*[pr->nbBootstrap];
    double* rho_bootstrap = new double[pr->nbBootstrap];
    double* rho_lower = new double[pr->nbBootstrap];
    double* rho_upper = new double[pr->nbBootstrap];
    double** other_rhos_bootstrap  = new double*[pr->nbBootstrap];
    int s = 0;
    if (io->inOutgroup){
        extrait_outgroup(io, pr, true);
    }
    bool constraintConsistent=true;
    if (diffTopology){
        cout<<"At least one bootstrap tree does not have the same topology as the input tree."<<endl;
        std::ostringstream oss;
        oss<<"- The bootstrap trees do not have the same topology as the input tree, so\nonly the confidence intervals of rates and root dates are calculated.\n";
        pr->warningMessage.push_back(oss.str());
    }
    streampos pos = (*(io->inBootstrapTree)).tellg();
    io->inBootstrapTree->seekg(0);
    Pr* pr_bootstrap = new Pr(pr->nbINodes,pr->nbBranches);
    pr_bootstrap->copy(pr);
    double median_rate = pr->rho_min;
    double br=0;
    bool unique = true;
    for (int y = 0; y< pr->nbBootstrap; y++){
        pr_bootstrap->internalConstraints.clear();
        pr_bootstrap->rooted = false;
        Node** nodes_bootstrap=tree2data(*(io->inBootstrapTree),pr_bootstrap,s);
        readInputDate(io, pr_bootstrap,nodes_bootstrap,constraintConsistent);
        computeSuc_polytomy(pr_bootstrap,nodes_bootstrap);
        if (pr->c == -1){
            pr_bootstrap->b = max(median_branch_lengths(pr_bootstrap,nodes_bootstrap),10./pr->seqLength);
        } else {
            pr_bootstrap->b = pr_bootstrap->c;
        }
        if (diffTopology){
            double minB = nodes_bootstrap[1]->B;
            if (pr->minblen < 0){
                for (int i=2; i <= pr_bootstrap->nbBranches; i++){
                    if (nodes_bootstrap[i]->B < minB) minB = nodes_bootstrap[i]->B;
                }
            }
            if (pr_bootstrap->rooted && pr_bootstrap->estimate_root!="" && pr_bootstrap->estimate_root!="k"){
                rooted2unrooted(pr_bootstrap,nodes_bootstrap);
            }
            collapseUnInformativeBranches(pr_bootstrap,nodes_bootstrap,false);
            if (!pr_bootstrap->rooted){
                unrooted2rooted(pr_bootstrap,nodes_bootstrap);
            }
            computeVariance(pr_bootstrap,nodes_bootstrap);
            if (pr->estimate_root=="" || pr->estimate_root=="k") constraintConsistent = initConstraint(pr_bootstrap, nodes_bootstrap);
            if (pr->e>0) calculateOutliers(pr_bootstrap,nodes_bootstrap,median_rate,false);
            if (pr->splitExternal) splitExternalBranches(pr_bootstrap,nodes_bootstrap);
            if (pr->constraint) imposeMinBlen(*(io->outResult),pr_bootstrap,nodes_bootstrap,minB,false);
        } else{
            if (!pr_bootstrap->rooted){
                unrooted2rooted(pr_bootstrap,nodes_bootstrap);
                reroot_rootedtree(br,r,pr_bootstrap,nodes_bootstrap);
            }
            if (pr->ratePartition.size()>0) assignRateGroupToTree(pr,nodes_bootstrap);
            computeVariance(pr_bootstrap,nodes_bootstrap);
            if (pr->splitExternal) splitExternalBranches(pr_bootstrap,nodes_bootstrap);
            initConstraint(pr_bootstrap, nodes_bootstrap);
            if (pr->e>0) remove_outlier_nodes(pr,nodes_bootstrap);
            if (pr->constraint){
                for (int i=0; i<= pr_bootstrap->nbBranches;i++){
                    nodes_bootstrap[i]->minblen = nodes[i]->minblen;
                }
            }
        }
        if (!pr->constraint){//LD without constraints
            if (pr->estimate_root==""){
                without_constraint_multirates(pr_bootstrap,nodes_bootstrap,true);
            }
            else if (!diffTopology){
                without_constraint_active_set_lambda_multirates(true,br,pr_bootstrap,nodes_bootstrap,true);
            } else {
                if (pr->estimate_root.compare("l")==0){
                    r=estimate_root_without_constraint_local_rooted(pr_bootstrap,nodes_bootstrap);
                } else{
                    r=estimate_root_without_constraint_rooted(pr_bootstrap,nodes_bootstrap);
                }
                if (r>0){
                    reroot_rootedtree(br,r,pr_bootstrap,nodes_bootstrap);
                    without_constraint_active_set_lambda_multirates(true,br,pr_bootstrap,nodes_bootstrap,true);
                }
            }
        } else {//QPD with temporal constrains
            if (pr->estimate_root==""){
                with_constraint_multirates(pr_bootstrap,nodes_bootstrap,true);
            }
            else  if (!diffTopology){
                with_constraint_active_set_lambda_multirates(true,br,pr_bootstrap,nodes_bootstrap,true);
            } else{
                if (pr->estimate_root.compare("l")==0){
                    r=estimate_root_with_constraint_local_rooted(pr_bootstrap,nodes_bootstrap);
                } else if (pr->estimate_root.compare("a")==0){
                    r=estimate_root_with_constraint_fast_rooted(pr_bootstrap,nodes_bootstrap);
                } else{
                    r=estimate_root_with_constraint_rooted(pr_bootstrap,nodes_bootstrap);
                }
                if (r>0){
                    reroot_rootedtree(br,r,pr_bootstrap,nodes_bootstrap);
                    with_constraint_active_set_lambda_multirates(true,br,pr_bootstrap,nodes_bootstrap,true);
                }
            }
        }
        if (pr->rhoLower != 0 && pr->rhoUpper!=0) unique  = false;
        if (pr->verbose){
            cout<<"Tree "<<y+1<<" rate: "<<pr_bootstrap->rho<<", tMRCA: "<<nodes_bootstrap[0]->D<<endl;
        }
        if (!diffTopology){
            T_bootstrap[y] = new double[pr_bootstrap->nbBranches+1];
            H_bootstrap[y] = new double[pr_bootstrap->nbBranches+1];
            HD_bootstrap[y] = new double[pr_bootstrap->nbBranches+1];
            calculate_tree_height(pr_bootstrap,nodes_bootstrap);
            for (int i=0;i<=pr_bootstrap->nbBranches;i++){
                T_bootstrap[y][i] = nodes_bootstrap[i]->D;
                H_bootstrap[y][i] = nodes_bootstrap[i]->H;
                HD_bootstrap[y][i] = nodes_bootstrap[i]->HD;
            }
        } else{
            T_bootstrap[y] = new double[1];
            T_bootstrap[y][0] = nodes_bootstrap[0]->D;
        }
        rho_bootstrap[y] = pr_bootstrap->rho;
        rho_lower[r]  = pr->rhoLower;
        rho_upper[r] = pr->rhoUpper;
        other_rhos_bootstrap[y] = new double[pr->ratePartition.size()];
        for (int g=1; g<=pr->ratePartition.size(); g++) {
            other_rhos_bootstrap[y][g-1] = pr_bootstrap->rho*pr_bootstrap->multiplierRate[g];
        }
        for (int i=0;i<=pr_bootstrap->nbBranches;i++) delete nodes_bootstrap[i];
        delete[] nodes_bootstrap;
    }
    (*(io->inBootstrapTree)).clear();
    (*(io->inBootstrapTree)).seekg(pos);
    if (unique){
        if (!diffTopology){
            double* T_sort = new double[pr->nbBootstrap];
            double* H_sort = new double[pr->nbBootstrap];
            double* HD_sort = new double[pr->nbBootstrap];
            for (int i=0;i<=pr->nbBranches;i++){
                for (int j=0;j<pr->nbBootstrap;j++) {
                    T_sort[j]=T_bootstrap[j][i];
                    H_sort[j]=H_bootstrap[j][i];
                    HD_sort[j]=HD_bootstrap[j][i];
                }
                sort(T_sort,pr->nbBootstrap);
                sort(H_sort,pr->nbBootstrap);
                sort(HD_sort,pr->nbBootstrap);
                T_left[i]=T_sort[int(0.025*pr->nbBootstrap)];
                if (T_left[i]>nodes[i]->D) T_left[i]=nodes[i]->D;
                T_right[i]=T_sort[pr->nbBootstrap-int(0.025*pr->nbBootstrap)-1];
                if (T_right[i]<nodes[i]->D) T_right[i]=nodes[i]->D;
                
                H_left[i]=H_sort[int(0.025*pr->nbBootstrap)];
                if (H_left[i]>nodes[i]->H) H_left[i]=nodes[i]->H;
                H_right[i]=H_sort[pr->nbBootstrap-int(0.025*pr->nbBootstrap)-1];
                if (H_right[i]<nodes[i]->H) H_right[i]=nodes[i]->H;
                
                HD_left[i]=HD_sort[int(0.025*pr->nbBootstrap)];
                if (HD_left[i]>nodes[i]->HD) HD_left[i]=nodes[i]->HD;
                HD_right[i]=HD_sort[pr->nbBootstrap-int(0.025*pr->nbBootstrap)-1];
                if (HD_right[i]<nodes[i]->HD) HD_right[i]=nodes[i]->HD;
            }
            delete[] T_sort;
            delete[] H_sort;
            delete[] HD_sort;
        } else{
            double* T_sort = new double[pr->nbBootstrap];
            for (int i=0; i<pr->nbBootstrap; i++){
                T_sort[i] = T_bootstrap[i][0];
            }
            sort(T_sort,pr->nbBootstrap);
            T_left[0]=T_sort[int(0.025*pr->nbBootstrap)];
            if (T_left[0]>nodes[0]->D) T_left[0]=nodes[0]->D;
            T_right[0]=T_sort[pr->nbBootstrap-int(0.025*pr->nbBootstrap)-1];
            if (T_right[0]<nodes[0]->D) T_right[0]=nodes[0]->D;
            delete[] T_sort;
        }
        sort(rho_bootstrap,pr->nbBootstrap);
        rho_left=rho_bootstrap[int(0.025*pr->nbBootstrap)];
        rho_right=rho_bootstrap[pr->nbSampling-int(0.025*pr->nbBootstrap)-1];
        if (pr->rho<rho_left) rho_left=pr->rho;
        if (pr->rho>rho_right) rho_right=pr->rho;
        double* other_rhos_bootstrap_sort = new double[pr->nbBootstrap];
        for (int g=1;g<=pr->ratePartition.size();g++){
            for (int r=0;r<pr->nbBootstrap;r++) {
                other_rhos_bootstrap_sort[r]=other_rhos_bootstrap[r][g-1];
            }
            sort(other_rhos_bootstrap_sort,pr->nbBootstrap);
            other_rhos_left[g]=other_rhos_bootstrap_sort[int(0.025*pr->nbBootstrap)];
            if (other_rhos_left[g]>pr->rho*pr->multiplierRate[g]) other_rhos_left[g]=pr->rho*pr->multiplierRate[g];
            other_rhos_right[g]=other_rhos_bootstrap_sort[pr->nbBootstrap-int(0.025*pr->nbSampling)-1];
            if (other_rhos_right[g]<pr->rho*pr->multiplierRate[g]) other_rhos_right[g]=pr->rho*pr->multiplierRate[g];
        }
    } else {
        if (!diffTopology){
            double* T_sort = new double[2*pr->nbBootstrap];
            double* H_sort = new double[2*pr->nbBootstrap];
            double* HD_sort = new double[2*pr->nbBootstrap];
            for (int i=0;i<=pr->nbBranches;i++){
                for (int j=0;j<pr->nbBootstrap;j++) {
                    T_sort[j]=T_bootstrap[j][i]*rho_bootstrap[j] / rho_lower[j];
                    T_sort[2*pr->nbBootstrap-j-1]=T_bootstrap[j][i]*rho_bootstrap[j] / rho_upper[j];
                    H_sort[j]=H_bootstrap[j][i]*rho_bootstrap[j] / rho_lower[j];
                    H_sort[2*pr->nbBootstrap-j-1]=H_bootstrap[j][i]*rho_bootstrap[j] / rho_upper[j];
                    HD_sort[j]=HD_bootstrap[j][i]*rho_bootstrap[j] / rho_lower[j];
                    HD_sort[2*pr->nbBootstrap-j-1]=HD_bootstrap[j][i]*rho_bootstrap[j] / rho_upper[j];
                }
                sort(T_sort,2*pr->nbBootstrap);
                sort(H_sort,2*pr->nbBootstrap);
                sort(HD_sort,2*pr->nbBootstrap);
    
                T_left[i]=T_sort[int(0.025*pr->nbBootstrap*2)];
                if (T_left[i]>nodes[i]->D) T_left[i]=nodes[i]->D;
                T_right[i]=T_sort[2*pr->nbBootstrap-int(0.025*pr->nbBootstrap*2)-1];
                if (T_right[i]<nodes[i]->D) T_right[i]=nodes[i]->D;
                
                H_left[i]=H_sort[int(0.025*pr->nbBootstrap*2)];
                if (H_left[i]>nodes[i]->H) H_left[i]=nodes[i]->H;
                H_right[i]=H_sort[2*pr->nbBootstrap-int(0.025*pr->nbBootstrap*2)-1];
                if (H_right[i]<nodes[i]->H) H_right[i]=nodes[i]->H;
                
                HD_left[i]=HD_sort[int(0.025*pr->nbBootstrap*2)];
                if (HD_left[i]>nodes[i]->HD) HD_left[i]=nodes[i]->HD;
                HD_right[i]=HD_sort[2*pr->nbBootstrap-int(0.025*pr->nbBootstrap*2)-1];
                if (HD_right[i]<nodes[i]->HD) HD_right[i]=nodes[i]->HD;
            }
            delete[] T_sort;
            delete[] H_sort;
            delete[] HD_sort;
        } else{
            double* T_sort = new double[pr->nbBootstrap*2];
            for (int i=0; i<pr->nbBootstrap; i++){
                T_sort[i] = T_bootstrap[i][0]*rho_bootstrap[i] / rho_lower[i];
                T_sort[2*pr->nbBootstrap-i-1] = T_bootstrap[i][0]*rho_bootstrap[i] / rho_upper[i];
            }
            sort(T_sort,2*pr->nbBootstrap);
            T_left[0]=T_sort[int(0.025*pr->nbBootstrap*2)];
            if (T_left[0]>nodes[0]->D) T_left[0]=nodes[0]->D;
            T_right[0]=T_sort[2*pr->nbBootstrap-int(0.025*pr->nbBootstrap*2)-1];
            if (T_right[0]<nodes[0]->D) T_right[0]=nodes[0]->D;
            delete[] T_sort;
        }
        double* rho_sort = new double[2*pr->nbBootstrap];
        for (int i=0; i<pr->nbBootstrap; i++){
            rho_sort[i] = rho_lower[i];
            rho_sort[2*pr->nbBootstrap-i-1] = rho_upper[i];
        }
        sort(rho_sort,2*pr->nbBootstrap);
        rho_left=rho_sort[int(0.025*pr->nbBootstrap*2)];
        rho_right=rho_sort[2*pr->nbSampling-int(0.025*pr->nbBootstrap*2)-1];
        if (pr->rho<rho_left) rho_left=pr->rho;
        if (pr->rho>rho_right) rho_right=pr->rho;
        double* other_rhos_bootstrap_sort = new double[2*pr->nbBootstrap];
        for (int g=1;g<=pr->ratePartition.size();g++){
            for (int r=0;r<pr->nbBootstrap;r++) {
                other_rhos_bootstrap_sort[r]=other_rhos_bootstrap[r][g-1]*rho_bootstrap[r] / rho_lower[r];
                other_rhos_bootstrap_sort[2*pr->nbBootstrap-r-1]=other_rhos_bootstrap[r][g-1]*rho_bootstrap[r] / rho_upper[r];
            }
            sort(other_rhos_bootstrap_sort,2*pr->nbBootstrap);
            other_rhos_left[g]=other_rhos_bootstrap_sort[int(0.025*pr->nbBootstrap*2)];
            if (other_rhos_left[g]>pr->rhoLower*pr->multiplierRate[g]) other_rhos_left[g]=pr->rhoLower*pr->multiplierRate[g];
            other_rhos_right[g]=other_rhos_bootstrap_sort[2*pr->nbBootstrap-int(0.025*pr->nbSampling*2)-1];
            if (other_rhos_right[g]<pr->rhoUpper*pr->multiplierRate[g]) other_rhos_right[g]=pr->rhoUpper*pr->multiplierRate[g];
        }
        delete[] rho_sort;
        delete[] other_rhos_bootstrap_sort;
    }
    delete[] rho_bootstrap;
    delete[] other_rhos_bootstrap;
    if (!diffTopology){
        for (int i=0;i<pr->nbBootstrap;i++){
            delete[] H_bootstrap[i];
            delete[] HD_bootstrap[i];
        }
    }
    for (int i=0;i<pr->nbBootstrap;i++){
        delete[] T_bootstrap[i];
    }
    delete[] T_bootstrap;
    delete[] H_bootstrap;
    delete[] HD_bootstrap;
    delete pr_bootstrap;
    return diffTopology;
}

void output(double br,int y, InputOutputStream *io, Pr* pr,Node** nodes,ostream& f,ostream& tree2,ostream& tree3,int r, bool diffTopology){
    if (!pr->constraint && pr->ci) {
        std::ostringstream oss;
        oss<<"- Confidence intervals are not warranted under non-constraint mode.\n";
        pr->warningMessage.push_back(oss.str());
    }
    if (pr->rhoLower != pr->rhoUpper){
        std::ostringstream oss;
        double tip_date = 0;
        for (int i=pr->nbINodes;i<=pr->nbBranches;i++){
            if (nodes[i]->type=='p'){
                tip_date=nodes[i]->D;
                break;
            }
        }
        oss<<"- There's not enough input date information to have a unique solution. \n";
        if  (pr->rhoUpper == 0){
            oss<<"The reported results show the lower solution: lower_rate and lower_dates.\n";
            oss<<"For any coefficient k>=1, there's a corresponding solution where the rate is lower_rate*k,\n";
            oss<<"and the date of any node i is lower_dates[i]/k";
            if (tip_date!=0) oss<<"+"<<tip_date<<"*(k-1)/k\n";
            else oss<<"\n";
        } else if (pr->rhoLower == 0) {
            oss<<"The reported results show the upper solution: upper_rate and upper_dates.\n";
            oss<<"For any coefficient k<=1, there's a corresponding solution where the rate is upper_rate*k,\n";
            oss<<"and the date of any node i is upper_dates[i]/k";
            if (tip_date!=0) oss<<"+"<<tip_date<<"*(k-1)/k\n";
            else oss<<"\n";
        } else {
            oss<<"The reported results show the lower and upper solutions: lower_rate:upper_rate,\n";
            oss<<"lower_dates:upper_dates. The output date tree reports the average date from the 2 \n";
            oss<<"boundaries for any node. ";
            if (!pr->ci) oss<<"And the 2 boundaries are writen as confidence interval.\n";
            else oss<<"\n";
            oss<<"For any coefficient k such that 1 <= k <= upper_rate/lower_rate, there is a \n";
            oss<<"corresponding solution where the rate is lower_rate*k and the date of any node i is\n";
            oss<<"lower_dates[i]/k";
            if (tip_date!=0) oss<<"+"<<tip_date<<"*(k-1)/k\n";
            else oss<<"\n";
        }
        pr->warningMessage.push_back(oss.str());
    }
    ostringstream tMRCA,tMRCALower, tMRCAUpper;
    if (pr->outDateFormat==2){
        tMRCA<<realToYearMonthDay(nodes[0]->D);
        if (pr->rhoLower != 0) tMRCALower<<realToYearMonthDay(nodes[0]->D*pr->rho/pr->rhoLower);
        if (pr->rhoUpper != 0) tMRCAUpper<<realToYearMonthDay(nodes[0]->D*pr->rho/pr->rhoUpper);
    } else if (pr->outDateFormat==3){
        tMRCA<<realToYearMonth(nodes[0]->D);
        if (pr->rhoLower != 0) tMRCALower<<realToYearMonth(nodes[0]->D*pr->rho/pr->rhoLower);
        if (pr->rhoUpper != 0) tMRCAUpper<<realToYearMonth(nodes[0]->D*pr->rho/pr->rhoUpper);
    } else {
        tMRCA<<nodes[0]->D;
        if (pr->rhoLower != 0) tMRCALower<<nodes[0]->D*pr->rho/pr->rhoLower;
        if (pr->rhoUpper != 0) tMRCAUpper<<nodes[0]->D*pr->rho/pr->rhoUpper;
    }
    std::ostringstream oss;
    oss<<"- Dating results:\n";
    if (pr->ratePartition.size()==0) {
        if (pr->rhoLower  == pr->rhoUpper) oss<<" rate "<<pr->rho<<", ";
        else if (pr->rhoLower ==  0) oss<<" rate NA:"<<pr->rhoUpper<<", ";
        else if (pr->rhoUpper == 0) oss<<" rate "<<pr->rhoLower<<":NA, ";
        else oss<<" rate "<<pr->rhoLower<<":"<<pr->rhoUpper<<", ";
    } else {
        if (pr->multiplierRate[0]!=-1){
            if (pr->splitExternal) oss<<" rate internal Branches "<<pr->rho<<", ";
            else oss<<" rate "<<pr->rho<<", ";
        }
        for (int i=1; i<=pr->ratePartition.size(); i++) {
            if (pr->multiplierRate[i]>0)
                oss<<"rate "<<pr->ratePartition[i-1]->name.c_str()<<" "<<pr->rho*pr->multiplierRate[i]<<", ";
        }
    }
    if (pr->rhoLower  == pr->rhoUpper) oss<<"tMRCA "<<tMRCA.str()<<" , objective function "<<pr->objective<<"\n";
    else if (pr->rhoLower ==  0) oss<<"tMRCA NA:"<<tMRCAUpper.str()<<" , objective function "<<pr->objective<<"\n";
    else if (pr->rhoUpper == 0) oss<<"tMRCA "<<tMRCALower.str()<<":NA , objective function "<<pr->objective<<"\n";
    else oss<<"tMRCA "<<tMRCALower.str()<<":"<<tMRCAUpper.str()<<" , objective function "<<pr->objective<<"\n";
    pr->resultMessage.push_back(oss.str());
    //variance = 2
    if (pr->variance==2){
        cout<<"Re-estimating using variances based on the branch lengths of the first run ...\n";
        if (pr->estimate_root=="") {
            computeNewVariance(pr,nodes);
            if (pr->constraint) {
                with_constraint_multirates(pr,nodes,false);
            }
            else{
                without_constraint_multirates(pr,nodes,false);
            }
        }
        else{
            computeNewVarianceEstimateRoot(pr,nodes);
            if (pr->constraint){
                with_constraint_active_set_lambda_multirates(true,br,pr, nodes,true);
            }
            else{
                without_constraint_active_set_lambda_multirates(true,br,pr, nodes,true);
            }
        }
        std::ostringstream oss;
        oss<<"- Results of the second run with variances based on results of the first run:\n";
        ostringstream tMRCA,tMRCALower, tMRCAUpper;
        if (pr->outDateFormat==2){
            tMRCA<<realToYearMonthDay(nodes[0]->D);
            if (pr->rhoLower != 0) tMRCALower<<realToYearMonthDay(nodes[0]->D*pr->rho/pr->rhoLower);
            if (pr->rhoUpper != 0) tMRCAUpper<<realToYearMonthDay(nodes[0]->D*pr->rho/pr->rhoUpper);
        } else if (pr->outDateFormat==3){
            tMRCA<<realToYearMonth(nodes[0]->D);
            if (pr->rhoLower != 0) tMRCALower<<realToYearMonth(nodes[0]->D*pr->rho/pr->rhoLower);
            if (pr->rhoUpper != 0) tMRCAUpper<<realToYearMonth(nodes[0]->D*pr->rho/pr->rhoUpper);
        } else {
            tMRCA<<nodes[0]->D;
            if (pr->rhoLower != 0) tMRCALower<<nodes[0]->D*pr->rho/pr->rhoLower;
            if (pr->rhoUpper != 0) tMRCAUpper<<nodes[0]->D*pr->rho/pr->rhoUpper;
        }
        if (pr->ratePartition.size()==0) {
            if (pr->rhoLower  == pr->rhoUpper) oss<<" rate "<<pr->rho<<", ";
            else if (pr->rhoLower ==  0) oss<<" rate NA:"<<pr->rhoUpper<<", ";
            else if (pr->rhoUpper == 0) oss<<" rate "<<pr->rhoLower<<":NA, ";
            else oss<<" rate "<<pr->rhoLower<<":"<<pr->rhoUpper<<", ";
        } else {
            if (pr->multiplierRate[0]!=-1){
                if (pr->splitExternal) oss<<" rate internal Branches "<<pr->rho<<", ";
                else oss<<" rate "<<pr->rho<<", ";
            }
            for (int i=1; i<=pr->ratePartition.size(); i++) {
                if (pr->multiplierRate[i]>0)
                    oss<<"rate "<<pr->ratePartition[i-1]->name.c_str()<<" "<<pr->rho*pr->multiplierRate[i]<<", ";
            }
        }
        if (pr->rhoLower  == pr->rhoUpper) oss<<"tMRCA "<<tMRCA.str()<<" , objective function "<<pr->objective<<"\n";
        else if (pr->rhoLower ==  0) oss<<"tMRCA NA:"<<tMRCAUpper.str()<<" , objective function "<<pr->objective<<"\n";
        else if (pr->rhoUpper == 0) oss<<"tMRCA "<<tMRCALower.str()<<":NA , objective function "<<pr->objective<<"\n";
        else oss<<"tMRCA "<<tMRCALower.str()<<":"<<tMRCAUpper.str()<<" , objective function "<<pr->objective<<"\n";
        pr->resultMessage.push_back(oss.str());
    }
    //Write input trees out after all pre-process steps (collapse branches, remove outgroups, rooting, remove outliers etc)
    if (pr->verbose){
        std::ostringstream oss;
        oss<<"- Evolution distance & date of each node:\n";
        oss<<"Node  Branch_length Branch_elapsed_time   RttDistance NodeDate\n";
        double* r2t = rtt(pr,nodes);
        for (int i=0;i<=pr->nbBranches;i++){
            if (i<pr->nbINodes) oss<<"node_"<<i<<"   ";
            else oss<<nodes[i]->L<<"    ";
            if (i==0) oss<<"NA  NA  ";
            else oss<<nodes[i]->B<<" "<<pr->rho*(nodes[i]->D-nodes[nodes[i]->P]->D)<<"   ";
            oss<<r2t[i]<<"    "<<nodes[i]->D<<"\n";
        }
        pr->resultMessage.push_back(oss.str());
    }
    //Replace branche length by time-scaled lengths
    if (pr->ratePartition.size()==0){
        for (int i=1;i<=pr->nbBranches;i++){
            nodes[i]->B=pr->rho*(nodes[i]->D-nodes[nodes[i]->P]->D);
        }
    }
    else {
        for (int i=1;i<=pr->nbBranches;i++){
            int g = nodes[i]->rateGroup;
            if (g==0) {
                nodes[i]->B=pr->rho*(nodes[i]->D-nodes[nodes[i]->P]->D);
            }
            else{
                nodes[i]->B=pr->rho*pr->multiplierRate[g]*(nodes[i]->D-nodes[nodes[i]->P]->D);
            }
        }
    }
    //Warning message without option -c
    
    if (!pr->constraint){
        int count=0;
        for (int i=1;i<=pr->nbBranches;i++){
            if (nodes[i]->D<nodes[nodes[i]->P]->D) count++;
        }
        
        if (count>0) {
            std::ostringstream oss;
            oss<<"- Number of violated temporal constraints (nodes having date smaller than the one of its parent):  "<<count<<" ("<<(count*100)/(double)pr->nbBranches<<"%%). Try without option -F to impose temporal constraints on the estimated trees.\n";
            pr->warningMessage.push_back(oss.str());
        }
        
    }
    calculate_tree_height(pr,nodes);
    //Calculate confidence intervals
    if (!pr->ci){
        tree2<<"tree "<<y<<" = ";
        if (pr->rhoLower != pr->rhoUpper && pr->rhoLower !=0 && pr->rhoUpper !=0){
            double* T_min=new double[pr->nbBranches+1];
            double* T_max=new double[pr->nbBranches+1];
            double* HD_min=new double[pr->nbBranches+1];
            double* HD_max=new double[pr->nbBranches+1];
            for (int j=0;j<=pr->nbBranches;j++){
                T_min[j] = nodes[j]->D * pr->rho / pr->rhoLower;
                HD_min[j] = nodes[j]->HD * pr->rho / pr->rhoLower;
                T_max[j] = nodes[j]->D * pr->rho / pr->rhoUpper;
                HD_max[j] = nodes[j]->HD * pr->rho / pr->rhoUpper;
                nodes[j]->D = (T_min[j]+T_max[j])/2;
            }
            tree2<<nexusICDate(0,pr,nodes,T_min,T_max,HD_min,HD_max).c_str();
        } else {
            tree2<<nexusDate(0,pr,nodes).c_str();
        }
        int n=0;
        tree3<<newick(0,0,pr,nodes,n).c_str();
    }
    else if (pr->haveUnique || (pr->haveLower && pr->haveUpper)){
        double* T_min = new double[pr->nbBranches+1];
        double* T_max = new double[pr->nbBranches+1];
        double* H_min = new double[pr->nbBranches+1];
        double* H_max = new double[pr->nbBranches+1];
        double* HD_min = new double[pr->nbBranches+1];
        double* HD_max = new double[pr->nbBranches+1];
        double rho_left,rho_right;
        double* other_rhos_left = new double[pr->ratePartition.size()+1];
        double* other_rhos_right = new double[pr->ratePartition.size()+1];
        if (pr->bootstraps_file==""){
            cout<<"Computing confidence intervals using sequence length "<<pr->seqLength<<" and a lognormal\n relaxed clock with mean 1, standard deviation "<<pr->q<<" (settable via option -q)"<<endl;
            computeIC(br,pr,nodes,T_min,T_max,H_min,H_max,HD_min,HD_max,rho_left,rho_right,other_rhos_left,other_rhos_right);
        } else {
            cout<<"Computing confidence intervals using input bootstrap trees ..."<<endl;
            computeIC_bootstraps(io,pr,nodes,T_min,T_max,H_min,H_max,HD_min,HD_max,rho_left,rho_right,other_rhos_left,other_rhos_right,r,diffTopology);
        }
        std::ostringstream oss;
        oss<<"- Results with confidence intervals:\n";
        ostringstream tMRCA,tMRCALower, tMRCAUpper,tmin,tmax;
        if (pr->outDateFormat==2){
            tMRCA<<realToYearMonthDay(nodes[0]->D);
            if (pr->rhoLower != 0) tMRCALower<<realToYearMonthDay(nodes[0]->D*pr->rho/pr->rhoLower);
            if (pr->rhoUpper != 0) tMRCAUpper<<realToYearMonthDay(nodes[0]->D*pr->rho/pr->rhoUpper);
            tmin<<realToYearMonthDay(T_min[0]);
            tmax<<realToYearMonthDay(T_max[0]);
        } else if (pr->outDateFormat==3){
            tMRCA<<realToYearMonth(nodes[0]->D);
            if (pr->rhoLower != 0) tMRCALower<<realToYearMonth(nodes[0]->D*pr->rho/pr->rhoLower);
            if (pr->rhoUpper != 0) tMRCAUpper<<realToYearMonth(nodes[0]->D*pr->rho/pr->rhoUpper);
            tmin<<realToYearMonth(T_min[0]);
            tmax<<realToYearMonth(T_max[0]);
        } else {
            tMRCA<<nodes[0]->D;
            if (pr->rhoLower != 0) tMRCALower<<nodes[0]->D*pr->rho/pr->rhoLower;
            if (pr->rhoUpper != 0) tMRCAUpper<<nodes[0]->D*pr->rho/pr->rhoUpper;
            tmin<<T_min[0];
            tmax<<T_max[0];
        }
        if (pr->ratePartition.size()==0 || !pr->splitExternal) oss<<" rate ";
        else oss<<" rate internal Branches ";
        if (pr->rhoLower  == pr->rhoUpper){
            oss<<pr->rho<<" ["<<rho_left<<"; "<<rho_right<<"], ";
            for (int i=1; i<=pr->ratePartition.size(); i++) {
                if (pr->multiplierRate[i]>0) oss<<"rate "<<pr->ratePartition[i-1]->name.c_str()<<" "<<pr->rho*pr->multiplierRate[i]<<" ["<<other_rhos_left[i]<<"; "<<other_rhos_right[i]<<"], ";
            }
        } else{
            oss<<pr->rhoLower<<":"<<pr->rhoUpper<<" ["<<rho_left<<"; "<<rho_right<<"], ";
            for (int i=1; i<=pr->ratePartition.size(); i++) {
                if (pr->multiplierRate[i]>0) oss<<"rate "<<pr->ratePartition[i-1]->name.c_str()<<" "<<pr->rhoLower*pr->multiplierRate[i]<<":"<<pr->rhoUpper*pr->multiplierRate[i]<<" ["<<other_rhos_left[i]<<"; "<<other_rhos_right[i]<<"], ";
            }
        }
        if (pr->rhoLower  == pr->rhoUpper){
            oss<<"tMRCA "<<tMRCA.str()<<" ["<<tmin.str()<<"; "<<tmax.str()<<"], objective function "<<pr->objective<<"\n";
        } else{
            oss<<"tMRCA "<<tMRCALower.str()<<":"<<tMRCAUpper.str()<<" ["<<tmin.str()<<"; "<<tmax.str()<<"], objective function "<<pr->objective<<"\n";
        }
        pr->resultMessage.push_back(oss.str());
        
        delete[] other_rhos_left;
        delete[] other_rhos_right;
        
        tree2<<"tree "<<y<<" = ";
        if (!diffTopology){
            if (pr->rhoLower!=0 && pr->rhoUpper!=0){
                for (int j=0;j<=pr->nbBranches;j++){
                    nodes[j]->D = (nodes[j]->D * pr->rho)*(1 / pr->rhoLower + 1/ pr->rhoUpper)/2;
                }    
            }
            tree2<<nexusICDate(0,pr,nodes,T_min,T_max,HD_min,HD_max).c_str();
        } else{
            if (pr->rhoLower!=0 && pr->rhoUpper!=0){
                nodes[0]->D = (nodes[0]->D * pr->rho)*(1 / pr->rhoLower + 1 / pr->rhoUpper)/2;
            }
            tree2<<nexusDate(0,pr,nodes).c_str();
        }
        int n=0;
        tree3<<newick(0,0,pr,nodes,n).c_str();
        
        delete[] T_min;
        delete[] T_max;
        delete[] H_min;
        delete[] H_max;
    }
    
    if (pr->rho==pr->rho_min) {
        std::ostringstream oss;
        oss<<"- The estimated rate reaches the given lower bound. To change the lower bound, use option -t.\n";
        pr->warningMessage.push_back(oss.str());
    }
    if ((pr->warningMessage).size()>0){
        f<<"\n*WARNINGS:\n";
        cout<<"\n*WARNINGS:\n";
        for (int i=0;i<pr->warningMessage.size();i++){
            f<<string(pr->warningMessage[i]).c_str();
            cout<<string(pr->warningMessage[i]).c_str();
        }
        
    }
    f<<"\n*RESULTS:\n";
    cout<<"\n*RESULTS:\n";
    for (int i=0;i<pr->resultMessage.size();i++){
        f<<string(pr->resultMessage[i]).c_str();
        cout<<string(pr->resultMessage[i]).c_str();
    }
}
