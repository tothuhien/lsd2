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

#include "pr.h"
#include "options.h"
#include "readData.h"
#include "dating.h"
#include "utils.h"
#include "node.h"
#include "date.h"
#include "estimate_root.h"
#include "confidence_interval.h"
#include "outliers.h"
#include "lsd.h"
#include <ctime>


using namespace std;
using namespace lsd;

int lsd::buildTimeTree( int argc, char** argv, InputOutputStream *inputOutput)
{
    Pr* opt = getOptions( argc, argv);
    // initialise input and output streams
    InputOutputStream *io;
    if (inputOutput)
        io = inputOutput;
    else
        io = new InputOutputFile(opt);
    printInterface(*(io->outResult), opt);
    clock_t start = clock();
    double elapsed_time;
    if (io->inOutgroup){
        extrait_outgroup(io, opt, false);
    }
    *(io->outTree2)<<"#NEXUS\n";
    bool constraintConsistent=true;
    int s=0;
    double median_rate = opt->rho_min;
    if (opt->partitionFile!="") {
        readPartitionFile(*(io->inPartition), opt);
        for (int i=0;i<=opt->ratePartition.size();i++){
            opt->givenRate.push_back(false);
        }
    }
    for (int y=1;y<=opt->nbData;y++){
        *(io->outResult)<<"\nTree "<<y<<" \n";
        cout<<"\nTREE "<<y<<endl;
        cout<<"*PROCESSING:"<<endl;
        cout<<"Reading the tree ... "<<endl;
        opt->init();
        Node** nodes=tree2data(*(io->inTree),opt,s);
        readInputDate(io, opt,nodes,constraintConsistent);
        if (!constraintConsistent){
            ostringstream oss;
            oss<<"- There's conflict in the input temporal constraints.\n";
            opt->warningMessage.push_back(oss.str());
        }
        computeSuc_polytomy(opt,nodes);
        double minB = nodes[(!opt->rooted)+1]->B;
        if (opt->minblen < 0){
            for (int i=2; i <= opt->nbBranches; i++){
                if (nodes[i]->B < minB) minB = nodes[i]->B;
            }
        }
        bool diffTopology = false;
        if (io->inBootstrapTree){
            streampos pos = (*(io->inBootstrapTree)).tellg();
            io->inBootstrapTree->seekg(0);
            int lineNb=getLineNumber(*(io->inBootstrapTree));
            Pr* pr_bootstrap = new Pr(opt->nbINodes,opt->nbBranches);
            pr_bootstrap->copy(opt);
            opt->nbBootstrap = lineNb;
            for (int k = 0; k< opt->nbBootstrap; k++){
                int sk;
                pr_bootstrap->rooted = false;
                Node** nodes_bootstrap=tree2data(*(io->inBootstrapTree),pr_bootstrap,sk);
                computeSuc_polytomy(pr_bootstrap,nodes_bootstrap);
                if ((pr_bootstrap->nbINodes != opt->nbINodes) || (pr_bootstrap->nbBranches != opt->nbBranches) || !checkTopology(opt,nodes,nodes_bootstrap)){
                    diffTopology = true;
                    break;
                }
            }
            (*(io->inBootstrapTree)).clear();
            (*(io->inBootstrapTree)).seekg(pos);
        }
        if (diffTopology || !(io->inBootstrapTree)){
            if (opt->rooted && opt->estimate_root!="" && opt->estimate_root!="k"){
                rooted2unrooted(opt,nodes);
            }
            collapseUnInformativeBranches(opt,nodes,true);
        }
        if (!opt->rooted){
            unrooted2rooted(opt,nodes);
        }
        if (y==1){
            *(io->outTree2)<<"Begin trees;\n";
        }
        if (opt->c == -1){
            opt->b = max(median_branch_lengths(opt,nodes),10./opt->seqLength);
            if (opt->variance>0){
                cout<<"Parameter to adjust variances was set to "<<opt->b<<" (settable via option -b)"<<endl;
                *(io->outResult)<<"Parameter to adjust variances was set to "<<opt->b<<" (settable via option -b)\n";
            }
        } else {
            opt->b = opt->c;
        }
        computeVariance(opt,nodes);
        double br=0;
        int r=0;
        if (opt->givenRate[0]){
            string line;
            if( getline(*(io->inRate), line)) {
                vector<double> all_rates = read_double_from_line(line);
                if (all_rates.size()==0){
                    cerr<<"Error from the input rate file, expect a real number"<<endl;
                    exit(EXIT_FAILURE);
                }
                opt->rho = all_rates[0];
                if (all_rates.size() > 1){
                    int i = 1;
                    while (i<=opt->ratePartition.size() && i<all_rates.size()){
                        opt->multiplierRate.push_back(all_rates[i]/opt->rho);
                        opt->givenRate[i] = true;
                        i++;
                    }
                }
            } else{
                std::ostringstream oss;
                oss<<"- There are less given rates than the number of given trees.\n";
                opt->warningMessage.push_back(oss.str());
                opt->givenRate[0] = false;
            }
        }
        if (opt->estimate_root=="" || opt->estimate_root=="k"){
            constraintConsistent = initConstraint(opt, nodes);
        }
        if (opt->e>0) calculateOutliers(opt,nodes,median_rate,true);
        if (opt->splitExternal) splitExternalBranches(opt,nodes);
        if (!opt->constraint){//LD without constraints
            if (!constraintConsistent){
                ostringstream oss;
                oss<<"- There's conflict in the input temporal constraints.\n";
                opt->warningMessage.push_back(oss.str());
            }
            if (opt->estimate_root==""){//keep the given root
                cout<<"Dating without temporal constraints ..."<<endl;
                without_constraint_multirates(opt,nodes,true);
                output(br,y,io, opt,nodes,*(io->outResult),*(io->outTree2),*(io->outTree3),r,diffTopology);
            }
            else if (opt->estimate_root=="k"){
                cout<<"Estimating the root position on the branch defined by given outgroups ..."<<endl;
                vector<int>::iterator iter=nodes[0]->suc.begin();
                int s1=(*iter);
                iter++;
                int s2=(*iter);
                br=nodes[s1]->B+nodes[s2]->B;
                nodes[s1]->V=variance(opt,br);
                nodes[s2]->V=nodes[s1]->V;
                constraintConsistent = without_constraint_active_set_lambda_multirates(true,br,opt,nodes,true);
                if (constraintConsistent){
                    output(br,y,io, opt,nodes,*(io->outResult),*(io->outTree2),*(io->outTree3),r,diffTopology);
                }
                else{
                    myExit("There's conflict or not enough signal in the input temporal constraints.\n");
                }
            }
            else{//estimate the root
                if (opt->estimate_root.compare("l")==0){//improve locally the root around the given root
                    cout<<"Estimating the root position locally around the given root ..."<<endl;
                    r=estimate_root_without_constraint_local_rooted(opt,nodes);
                }
                else{ //forget the given root and re-estimate the position of the root over all branhces
                    cout<<"Estimating the root position on all branches ..."<<endl;
                    r=estimate_root_without_constraint_rooted(opt,nodes);
                }
                if (r>0){
                    Node** nodes_new = cloneLeaves(opt,nodes,0);
                    vector<int>::iterator iter=nodes[0]->suc.begin();
                    int s1=(*iter);
                    iter++;
                    int s2=(*iter);
                    for (int i=opt->nbINodes; i<=opt->nbBranches; i++) {
                        nodes_new[i]->status=nodes[i]->status;
                    }
                    reroot_rootedtree(br,r,s1,s2,opt,nodes,nodes_new);
                    without_constraint_active_set_lambda_multirates(true,br,opt,nodes_new,true);
                    output(br,y,io, opt,nodes_new,*(io->outResult),*(io->outTree2),*(io->outTree3),r,diffTopology);
                    for (int i=0;i<opt->nbBranches+1;i++) delete nodes_new[i];
                    delete[] nodes_new;
                }
            }
        }
        else {//QPD with temporal constrains
            imposeMinBlen(*(io->outResult),opt,nodes,minB,true);
            if (constraintConsistent || (opt->estimate_root!="" && opt->estimate_root!="k")){
                if (constraintConsistent){
                    if (opt->estimate_root==""){//keep the given root
                        if (constraintConsistent){
                            cout<<"Dating under temporal constraints mode ..."<<endl;
                            constraintConsistent = with_constraint_multirates(opt,nodes,true);
                        }
                        if (constraintConsistent) {
                            output(br,y,io, opt,nodes,*(io->outResult),*(io->outTree2),*(io->outTree3),r,diffTopology);
                        }
                        else{
                            myExit("There's conflict or not enough signal in the input temporal constraints.\n");
                        }
                    }
                    else if (opt->estimate_root=="k"){
                        
                        if (constraintConsistent){
                            cout<<"Estimating the root position on the branch defined by given outgroups ..."<<endl;
                            vector<int>::iterator iter=nodes[0]->suc.begin();
                            int s1=(*iter);
                            iter++;
                            int s2=(*iter);
                            br=nodes[s1]->B+nodes[s2]->B;
                            nodes[s1]->V=variance(opt,br);
                            nodes[s2]->V=nodes[s1]->V;
                            constraintConsistent = with_constraint_active_set_lambda_multirates(true,br,opt,nodes,true);
                        }
                        if (constraintConsistent) {
                            output(br,y,io, opt,nodes,*(io->outResult),*(io->outTree2),*(io->outTree3),r,diffTopology);
                        } else {
                            myExit("There's conflict or not enough signal in the input temporal constraints.\n");
                        }
                    }
                    else{//estimate the root
                        if (opt->estimate_root.compare("l")==0){//improve locally the root around the given root
                            cout<<"Estimating the root position locally around the given root ..."<<endl;
                            r=estimate_root_with_constraint_local_rooted(opt,nodes);
                        }
                        else if (opt->estimate_root.compare("a")==0){
                            //forget the given root and re-estimate the position of the root over all branhces using fast method
                            cout<<"Estimating the root position on all branches using fast method ..."<<endl;
                            r=estimate_root_with_constraint_fast_rooted(opt,nodes);
                        }
                        else{ //forget the given root and re-estimate the position of the root over all branches
                            cout<<"Estimating the root position on all branches ..."<<endl;
                            r=estimate_root_with_constraint_rooted(opt,nodes);
                        }
                        Node** nodes_new = cloneLeaves(opt,nodes,0);
                        vector<int>::iterator iter=nodes[0]->suc.begin();
                        int s1=(*iter);
                        iter++;
                        int s2=(*iter);
                        for (int i=opt->nbINodes; i<=opt->nbBranches; i++) {
                            nodes_new[i]->status=nodes[i]->status;
                        }
                        reroot_rootedtree(br,r,s1,s2,opt,nodes,nodes_new);
                        with_constraint_active_set_lambda_multirates(true,br,opt,nodes_new,true);
                        output(br,y,io, opt,nodes_new,*(io->outResult),*(io->outTree2),*(io->outTree3),r,diffTopology);
                        for (int i=0;i<opt->nbBranches+1;i++) delete nodes_new[i];
                        delete[] nodes_new;
                    }
                } else {
                    myExit("There's conflict or not enough signal in the input temporal constraints.\n");
                }
            } else {
                myExit("There's conflict or not enough signal in the input temporal constraints.\n");
            }
        }
        for (int i=0;i<=opt->nbBranches;i++) delete nodes[i];
        delete[] nodes;
    }
    delete opt;
    elapsed_time = (double)(clock()-start)/CLOCKS_PER_SEC;
    *(io->outResult)<<"\n*********************************************************\n";
    *(io->outResult)<<"\nTOTAL ELAPSED TIME: "<<elapsed_time<<" seconds\n";
    *(io->outTree2)<<"End;\n";
    cout<<"\nTOTAL ELAPSED TIME: "<<elapsed_time<<" seconds"<<endl;
    if (!inputOutput){
        delete io;
    }
    return EXIT_SUCCESS;
}

#ifndef USE_LSD2
int main( int argc, char** argv ) {
    return buildTimeTree(argc, argv);
}
#endif
