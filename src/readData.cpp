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
#include "readData.h"

Node** tree2data(FILE * tree,Pr* pr,int & s){
    //checkRooted(pr);
    int inode = 1;//id of internal nodes;
    int countleaf=1;//n;
    stack<int> pileNode;
    char c = readBracket(tree,"input tree");
    int a=1;
    s=1;
    int nbChild=0;
    vector<Node*> internal_nodes;
    vector<Node*> leaves;
    internal_nodes.push_back(new Node());
    leaves.push_back(new Node());
    do{
        c = readChar(tree,"input tree");
        if (c==')'){
            a--;
            nbChild=0;
            s=1;
            Node* new_inode = new Node();//nodes[pr->nbINodes-inode]=new Node();
            new_inode->L = readSupport(tree,"input tree");//nodes[pr->nbINodes-inode]->L=readSupport(tree,"input tree");
            //stack<int> listSuc;
            while (!pileNode.empty() && s!=0){
                s=pileNode.top();pileNode.pop();
                //cout<<s<<" "<<inode<<endl;
                if (s>0){//s!=-1){
                    internal_nodes[s]->P=inode;//nodes[s]->P=pr->nbINodes-inode;
                } else{
                    leaves[-s]->P=inode;
                }
                //listSuc.push(s);
                nbChild++;
            }
            /*while (!listSuc.empty()) {
                s=listSuc.top();
                listSuc.pop();
                nodes[inode]->suc.push_back(s);//nodes[pr->nbINodes-inode]->suc.push_back(s);
            }*/
            if (a>0) new_inode->B=readdouble(tree,"input tree");//nodes[pr->nbINodes-inode]->B=readdouble(tree,"input tree");
            pileNode.push(inode);//pileNode.push(pr->nbINodes-inode);
            inode++;
            internal_nodes.push_back(new_inode);
        }
        else if (c!='(' && c!=',' && c!=-1 && c!=';' && c!='\n'){
            string lb=readLabel(c,tree,a);
            pileNode.push(-countleaf);//pileNode.push(pr->nbBranches-countleaf);
            Node* new_leaf = new Node();//nodes[pr->nbBranches-countleaf]=new Node();
            new_leaf->L=lb;//nodes[pr->nbBranches-countleaf]->L=lb;
            new_leaf->B=readdouble(tree,"input tree");//nodes[pr->nbBranches-countleaf]->B=readdouble(tree,"input tree");
            leaves.push_back(new_leaf);
            countleaf++;//countleaf++;
        }
        else if (c=='(') {a++;pileNode.push(0);}
        else if (c=='\n') {
            c=readChar(tree,"input tree");
        }
    }  while (a>0);
    
    //Re-organize internal & leaves indices
    Node** nodes;
    if (nbChild==2) {
        pr->rooted=true;
        pr->nbINodes = internal_nodes.size()-1;
        pr->nbBranches = leaves.size() + internal_nodes.size() - 3;
        nodes = new Node*[pr->nbBranches+1];
        for (int i=0;i< pr->nbINodes;i++){
            nodes[i] = internal_nodes[pr->nbINodes-i];
            nodes[i]->P = pr->nbINodes-internal_nodes[pr->nbINodes-i]->P;
        }
        for (int i=pr->nbINodes;i<=pr->nbBranches;i++){
            nodes[i] = leaves[pr->nbBranches-i+1];
            nodes[i]->P = pr->nbINodes-leaves[pr->nbBranches-i+1]->P;
        }
        nodes[0]->P=-1;
        nodes[0]->B=-1;
    }
    else{
        pr->rooted=false;
        pr->nbINodes = internal_nodes.size();
        pr->nbBranches = leaves.size() + internal_nodes.size() - 2;
        if (s>0) s = pr->nbINodes-s;
        else s = pr->nbBranches+s+1;
        nodes = new Node*[pr->nbBranches+1];
        for (int i=1;i< pr->nbINodes;i++){
            nodes[i] = internal_nodes[pr->nbINodes-i];
            nodes[i]->P = pr->nbINodes-internal_nodes[pr->nbINodes-i]->P;
        }
        for (int i=pr->nbINodes;i<=pr->nbBranches;i++){
            nodes[i] = leaves[pr->nbBranches-i+1];
            nodes[i]->P = pr->nbINodes-leaves[pr->nbBranches-i+1]->P;
        }
        nodes[1]->P=-1;
        nodes[1]->B=-1;
        nodes[0] = new Node();
    }
    delete internal_nodes[0];
    delete leaves[0];
    if (!pr->rooted){
        if (pr->estimate_root=="" && pr->fnOutgroup=="") {
            cout<<"The input trees are not rooted, use either option -g to specify the outgroups file or -r to estimate the root"<<endl;
            exit(EXIT_FAILURE);
        }
    }
    if (!pr->rooted && pr->fnOutgroup=="" && pr->estimate_root!="" && pr->estimate_root.compare("l")==0) {
        cout<<"The input trees are not rooted, so the \"l\" method for rooting function is turned to the \"a\" method."<<endl;
        pr->estimate_root="a";
    }
    return nodes;
}
void readDateFile(Pr* pr,Node** &nodes,bool& constraintConsistent){
    int lineNb=getLineNumber(pr->inDateFile);
    FILE * dateFile = fopen(pr->inDateFile.c_str(),"rt");
    int ino=readInt(dateFile,"Error in the date file, the file should begin with an integer (the number of temporal constrains)");
    if (lineNb-1<ino) {
        cout<<"The number of given constraints is small than the number of constraints to read. Please change the number of constraints to read at the first line of the input date file."<<endl;
        exit(EXIT_FAILURE);
    }
    string w1="";//store temporal constraints that are not in the tree;
    string w2="";//store nodes that have more than 1 temporal constraints.
    for (int i=0;i<ino;i++){
        string s=readWord(dateFile,pr->inDateFile);
        int type='n';
        double v1=0;
        double v2=0;
        int k = getPosition(nodes,s,0,pr->nbBranches+1);
        vector<int> mr;
        string ld=s;
        if (k==-1 && (s.compare("mrca")==0)){
            char c='(';
            ld="";
            while (c!=')'){
                string ss="";
                c=readCommaBracket(dateFile,pr->inDateFile,ss);
                int k1=getPosition(nodes,ss,0,pr->nbBranches+1);
                if (k1!=-1){
                    mr.push_back(k1);
                    if (ld=="") ld=ld+ss;
                    else ld=ld+","+ss;
                }
            }
            ld="mrca("+ld+")";
            if (mr.size()>0){
                k=mrca(nodes,mr);
            }
        }
        if (k!=-1){
            if (nodes[k]->type!='n'){
                w2=w2+" "+s;
            }
            char c = readChar(dateFile,pr->inDateFile);
            while (c<33 || c>126) c=readChar(dateFile,pr->inDateFile);
            if (c=='l' || c=='L' || c=='u' || c=='U' || c=='b' || c=='B'){//interval value
                char p = readChar(dateFile,pr->inDateFile);
                if (p=='('){
                    if (c=='l' || c=='L'){
                        type='l';
                        v1=readdouble(dateFile,pr->inDateFile);
                    }
                    else if (c=='u' || c=='U'){
                        type='u';
                        v1=readdouble(dateFile,pr->inDateFile);
                    }
                    else if (c=='b' || c=='B'){
                        type='b';
                        v1=readdouble(dateFile,pr->inDateFile);
                        if (readChar(dateFile,pr->inDateFile)==','){
                            v2=readdouble(dateFile,pr->inDateFile);
                        }
                        else{
                            cout<<"date constraint of type 'b' must have two values"<<endl;
                            exit(EXIT_FAILURE);
                        }
                        if (v1>v2) {
                            double t=v1;
                            v1=v2;
                            v2=t;
                        }
                        if (v1==v2) {
                            type='p';
                        }
                        else {
                            type='b';
                        }
                    }
                    c=readChar(dateFile,pr->inDateFile);
                    while (c<33 || c>126) c=readChar(dateFile,pr->inDateFile);
                }
                else{
                    cout<<"Error reading "<<pr->inDateFile<<" file: calibration point must be defined as either 'l(lower_bound)' or 'u(upper_bound)' or 'b(lower_bound,upper_bound)'"<<endl;
                    exit(EXIT_FAILURE);
                }
            }
            else {//precise value
                string wd="";
                wd+=c;
                while (fscanf(dateFile,"%c",&c)==1 && c>=33 && c<=126) {
                    wd+=c;
                }
                v1=atof(wd.c_str());
                type='p';
            }
            Date* newdate;
            if (mr.size()>0){
                newdate = new Date(ld,type,v1,v2,mr);
            }
            else{
                newdate = new Date(ld,type,v1,v2,k);
            }
            if (k<pr->nbINodes){
                pr->internalConstraints.push_back(newdate);
            }
            else{
                bool bl = nodes[k]->addConstraint(newdate);
                delete newdate;
                constraintConsistent=constraintConsistent && bl;
            }
        }
        else{
            w1=w1+" "+s;
            if (i<ino-1){
                char c=readChar(dateFile,pr->inDateFile);
                while (c!='\n') c=readChar(dateFile,pr->inDateFile);
            }
        }
    }
    if (w1!=""){
        std::ostringstream oss;
        if (pr->verbose) oss<<"- The nodes"+w1+" in the input date file are not present in the input tree.\n";
        else oss<<"- There are some nodes in the input date file that are not present in the input tree.\n";
        pr->warningMessage.push_back(oss.str());
    }
    if (w2!=""){
        std::ostringstream oss;
        if (pr->verbose) oss<<"- There nodes"+w2+" have more than one temporal constraint.\n";
        else oss<<"- There nodes that have more than one temporal constraint.\n";
        pr->warningMessage.push_back(oss.str());
    }
    fclose(dateFile);
}


void readPartitionFile(Pr* pr){
    //Read partition file
    ifstream partFile(pr->partitionFile.c_str());
    string line;
    pr->multiplierRate.push_back(1);
    while (getline(partFile,line)) {
        int pos=0;
        string groupName = readWord(line,pos);
        double m = readDouble(line,pos);
        if (m<=0){
            cout<<"Error in the partition file: after the group name must be a positive real which is the prior proportion of the group rate compared to the main rate. Put 1 if you don't known about this value. Selecting appropriate value helps to converge faster."<<endl;
            exit(EXIT_FAILURE);
        }
        pr->multiplierRate.push_back(m);
        int s = line.find_first_of("{",0);
        int e = line.find_first_of ("}",s+1);
        Part* part = new Part(groupName);
        while (s<e && s!=-1) {
            string term = line.substr(s+1,(e-s-1));
            int ss = 0;
            int ee = ss;
            Subtree* subtree;
            bool first = true;
            while (ee!=term.length()) {
                Pair* node;
                ss = firstCharacter(term,ss);
                ee = lastCharacter(term,ss);
                if (term.at(ss)=='[' && term.at(ee-1)==']') {
                    node = new Pair(true,term.substr(ss+1,ee-ss-2));
                }
                else {
                    node = new Pair(false,term.substr(ss,ee-ss));
                }
                if (first) {
                    subtree = new Subtree(node);
                    first=false;
                }
                else{
                    subtree->tips.push_back(node);
                }
                ss=ee+1;
            }
            part->subtrees.push_back(subtree);
            s = line.find_first_of("{",e+1);
            e = line.find_first_of("}",s+1);
        }
        pr->ratePartition.push_back(part);
    }
}


list<int> path(Pr* pr, Node** nodes,int s,int t){
    int cs=s;
    int ct=t;
    list<int> froms;
    stack<int> fromt;
    bool stops=false;
    bool stopt=false;
    while (!stops || !stopt){
        if (!stops) {
            froms.push_back(cs);
            cs=nodes[cs]->P;
            if (isAncestor(nodes, cs, t)) {
                stops=true;
            }
        }
        if (!stopt) {
            fromt.push(ct);
            ct=nodes[ct]->P;
            if (isAncestor(nodes, ct, s)) {
                stopt=true;
            }
        }
    }
    concat(froms,fromt);
    return froms;
}


int getBranchOut(Pr* pr,Node** nodes,list<string> &outgroups,bool &keepBelow){
    list<string>::iterator iter=outgroups.begin();
    if (outgroups.size()==1) {
        keepBelow=false;
        return getPosition(nodes,*iter,0,pr->nbBranches+1);
    }
    else{
        list<int> out;
        while (iter!=outgroups.end()) {
            int t=getPosition(nodes,*iter,0,pr->nbBranches+1);
            if (t!=-1) {
                out.push_back(t);
            }
            iter++;
        }
        if (out.size()>0) {
            int s=pr->nbINodes;
            while (s<=pr->nbBranches && contain(s , out)) {
                s++;
            }
            list<int>::iterator t=out.begin();
            list<int> common = path(pr,nodes,s,*t);
            t++;
            while (t!=out.end()) {
                list<int> p = path(pr,nodes,s,*t);
                common = intersect(common,p);
                t++;
            }
            int nbPass0=0;
            list<int>::iterator iter=common.begin();
            int r=*iter;
            while (iter!=common.end()){
                r=*iter;
                if (nodes[r]->P==0){
                    nbPass0++;
                }
                iter++;
            }
            if (nbPass0==2) {
                keepBelow=false;
            }
            else{
                keepBelow=true;
            }
            return r;
        }
        else{
            return -1;
        }
    }
}

void extrait_outgroup(Pr* pr,list<string> &outgroups){
    FILE * tree = fopen(pr->inFile.c_str(),"rt");
    int s =0;
    if (tree==NULL) cout<<"Can not open the tree file"<<endl;
    else{
        string newFile = pr->inFile;
        if (pr->keepOutgroup) newFile+=".reroot";
        else newFile+=".ingroup";
        FILE* w=fopen(newFile.c_str(),"wt");
        for (int y=0;y<pr->nbData;y++){
            printf("Removing outgroups of tree %d ...\n",y+1);
            Node** nodes = tree2data(tree,pr,s);
            if (!pr->rooted) {
                nodes=unrooted2rootedS(pr, nodes, s);
            } else{
                computeSuc_polytomy(pr, nodes);
            }
            bool keepBelow;
            int r=getBranchOut(pr, nodes, outgroups,keepBelow);
            if (r!=-1) {
                Node** nodes_new = cloneLeaves(pr,nodes,0);
                int p_r=reroot_rootedtree(r, pr, nodes, nodes_new);
                computeSuc_polytomy(pr, nodes_new);
                if (pr->keepOutgroup) {
                    fprintf(w,"%s",newick(0,0,pr,nodes_new).c_str());
                }
                else{
                    if (keepBelow) {
                        fprintf(w,"%s",newick(r, r, pr, nodes_new).c_str());
                    }
                    else{
                        fprintf(w,"%s",newick(p_r, p_r,pr, nodes_new).c_str());
                    }
                }
                for (int i=0;i<=pr->nbBranches;i++) delete nodes_new[i];
                delete[] nodes_new;
            }
            else{
                cout<<"The outgroups are not in the tree "<<y+1<<endl;
                fprintf(w,"%s",newick(0,0,pr,nodes).c_str());
            }
            for (int i=0;i<=pr->nbBranches;i++) delete nodes[i];
            delete[] nodes;
        }
        cout<<"The new input trees are writen to the file "<<newFile<<endl;
        fclose(w);
        pr->inFile=newFile;
    }
    fclose(tree);
}

