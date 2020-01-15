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
#include "options.h"

Pr* getOptions( int argc, char** argv )
{
    if( argc>1 )
        return getCommandLine( argc, argv );
    else
        return getInterface();
}

Pr* getCommandLine( int argc, char** argv)
{
    Pr* opt = new Pr();
    int c;
    string s;
    bool iflag = false,
    dflag = false,
    flagA=false,
    flagZ=false,
    fflag=false;
    while ( (c = getopt(argc, argv, ":i:d:o:s:n:g:r:v:ct:w:b:ha:z:f:kje:m:p:V")) != -1 )
    {
        switch (c)
        {
            case 'i':
                if( access( optarg, R_OK )!=0 )
                    myExit( "Cannot read the file named \"%s\"\n", optarg );
                opt->inFile = optarg;
                iflag = true;
                break;
            case 'd':
                if( access( optarg, R_OK )!=0 )
                    myExit( "Cannot read the file named \"%s\"\n", optarg );
                opt->inDateFile = optarg;
                dflag = true;
                break;
            case 'p':
                if( access( optarg, R_OK )!=0 )
                    myExit( "Cannot read the file named \"%s\"\n", optarg );
                opt->partitionFile = optarg;
                break;
            case 'o':
                opt->outFile = optarg;
                break;
            case 's':
                if( !isInteger(optarg) )
                    myExit("Argument of option -s must be an integer.\n");
                opt->seqLength = atoi(optarg);
                if( opt->seqLength<1 )
                    myExit("Argument of option -s must be strictly positive.\n");
                break;
            case 'n':
                if( !isInteger(optarg) )
                    myExit("Argument of option -n must be an integer.\n");
                opt->nbData = atoi(optarg);
                if( opt->nbData<1 )
                    myExit("Argument of option -n must be strictly positive.\n");
                break;
            case 'g':
                if( access( optarg, R_OK )!=0 )
                    myExit( "Cannot read the file named \"%s\"\n", optarg );
                opt->fnOutgroup = optarg;
                break;
            case 'k':
                if (opt->fnOutgroup!=""){
                	opt->keepOutgroup=true;
                	opt->estimate_root="k";
                }
                break;
            case 'r':
                opt->estimate_root = optarg;
                if (opt->estimate_root.compare("l")!=0 && opt->estimate_root.compare("a")!=0 && opt->estimate_root.compare("as")!=0 && opt->estimate_root.compare("k")!=0)
                    myExit("Argument of option -r must be either \"l\" or \"a\" or \"as\".\n");
                break;
            case 'v':
                if( !isInteger(optarg) )
                    myExit("Argument of option -v must be either 1 or 2.\n");
                opt->variance = atoi(optarg);
                if (opt->variance!=1 && opt->variance!=2){
                    myExit("Argument of option -v must be either 1 (for using orginal branches to compute variance) or 2 (LSD will be run twice, the second time uses the variances based on the estimated branch lengths of the first time).\n");
                }
                break;
            case 'c':
                opt->constraint = true;
                break;
            case 'b':
                if( !isInteger(optarg) )
                    myExit("Argument of option -b must be an integer.\n");
                opt->c = atoi( optarg );
                if (opt->c<0)
                    myExit("Argument of option -b must not be negative.\n");
                break;
            case 'e':
                if( !isReal(optarg) )
                    myExit("Argument of option -e must be a real.\n");
                opt->e = atof(optarg);
                if (opt->e<0){
                    std::ostringstream oss;
                    oss<<"- The specified argument of option e was negative, so outliers detection option was not processed\n";
                    opt->warningMessage.push_back(oss.str());
                }
                break;
            case 'm':
                if( !isInteger(optarg) )
                    myExit("Argument of option -m must be an integer.\n");
                opt->m = atof( optarg );
                if (opt->m<0) myExit("Argument of option -m can not be negative.\n");
                break;
            case 't':
                if( !isReal(optarg) )
                    myExit("Argument of option -t must be a real.\n");
                opt->rho_min = atof(optarg);
                if (opt->rho_min<=0)
                    myExit("Argument of option -t must be strictly positive.\n");
                break;
            case 'w':
                if( access( optarg, R_OK )!=0 )
                    myExit( "Cannot read the file named \"%s\"\n", optarg );
                opt->rate = optarg;
                opt->givenRate[0] = true;
                break;
            case 'a':
                if( !isInteger(optarg) )
                    myExit("Argument of option -a must be an integer.\n");
                opt->mrca=atof(optarg);
                flagA=true;
                break;
            case 'z':
                if( !isInteger(optarg) )
                    myExit("Argument of option -z must be an integer.\n");
                opt->leaves=atof(optarg);
                flagZ=true;
                break;
            case 'h':
                printHelp();
                exit( EXIT_SUCCESS );
                break;
            case 'j':
                opt->verbose = true;
                break;
            case 'f':
                if( !isInteger(optarg) )
                    myExit("Argument of option -f must be the number of samplings to compute confidence intervals, e.g. 100 ...\n");
                opt->nbSampling = atoi( optarg );
                if (opt->nbSampling<0)
                    myExit("Argument of option -f must be a positive integer.\n");
                opt->ci=true;
                fflag = true;
                break;
            case 'V':
                printf("lsd2 " VERSION"\n");
                exit( EXIT_SUCCESS );
            case '?':
                myExit("Unrecognized option: -%c\n", optopt);
            case ':':
                if (optopt=='v') myExit("Argument of option -v must be either 1 (for using orginal branches to compute variance) or 2 (for using branches estimated by LSD to compute variance, i.e LSD will be run 2 times: the first time is just used to compute the variance).\n");
                else myExit("Option -%c requires an operand\n", optopt );
            default:
                myExit("?? getopt returned character code 0%o ??\n", c );
        }
    }
    if( !(iflag) )
        myExit("Argument -i is necessary to continue.\n");
    if (!dflag && (flagA && !flagZ)){
        myExit("The input date file is not provided, so option -z is required to use with option -a to estimate relative dates.\n");
    }
    if (!dflag && (!flagA && flagZ)){
        myExit("The input date file is not provided, so option -a is required to use with option -z to estimate relative dates.\n");
    }
    if (!dflag && (!flagA && !flagZ)){
        opt->relative=true;
        opt->mrca=0;
        opt->leaves=1;
    }
    if (!dflag && (flagA && flagZ)){
        if (opt->mrca >= opt->leaves)
            myExit("The root date must be strictly smaller than the tips date.\n");
        opt->relative=true;
    }
    if (dflag)
        opt->relative=false;
    
    if (!opt->constraint && opt->estimate_root.compare("as")==0){
        cout<<"The non constrained mode is chosen, so the \"as\" method for rooting function is the same as the \"a\" method."<<endl;
        opt->estimate_root="a";
    }
    if( opt->outFile=="") opt->outFile = opt->inFile + ".result";
    
    return opt;
}

Pr* getInterface()
{
    Pr* opt = new Pr();
    opt->inFile = getInputFileName("Enter your Input Tree File name> ");
    checkRooted(opt);
    cout<<"Do you have a date file? y/n "<<endl;
    char letter[3];
    do{
        fgets(letter,3,stdin);
        if (*letter=='n' || *letter=='N') {
            cout<<"There is no date file, so the program will estimate relative dates with root date = 0 and tips date = 1. Type 'y' to continue or 'n' to modify the root date and the tips date"<<endl;
            char letter1[3];
            do {
                fgets( letter1, 3, stdin );
                if (*letter1=='n' || *letter1=='N'){
                    do {
                        opt->mrca = getInputReal("Enter the root date (default=0)> ");
                        opt->leaves = getInputReal("Enter the tips date (default=1)> ");
                        if (opt->leaves <= opt->mrca) cout<<"Root date must be smaller than the tips date."<<endl;
                    } while (opt->leaves <= opt->mrca);
                }
                else if (*letter1=='y' || *letter1=='Y'){
                    opt->mrca=0;
                    opt->leaves=1;
                }
                else {
                    cout<<"Type 'y' to continue or 'n' to modify the root date and tips date"<<endl;
                }
                opt->relative=true;
            } while (*letter1!='n' && *letter1!='N' && *letter1!='y' && *letter1!='Y');
        }
        else if (*letter=='y' || *letter=='Y'){
            opt->inDateFile = getInputFileName("Enter you input date file name>");
            opt->relative=false;
        }
    } while (*letter!='n' && *letter!='N' && *letter!='y' && *letter!='Y');
    opt->outFile = opt->inFile+".result";
    if (!opt->rooted) {
        opt->estimate_root="a";
    }
    do
    {
        printInterface( stdout, opt);
        cout<<endl;
        fgets( letter, 3, stdin );
        if( isOptionActivate( opt, *letter ) ) setOptionsWithLetter( opt, *letter);
    } while( *letter!='y' && *letter!='Y' );
    return opt;
}


void printInterface( FILE* in, Pr* opt)
{
    fprintf(in,"\nLEAST-SQUARE METHODS TO ESTIMATE RATES AND DATES - " VERSION" \n\n");
    fprintf(in,"\nInput files:\n");
    fprintf(in,"  I                                 Input tree file : %s\n",opt->inFile.c_str());
    if (opt->relative==true)
        fprintf(in,"  D                         Estimate relative dates : mrca date = %.2f, tips date =%.2f\n",opt->mrca,opt->leaves);
    else
        fprintf(in,"  D                                 Input date file : %s\n",opt->inDateFile.c_str());
    fprintf(in,"  P                                  Partition file : ");
    if (opt->partitionFile=="")        fprintf(in,"No\n");
    else fprintf(in,"%s\n",opt->partitionFile.c_str());
    fprintf(in,"Output file:\n");
    fprintf(in,"  O                                    Output file  : %s\n",opt->outFile.c_str());
    fprintf(in,"Parameters:\n");
    fprintf(in,"  C                                With constraints : ");
    if (!opt->constraint) fprintf(in,"No\n");
    else {
        fprintf(in,"Yes\n");
    }
    fprintf(in,"  T                        Lower bound for the rate : %.3e\n",opt->rho_min);
    fprintf(in,"  V                                  With variances : ");
    if (opt->variance==0) fprintf(in,"No\n");
    else {
        if (opt->variance==1) fprintf(in,"Yes, use variances based on original branches\n");
        else if (opt->variance==2) fprintf(in,"Yes, use variances based on branches estimated by LSD\n");
        fprintf(in,"  B                Adjusted parameter for variances : %d\n",opt->c);
    }
    fprintf(in,"  R                               Estimate the root : ");
    if (opt->estimate_root==""){
        fprintf(in,"No\n");
    }
    else if (opt->estimate_root.compare("l")==0){
        fprintf(in,"Around the given root\n");
    }
    else if (opt->estimate_root.compare("a")==0 && opt->constraint){
        fprintf(in,"Use fast method to search on all branches\n");
    }
    else if (opt->estimate_root.compare("a")==0 && !opt->constraint){
        fprintf(in,"Search on all branches\n");
    }
    else if (opt->estimate_root.compare("as")==0){
        fprintf(in,"Use constrained mode on all branches\n");
    }
    fprintf(in,"  W                         Given substitution rate : ");
    if (opt->rate=="") fprintf(in,"No\n");
    else fprintf(in,"%s\n",opt->rate.c_str());
    if (opt->fnOutgroup=="")
        fprintf(in,"  G                                 Given outgroups : No\n");
    else {
        fprintf(in,"  G                         File contains outgroups : %s\n",opt->fnOutgroup.c_str());
        if (opt->keepOutgroup) {
            fprintf(in,"  K           Keep outgroups in the estimating tree :  Yes\n");
        }
        else{
            fprintf(in,"  K           Keep outgroups in the estimating tree : No\n");
        }
    }
    fprintf(in,"  N                               Multiple data set : ");
    if( opt->nbData< 2 )
        fprintf(in,"No\n");
    else
        fprintf(in,"Yes, %i data sets\n",opt->nbData);
    fprintf(in,"  F                    Compute confidence intervals : ");
    if (opt->ci){
        fprintf(in,"Yes, sampling %d times\n",opt->nbSampling);
        fprintf(in,"  S                                 Sequence Length : %i\n",opt->seqLength);
    }
    else
        fprintf(in,"No\n");
    fprintf(in,"  E                            Exclude outlier tips : ");
    if (opt->e>0){
        fprintf(in,"Yes, detecting and excluding outliers from the analysis\n");
        fprintf(in,"  M     Number of sampling nodes to detect outliers : %i\n",opt->m);
        fprintf(in,"  E         The Zscore threshold to detect outliers : %f\n",opt->e);
    }
    else
        fprintf(in,"No\n");
    
    fprintf(in,"\n  H to print Help ");
    fprintf(in,"\n  Y to accept or type a letter to change an option (X = Exit) ");
}

void printHelp( void )
{
    printf(BOLD"LSD: LEAST-SQUARES METHODS TO ESTIMATE RATES AND DATES - " VERSION"\n\n");
    printf(BOLD"DESCRIPTION\n"
           FLAT"\tThis program estimates the rate and the dates of the input phylogenies given some temporal constraints.\n"
           FLAT"\tIt minimizes the square errors of the branch lengths under normal distribution model.\n\n"
           );
    printf(BOLD"SYNOPSIS\n"
           FLAT"\t" BOLD"./lsd " FLAT"[" BOLD"-i " LINE"inputFile" FLAT"] "
           FLAT"[" BOLD"-d " LINE"inputDateFile" FLAT"] "
           FLAT"[" BOLD"-o " LINE"outputFile" FLAT"] "
           FLAT"[" BOLD"-c" FLAT"] "
           FLAT"[" BOLD"-v " LINE"mode" FLAT"] "
           FLAT"[" BOLD"-s " LINE"sequenceLength" FLAT"] "
           FLAT"[" BOLD"-n " LINE"datasetNumber" FLAT"]\n"
           FLAT"\t     [" BOLD"-t " LINE"lowerBoundRate" FLAT"] "
           FLAT"[" BOLD"-r " LINE"rootingMethod" FLAT"] "
           FLAT"[" BOLD"-b " LINE"varianceParameter" FLAT"] "
           FLAT"[" BOLD"-w " LINE"givenRateFile" FLAT"] "
           FLAT"[" BOLD"-g " LINE"outgroupFile" FLAT"] "
	   FLAT"[" BOLD"-f " LINE"nbSamplings" FLAT"] "
           FLAT"[" BOLD"-h" FLAT"]\n"
           FLAT"\n");
    
    printf(BOLD"OPTIONS\n"
           FLAT"\t" BOLD"-a " LINE"rootDate\n"
           FLAT"\t   If the dates of all tips are equal (which is given by option -z), you must use this option to provide the root date.\n"
           FLAT"\t   In this case, the input date file can be omitted, and the program estimates only the relative dates based on the given\n"
           FLAT"\t   root date and tips date. By default, T[root]=0 and T[tips]=1.\n"
           FLAT"\t" BOLD"-b " LINE"varianceParameter\n"
           FLAT"\t   The parameter (between 1 and 100 - which is the proportion to the sequence length) to compute the variances in option -v.\n"
           FLAT"\t   By default b=10. The smaller it is the more variance would be linear to branch length, which is relevant for strict clock.\n"
           FLAT"\t   The bigger it is the less effect of branch length on variance which might be better for relaxed clock.\n"
           FLAT"\t" BOLD"-c \n"
           FLAT"\t   By using this option, we impose the constraints that the date of every node is equal or smaller then\n"
           FLAT"\t   the dates of its descendants. Without constraints, the runtime is linear (LD). With constraints, the\n"
           FLAT"\t   problem is a quadratic programming and is solved efficiently (quasi-linear) by the active-set method.\n"
           FLAT"\t" BOLD"-d " LINE"inputDateFile\n"
           FLAT"\t   This options is used to read the name of the input date file which contains temporal constraints of internal nodes\n"
           FLAT"\t   or tips. An internal node can be defined either by its label (given in the input tree) or by a subset of tips that have it as \n"
           FLAT"\t   the most recent common ancestor (mrca).\n"
           FLAT"\t   The first line of this file is the number of temporal constraints. A temporal constraint can be a real, or a \n"
           FLAT"\t   lower bound " FLAT LINE"l(value)" FLAT", or an upper bound " FLAT LINE"u(value)" FLAT", or an interval " LINE"b(v1,v2)\n"
           FLAT"\t   For example, if the input tree has 4 taxa a,b,c,d, and an internal node named n, then following is a possible date file:\n"
           FLAT"\t    6\n"
           FLAT"\t    a l(2003)\n"
           FLAT"\t    b u(2007)\n"
           FLAT"\t    c 2005\n"
           FLAT"\t    d b(2001,2007)\n"
           FLAT"\t    mrca(a,b,c,d) b(2000,2001)\n"
           FLAT"\t    n l(2004)\n"
           FLAT"\t   If this option is omitted, the program will estimate relative dates by giving T[root]=0 and T[tips]=1.\n"
           FLAT"\t" BOLD"-e " LINE"ZscoreOutlier\n"
           FLAT"\t   This option is used to estimate and exclude outlier nodes before dating process.\n"
           FLAT"\t   LSD2 normalize the branch residus and decide a node is outlier if its related residus is great than the " LINE"ZscoreOutlier.\n"
           FLAT"\t   A normal value of " LINE"ZscoreOutlier" FLAT"could be 3, but you can adjust it bigger/smaller depending if you want to have\n"
           FLAT"\t   less/more outliers. Note that for now, some functionalities could not be combined with outliers estimation, for example \n"
           FLAT"\t   estimating multiple rates, imprecise date constraints.\n"
           FLAT"\t" BOLD"-f " LINE"samplingNumberCI\n"
           FLAT"\t   This option calculates the confidence intervals of the estimated rate and dates. The branch lengths of the esimated\n"
           FLAT"\t   tree are sampled " FLAT LINE"samplingNumberCI" FLAT " times to generate a set of simulated trees. We use Poisson distribution \n"
           FLAT"\t   to generate branch lengths so we need to multiple the branch lengths with the sequence length provided by option -s.\n"
           FLAT"\t   However, to avoid over-estimate the confidence intervals in the case of very long sequence length but not necessarily\n"
           FLAT"\t   strict molecular clock, we use min of 1000 and the real sequence length to generate simulated trees.\n"
	   FLAT"\t   confidence intervals are inferred from these simulated trees and are written in the nexus tree output with the name \"CI\".\n"
           FLAT"\t" BOLD"-g " LINE"outgroupFile\n"
           FLAT"\t   If your data contain outgroups, then specify the name of the outgroup file here.\n"
           FLAT"\t   The program will use the outgroups to root the trees. It will keep the outgroup in the trees if option -k is used with, \n"
           FLAT"\t   otherwise it will remove the outgroups. The format of this file should be:\n"
           FLAT"\t        n\n"
           FLAT"\t        OUTGROUP1\n"
           FLAT"\t        OUTGROUP2\n"
           FLAT"\t        ...\n"
           FLAT"\t        OUTGROUPn\n"
           FLAT"\t" BOLD"-h " LINE"help\n"
           FLAT"\t   Print this message.\n"
           FLAT"\t" BOLD"-i " LINE"inputTreesFile\n"
           FLAT"\t   The name of the input trees file. It contains tree(s) in newick format, each tree on one line. Note that the taxa sets of all\n"
           FLAT"\t   trees must be the same.\n"
           FLAT"\t" BOLD"-j\n"
           FLAT"\t   Verbose mode for output messages.\n"
           FLAT"\t" BOLD"-k\n"
           FLAT"\t   Use this option to keep the outgroups (given in option -g) in the estimated tree. The root position is then estimated on the\n"
           FLAT"\t   branch determined by the outgroups. If this option is not used, the outgroups will be removed.\n"
           FLAT"\t" BOLD"-m " LINE"samplingNumberOutlier\n"
           FLAT"\t   The number of dated nodes to be sampled when detecting outlier nodes. This should be smaller than the number of dated nodes,\n"
           FLAT"\t   and is 10 by default.\n"
           FLAT"\t" BOLD"-n " LINE"datasetNumber\n"
           FLAT"\t   The number of trees that you want to read and analyse.\n"
           FLAT"\t" BOLD"-o " LINE"outputFile\n"
           FLAT"\t   The base name of the output files to write the results and the time-scale trees.\n"
           FLAT"\t" BOLD"-p " LINE"partitionFile\n"
           FLAT"\t   The file that defines the partition of branches into multiple subsets in the case that you know each subset has a different rate.\n"
           FLAT"\t   In the partition file, each line contains the name of the group, the prior proportion of the group rate compared to the main rate\n"
           FLAT"\t   (selecting an appropriate value for this helps to converge faster), and a list of subtrees whose branches are supposed to have the \n"
           FLAT"\t   same substitution rate. All branches that are not assigned to any subtree form a group having another rate.\n"
           FLAT"\t   A subtree is defined between {}: its first node corresponds to the root of the subtree, and the following nodes (if there any) \n"
           FLAT"\t   correspond to the tips of the subtree. If the first node is a tip label then it takes the mrca of all tips as the root of the subtree.\n"
           FLAT"\t   If the tips of the subtree are not defined (so there's only the defined root), then by \n"
           FLAT"\t   default this subtree is extended down to the tips of the full tree. For example the input tree is \n"
           FLAT"\t   ((A:0.12,D:0.12)n1:0.3,((B:0.3,C:0.5)n2:0.4,(E:0.5,(F:0.2,G:0.3)n3:0.33)n4:0.22)n5:0.2)root;\n"
           FLAT"\t   and you have the following partition file:\n"
           FLAT"\t         group1 1 {n1} {n5 n4}\n"
           FLAT"\t         group2 1 {n3}\n"
           FLAT"\t   then there are 3 rates: the first one includes the branches (n1,A), (n1,D), (n5,n4), (n5,n2), (n2,B), (n2,C); the second one \n"
           FLAT"\t   includes the branches (n3,F), (n3,G), and the last one includes all the remaining branches. If the internal nodes don't have labels,\n"
           FLAT"\t   then they can be defined by mrca of at least two tips, for example n1 is mrca(A,D)\n"
           FLAT"\t" BOLD"-r " LINE"rootingMethod\n"
           FLAT"\t   This option is used to specify the rooting method to estimate the position of the root for unrooted trees, or\n"
           FLAT"\t   re-estimate the root for rooted trees. The principle is to search for the position of the root that minimizes\n"
           FLAT"\t   the objective function.\n"
           FLAT"\t   Use " BOLD"-r l" FLAT" if your tree is rooted, and you want to re-estimate the root locally around the given root.\n"
           FLAT"\t   Use " BOLD"-r a" FLAT" if you want to estimate the root on all branches (ignoring the given root if the tree is rooted).\n"
           FLAT"\t       In this case, if the constrained mode is chosen (option -c), method \"a\" first estimates the root without using the constraints.\n"
           FLAT"\t       After that, it uses the constrained mode to improve locally the position of the root around this pre-estimated root.\n"
           FLAT"\t   Use " BOLD"-r as" FLAT" if you want to estimate to root using constrained mode on all branches.\n"
           FLAT"\t   Use " BOLD"-r k" FLAT" if you want to re-estimate the root position on the same branche of the given root.\n"
           FLAT"\t       If combined with option -g, the root will be estimated on the branche defined by the outgroups.\n"
           FLAT"\t" BOLD"-s " LINE"sequenceLength\n"
           FLAT"\t   This option is used to specify the length of the multiple alignments that were used to build the input trees. It is used to  \n"
           FLAT"\t   compute the confidence intervals with option -f. By default it is 1000 and it is 1000 if the sequence length is > 1000.\n"
           FLAT"\t" BOLD"-t " LINE"rateLowerBound\n"
           FLAT"\t   This option corresponds to the lower bound for the estimating rate. It is 1e-10 by default.\n"
           FLAT"\t" BOLD"-v " LINE"variance\n"
           FLAT"\t   Use this option if you want to apply variances for the branch lengths in order to recompense big errors on long estimated branch lengths. \n"
           FLAT"\t   The variance of the branch Bi is Vi = (Bi+b/100) where b is specified by option -b.\n"
           FLAT"\t   If " FLAT LINE"variance" FLAT"=1, then LSD uses the input branch lengths to calculate variances. If " FLAT LINE"variance" FLAT"=2, then LSD\n"
           FLAT"\t   runs twice where the second time it calculates the variances based on the estimated branch lengths of the first run. However -v 2 only \n"
           FLAT"\t   improves the result in the case variances were well estimated in the first run, most of the case it means your data follows a strict clock. \n"
           FLAT"\t   If your tree is likely relaxed, don't use -v 2.\n"
           FLAT"\t" BOLD"-V \n"
           FLAT"\t   Get the actual version.\n"
           FLAT"\t" BOLD"-w " LINE"givenRte\n"
           FLAT"\t   This option is used to specify the name of the file containing the substitution rates.\n"
           FLAT"\t   In this case, the program will use the given rates to estimate the dates of the nodes.\n"
           FLAT"\t   This file should have the following format\n"
           FLAT"\t        RATE1\n"
           FLAT"\t        RATE2\n"
           FLAT"\t        ...\n"
           FLAT"\t  where RATEi is the rate of the tree i in the inputTreesFile.\n"
           FLAT"\t" BOLD"-z " LINE"tipsDate\n"
           FLAT"\t   This option is used to give the date of the tips when they are all equal. It must be used with option -a to give the\n"
           FLAT"\t   root date. In this case the input date file can be omitted, and the program estimates only the relative dates based on\n"
           FLAT"\t   the given root date and tips date. By default, T[root]=0 and T[tips]=1.\n"
           );
}
string getInputString(string msg)
{
    string s;
    cout<<msg<<endl;
    cin>>s;
    return s;
}

string getInputFileName( string msg)
{
    string outfile;
    do
    {
        outfile = getInputString(msg);
        if( access(outfile.c_str(), F_OK)==0 || string(outfile).length()==0)
            break;
        printf( "The file \"%s\" does not exist.\n", outfile.c_str() );
    } while( true );
    if( access(outfile.c_str(), R_OK)!=0 && string(outfile).length()!=0)
        myExit("Could not access to the file named \"%s\" in reading.\n", outfile.c_str() );
    return outfile;
}

string getOutgroupFileName( string msg)
{
    string outfile;
    do
    {
        outfile = getInputString(msg);
        if( access(outfile.c_str(), F_OK)==0 || outfile.compare("")==0)
            break;
        printf( "The file \"%s\" does not exist.\n", outfile.c_str() );
    } while( true );
    if (outfile.compare("")==0) {
        return outfile;
    }
    if( access(outfile.c_str(), R_OK)!=0 )
        myExit("Could not access to the file named \"%s\" in reading.\n", outfile.c_str() );
    return outfile;
}

double getInputReal( string msg )
{
    string word;
    do
    {
        word = getInputString( msg );
        if( isReal(word.c_str()) )
            break;
        myErrorMsg("Your word is not recognized as a real.\n");
    } while( true );
    return atof( word.c_str() );
}

double getInputPositiveReal( string msg )
{
    string word;
    do
    {
        word = getInputString( msg );
        if( isReal(word.c_str()) && atof( word.c_str())>0)
	  break;
        myErrorMsg("Your word is not recognized as a strictly positive real.\n");
    } while( true );
    return atof( word.c_str() );
}

int getInputInteger( string msg )
{
    string word;
    do
    {
        word = getInputString( msg );
        if( isInteger(word.c_str()) )
            break;
        myErrorMsg("Your word is not recognized as an integer.\n");
    } while( true );
    return atoi( word.c_str() );
}

int getPositiveInputInteger( string msg )
{
    int i;
    do
    {
        i = getInputInteger(msg);
        if( i>0 )
            break;
        myErrorMsg("It must be a strictly positive integer.\n");
    } while( true );
    return i;
}

bool isOptionActivate( Pr* opt, char l )
{
    switch(l)
    {
        case 'i':
        case 'I':
        case 'd':
        case 'D':
        case 'o':
        case 'O':
        case 'p':
        case 'P':
        case 's':
        case 'S':
        case 'c':
        case 'C':
        case 'v':
        case 'V':
        case 'b':
        case 'B':
        case 'r':
        case 'R':
        case 'g':
        case 'G':
        case 'k':
        case 'K':
        case 't':
        case 'T':
        case 'w':
        case 'W':
        case 'f':
        case 'F':
        case 'n':
        case 'N':
        case 'y':
        case 'Y':
        case 'q':
        case 'Q':
        case 'e':
        case 'E':
        case 'm':
        case 'M':
        case 'j':
        case 'J':
        case 'x':
        case 'X':
        case 'h':
        case 'H':
        return true;
    }
    return false;
}

void setOptionsWithLetter( Pr* opt, char letter )
{
    //char* fnOut;
    switch( letter )
    {
        case 'x':
        case 'X':
            exit( EXIT_SUCCESS );
        case 'i':
        case 'I':
            opt->inFile = getInputFileName("Enter your Input File name> ");
            checkRooted(opt);
            break;
        case 'd':
        case 'D':
            cout<<"Do you have a date file? y/n "<<endl;
            char letter[3];
            do{
                fgets(letter,3,stdin);
                if (*letter=='n' || *letter=='N') {
                    cout<<"There is no date file, so the program will estimate relative dates with root date = 0 and tips date = 1. Type 'y' to continue or 'n' to modify the root date and the tips date"<<endl;
                    char letter1[3];
                    do {
                        fgets( letter1, 3, stdin );
                        if (*letter1=='n' || *letter1=='N'){
                            do {
                                opt->mrca = getInputReal("Enter the root date (default=0)> ");
                                opt->leaves = getInputReal("Enter the tips date (default=1)> ");
                                if (opt->leaves <= opt->mrca) cout<<"Root date must be smaller than the tips date."<<endl;
                            } while (opt->leaves <= opt->mrca);
                        }
                        else if (*letter1=='y' || *letter1=='Y'){
                            opt->mrca=0;
                            opt->leaves=1;
                        }
                        else {
                            cout<<"Type 'y' to continue or 'n' to modify the root date and tips date"<<endl;
                        }
                        opt->relative=true;
                    } while (*letter1!='n' && *letter1!='N' && *letter1!='y' && *letter1!='Y');
                }
                else if (*letter=='y' || *letter=='Y'){
                    opt->inDateFile = getInputFileName("Enter you input date file name>");
                    opt->relative=false;
                }
            } while (*letter!='n' && *letter!='N' && *letter!='y' && *letter!='Y');
            break;
        case 'p':
        case 'P':
            opt->partitionFile = getInputFileName("Enter your Partition File name> ");
            break;
        case 's':
        case 'S':
            opt->seqLength = getPositiveInputInteger("Enter your sequence length> ");
            break;
        case 'n':
        case 'N':
            opt->nbData = getPositiveInputInteger("Enter your number of dataset> ");
            break;
        case 'o':
        case 'O':
            opt->outFile=getInputString("Enter your output file name > ");
            while( access( opt->outFile.c_str(), F_OK )==0){
                cout<<"File "<<opt->outFile<<" already exists. Do you want to overwrite it? Y/N"<<endl;
                char letter[3];
                fgets( letter, 3, stdin );
                if (*letter=='N' || *letter=='n'){
                    opt->outFile = getInputString("Enter your output file name > ");
                }
                if (*letter=='Y' || *letter=='y') {
                    break;
                }
            }
            break;
        case 'c':
        case 'C':
            opt->constraint=!opt->constraint;
            break;
        case 'v':
        case 'V':
            if (opt->variance==0 || opt->variance==1) opt->variance+=1;
            else if (opt->variance==2) opt->variance=0;
            break;
        case 'b':
        case 'B':
            if (opt->variance) opt->c = getPositiveInputInteger("Enter the parameter for the variances> ");
            break;
        case 'r':
        case 'R':
            if (opt->rooted || opt->fnOutgroup!=""){
                if (opt->estimate_root==""){
                    opt->estimate_root="l";
                }
                else if (opt->estimate_root.compare("l")==0){
                    opt->estimate_root="a";
                }
                else if (opt->estimate_root.compare("a")==0 && opt->constraint){
                    opt->estimate_root="as";
                }
                else if (opt->estimate_root.compare("a")==0 && !opt->constraint) opt->estimate_root="";
                else if (opt->estimate_root.compare("as")==0) opt->estimate_root="";
            }
            else {
                if (opt->estimate_root==""){
                    opt->estimate_root="a";
                }
                else if (opt->estimate_root.compare("a")==0 && opt->constraint){
                    opt->estimate_root="as";
                }
                else if (opt->estimate_root.compare("a")==0 && !opt->constraint){
                    cout<<"The trees are not rooted, you must use either option -g to specify the outgroups file or -r to estimate the root"<<endl;
                }
                else if (opt->estimate_root.compare("as")==0){
                    opt->estimate_root="a";
                }
            }
            break;
        case 'g':
        case 'G':
            if (opt->fnOutgroup==""){
                string fo = getOutgroupFileName("Enter the name of the file that contains your outgroups> ");
                if (fo.compare("")!=0) {
                    opt->fnOutgroup=fo;
                    opt->estimate_root="";
                }
            }
            else{
                opt->fnOutgroup="";
                if (!opt->rooted && opt->estimate_root==""){
                    opt->estimate_root="a";
                }
            }
            break;
        case 'k':
        case 'K':
            if (opt->fnOutgroup!=""){
            	opt->keepOutgroup=!opt->keepOutgroup;
            	if (opt->keepOutgroup) {
                	opt->estimate_root="k";
            	}
	    }
            break;
        case 'w':
        case 'W':
            if (opt->rate=="") {
                opt->rate = getInputFileName("Enter the name of the file that contains the rates> ");
                if (string(opt->rate).length()==0) opt->rate="";
                else opt->givenRate[0]=true;
            }
            else{
                opt->rate="";
                opt->givenRate[0]=false;
            }
            break;
        case 't':
        case 'T':
            opt->rho_min = getInputPositiveReal("Enter the lower bound for the rate> ");
            break;
        case 'f':
        case 'F':
            if (!opt->ci) {
                opt->nbSampling = getInputInteger("Enter the number of sampling for calculating confidence intervals> ");
                opt->ci=true;
            }
            else{
                opt->ci = false;
            }
            break;
        case 'j':
        case 'J':
            opt->verbose = !opt->verbose;
            break;
        case 'm':
        case 'M':
            if (opt->e>0){
                opt->m = getInputInteger("Enter the number of sampling dated nodes for outliers detection> ");
            }
            break;
        case 'e':
        case 'E':
            opt->e = getInputReal("Enter the Zscore threshold for outliers detection> ");
        case 'h':
        case 'H':
            printHelp();
            break;
    }
}




