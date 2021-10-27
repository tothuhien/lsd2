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
    const string VERSION="v.2.3";
    Pr* opt = new Pr();
    int c;
    string s;
    bool iflag = false,
    dflag = false,
    flagA=false,
    flagZ=false,
    sflag=false,
    fflag=false,
    uflag=false,
    Uflag=false,
    vflag=false,
    lflag=false,
    validDate = true;
    while ( (c = getopt(argc, argv, ":i:d:D:o:s:n:g:r:v:Ft:w:b:ha:z:f:Gje:m:p:q:u:l:U:R:S:EV")) != -1 )
    {
        switch (c)
        {
            case 'i':
#ifndef USE_LSD2
                if( access( optarg, R_OK )!=0 )
                    myExit( "Cannot read the file named \"%s\"\n", optarg );
#endif
                opt->inFile = optarg;
                iflag = true;
                break;
            case 'd':
#ifndef USE_LSD2
                if( access( optarg, R_OK )!=0 )
                    myExit( "Cannot read the file named \"%s\"\n", optarg );
#endif
                opt->inDateFile = optarg;
                dflag = true;
                break;
            case 'D':
                if( !isInteger(optarg) )
                    myExit("Argument of option -D must be an integer.\n");
                opt->outDateFormat = atoi(optarg);
                if (opt->outDateFormat !=1 && opt->outDateFormat != 2 && opt->outDateFormat != 3){
                    myExit("Argument of option -D must be either 1 (date as real) or 2 (date as YY-MM-DD)\n or 3 (date as YY-MM).\n");
                }
                break;
            case 'p':
#ifndef USE_LSD2
                if( access( optarg, R_OK )!=0 )
                    myExit( "Cannot read the file named \"%s\"\n", optarg );
#endif
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
                sflag = true;
                break;
            case 'n':
                if( !isInteger(optarg) )
                    myExit("Argument of option -n must be an integer.\n");
                opt->nbData = atoi(optarg);
                if( opt->nbData<1 )
                    myExit("Argument of option -n must be strictly positive.\n");
                break;
            case 'g':
#ifndef USE_LSD2
                if( access( optarg, R_OK )!=0 )
                    myExit( "Cannot read the file named \"%s\"\n", optarg );
#endif
                opt->fnOutgroup = optarg;
                break;
            case 'G':
                opt->removeOutgroup=true;
                break;
            case 'r':
                opt->estimate_root = optarg;
                if (opt->estimate_root.compare("l")!=0 && opt->estimate_root.compare("a")!=0 && opt->estimate_root.compare("as")!=0 && opt->estimate_root.compare("k")!=0)
                    myExit("Argument of option -r must be either \"l\" or \"a\" or \"as\".\n");
                break;
            case 'v':
                if( !isInteger(optarg) )
                    myExit("Argument of option -v must be an integer either 0, 1, or 2.\n");
                opt->variance = atoi(optarg);
                if (opt->variance!=0 && opt->variance!=1 && opt->variance!=2){
                    myExit("Argument of option -v must be either 0 (do not use variance) 1 (default value,\n to using orginal branches to compute variance) or 2 (LSD will be run twice, the\n second time uses the variances based on the estimated branch lengths of the\n first time).\n");
                }
                vflag = true;
                break;
            case 'F':
                opt->constraint = false;
                break;
            case 'b':
                if( !isReal(optarg) )
                    myExit("Argument of option -b must be a real.\n");
                opt->c = atof( optarg );
                if (opt->c<=0)
                    myExit("Argument of option -b must be a positive number, see the help\n page for more information.\n");
                break;
            case 'e':
                if( !isReal(optarg) )
                    myExit("Argument of option -e must be a real.\n");
                opt->e = atof(optarg);
                if (opt->e<0){
                    std::ostringstream oss;
                    oss<<"- The specified argument of option e was negative, so\n outliers detection option was not processed\n";
                    opt->warningMessage.push_back(oss.str());
                }
                break;
            case 'm':
                if( !isInteger(optarg) )
                    myExit("Argument of option -m must be an integer.\n");
                opt->m = atof( optarg );
                if (opt->m<2) myExit("Argument of option -m must be >= 2.\n");
                break;
            case 't':
                if( !isReal(optarg) )
                    myExit("Argument of option -t must be a real.\n");
                opt->rho_min = atof(optarg);
                if (opt->rho_min<=0)
                    myExit("Argument of option -t must be strictly positive.\n");
                break;
            case 'q':
                if (!isReal(optarg))
                    myExit("Argument of option -q must be a real.\n");
                opt->q = atof(optarg);
                if (opt->q<0)
                    myExit("Argument of option -q could not be negative.\n");
                break;
            case 'w':
                if( access( optarg, R_OK )!=0 )
                    myExit( "Cannot read the file named \"%s\"\n", optarg );
                opt->rate = optarg;
                opt->givenRate[0] = true;
                break;
            case 'a':
                opt->MRCA = optarg;
                flagA=true;
                break;
            case 'z':
                opt->LEAVES = optarg;
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
                if ( !isInteger(optarg) ){
                    if( access( optarg, R_OK )!=0 ){
                        myExit( "Argument of option -f must be an integer or a file containing bootstraps trees\n" );
                    }
                    opt->bootstraps_file = optarg;
                } else{
                    opt->nbSampling = atoi( optarg );
                    if (opt->nbSampling <= 0){
                        myExit("Argument of option -f must be a positive integer.\n");
                    }
                    fflag = true;
                }
                opt->ci=true;
                break;
            case 'u':
                if ( !isReal(optarg)){
                    if (!strcmp(optarg,"e")){
                        opt->minblen = -1;
                        uflag = true;
                    } else {
                        myExit("Argument of option -u must be a real or letter e\n");
                    }
                } else{
                    opt->minblen = atof(optarg);
                    if (opt->minblen <0){
                        myExit("Argument of option -u must be >= 0\n");
                    }
                }
                break;
            case 'U':
                if ( !isReal(optarg)){
                    if (!strcmp(optarg,"e")){
                        opt->minblenL = -1;
                        uflag = true;
                    } else {
                        myExit("Argument of option -U must be a real or letter e\n");
                    }
                } else{
                    opt->minblenL = atof(optarg);
                    Uflag = true;
                    if (opt->minblenL <0){
                        myExit("Argument of option -U must be >= 0\n");
                    }
                }
                break;
            case 'l':
                if ( !isReal(optarg))
                    myExit("Argument of option -l must be a real\n");
                opt->nullblen = atof(optarg);
                lflag = true;
                break;
            case 'R':
                if ( !isReal(optarg))
                    myExit("Argument of option -R must be a real\n");
                opt->round_time = atof(optarg);
                if (opt->round_time<=0)
                    myExit("Argument of option -R must be positive\n");
                break;
            case 'S':
                if ( !isReal(optarg))
                    myExit("Argument of option -S must be a real\n");
                opt->support = atof(optarg);
                if (opt->round_time<0)
                    myExit("Argument of option -S must be non negative\n");
                break;
            case 'E':
                opt->splitExternal = true;
                break;
            case 'V':
                cout<<"lsd2 "<<VERSION<<endl;
                exit( EXIT_SUCCESS );
            case '?':
                if (optopt == 'c') myExit("Unrecognized option -c, temporal constraint now becomes by default without this\n option.\n");
                else myExit("Unrecognized option: -%c\n", optopt);
            case ':':
                if (optopt=='v') myExit("Argument of option -v must be either 1 (for using orginal branches to compute\n variance) or 2 (for using branches estimated by LSD to compute variance, i.e LSD\n will be run 2 times: the first time is just used to compute the variance).\n");
                else myExit("Option -%c requires an operand\n", optopt );
            default:
                myExit("?? getopt returned character code 0%o ??\n", c );
        }
    }
    if( !(iflag) )
        myExit("Argument -i is necessary to continue.\n");
    if ( !(sflag) && (fflag || uflag || !lflag || vflag))
        myExit("Argument -s is necessary to continue.\n");
    if (!dflag && !flagA && !flagZ){
        //opt->relative=true;
        opt->MRCA="0";
        opt->LEAVES="1";
        cout<<"Not any input date provided. The results correspond to the estimation of\n relative dates when T[mrca]=0 and T[tips]=1"<<endl;
    }
    if (!dflag && flagA && !flagZ){
        myExit("Tips date are required via option -z\n");
    }
    if (!dflag && !flagA && flagZ){
        myExit("Root date is required via option -a\n");
    }
    if (!lflag){
        opt->nullblen = 0.5/opt->seqLength;
    }
    if (opt->minblen < 0 && !Uflag){
        opt->minblenL = -1;
    }
    if (!opt->constraint && opt->estimate_root.compare("as")==0){
        cout<<"The non constrained mode is chosen, so the \"as\" method for rooting function\n is the same as the \"a\" method."<<endl;
        opt->estimate_root="a";
    }
    if (opt->estimate_root.compare("")==0 && !opt->removeOutgroup && opt->fnOutgroup != ""){
        opt->estimate_root="k";
    }
    if( opt->outFile=="") opt->outFile = opt->inFile + ".result";
    opt->treeFile2=opt->outFile+".date.nexus";
    opt->treeFile3=opt->outFile+".nwk";
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
            cout<<"There is no date file, so the program will estimate relative dates with root\n date = 0 and tips date = 1.\n Type 'y' to continue or 'n' to modify the root\n date and the tips date"<<endl;
            char letter1[3];
            do {
                fgets( letter1, 3, stdin );
                if (*letter1=='n' || *letter1=='N'){
                    opt->MRCA = getInputString("Enter the root date (default=0)> ");
                    opt->LEAVES = getInputString("Enter the tips date (default=1)> ");
                }
                else if (*letter1=='y' || *letter1=='Y'){
                    opt->MRCA = "0";
                    opt->LEAVES = "1";
                    opt->mrca=0;
                    opt->leaves=1;
                }
                else {
                    cout<<"Type 'y' to continue or 'n' to modify the root date and tips date"<<endl;
                }
            } while (*letter1!='n' && *letter1!='N' && *letter1!='y' && *letter1!='Y');
        }
        else if (*letter=='y' || *letter=='Y'){
            opt->inDateFile = getInputFileName("Enter you input date file name>");
        }
    } while (*letter!='n' && *letter!='N' && *letter!='y' && *letter!='Y');
    opt->outFile = opt->inFile+".result";
    opt->treeFile2=opt->outFile+".date.nexus";
    opt->treeFile3=opt->outFile+".nwk";
    if (!opt->rooted) {
        opt->estimate_root="a";
    }
    opt->seqLength = getPositiveInputInteger("Enter the length of sequences that were used to build your tree>");
    fgets( letter, 3, stdin );
    printInterface(std::cout, opt);
    do
    {
        fgets( letter, 3, stdin );
        if( isOptionActivate( opt, *letter ) ) setOptionsWithLetter( opt, *letter);
        printInterface(std::cout, opt);
        cout<<endl;
    } while(( *letter!='y' && *letter!='Y'));
    if (opt->nullblen<0){
        opt->nullblen = 0.5/opt->seqLength;
    }
    if (opt->inDateFormat==2 && opt->outDateFormat==0) opt->outDateFormat=2;
    return opt;
}


void printInterface(ostream& in, Pr* opt)
{
    const string VERSION = "v.2.3";

    in<<"\nLEAST-SQUARE METHODS TO ESTIMATE RATES AND DATES - "<<VERSION<<" \n\n";
    in<<"\nInput files:\n";
    in<<"  i                                               Input tree file : "<<opt->inFile.c_str()<<"\n";
    in<<"  d                                               Input date file : ";
    if (opt->inDateFile!=""){
        in<<opt->inDateFile.c_str()<<"\n";
    } else {
        in<<"No\n";
    }
    in<<"  p                                                Partition file : ";
    if (opt->partitionFile=="")        in<<"No\n";
    else in<<opt->partitionFile.c_str()<<"\n";
    if (opt->fnOutgroup=="")
        in<<"  g                                               Given outgroups : No\n";
    else {
        in<<"  g                                       File contains outgroups : "<<opt->fnOutgroup.c_str()<<"\n";
        if (opt->removeOutgroup) {
            in<<"  G                       Remove outgroups in the estimating tree : Yes\n";
        }
        else{
            in<<"  G                       Remove outgroups in the estimating tree : No\n";
        }
    }
    in<<"Output file:\n";
    in<<"  o                                                  Output file  : "<<opt->outFile.c_str()<<"\n";
    in<<"Parameters:\n";
    in<<"  a                                                     Root date : ";
    if (opt->MRCA!=""){
        in<<opt->MRCA.c_str()<<"\n";
    } else {
        in<<"No\n";
    }
    in<<"  z                                                     Tips date : ";
    if (opt->LEAVES!=""){
        in<<opt->LEAVES.c_str()<<"\n";
    } else {
        in<<"No\n";
    }
    in<<"  c                                              With constraints : ";
    if (!opt->constraint) in<<"No\n";
    else {
        in<<"Yes\n";
    }
    in<<"  t                                      Lower bound for the rate : "<<opt->rho_min<<"\n";
    in<<"  v                                                With variances : ";
    if (opt->variance==0) in<<"No\n";
    else {
        if (opt->variance==1) in<<"Yes, use variances based on input branch lengths\n";
        else if (opt->variance==2) in<<"Yes, use variances based on estimated branch lengths\n";
        in<<"  b                              Adjusted parameter for variances : ";
        if (opt->c==-1){
            in<<"To estimate\n";
        } else{
            in<<opt->c<<"\n";
        }
    }
    in<<"  r                                             Estimate the root : ";
    if (opt->estimate_root=="k"){
        in<<"On the branch given by the outgroups\n";
    }
    else if (opt->estimate_root.compare("l")==0){
        in<<"Around the given root\n";
    }
    else if (opt->estimate_root.compare("a")==0 && opt->constraint){
        in<<"Use fast method to search on all branches\n";
    }
    else if (opt->estimate_root.compare("a")==0 && !opt->constraint){
        in<<"Search on all branches\n";
    }
    else if (opt->estimate_root.compare("as")==0){
        in<<"Use constrained mode on all branches\n";
    }
    else{
        in<<"No\n";
    }
    in<<"  w                                       Given substitution rate : ";
    if (opt->rate=="") in<<"No\n";
    else in<<opt->rate.c_str()<<"\n";
    in<<"  n                                             Multiple data set : ";
    if( opt->nbData< 2 )
        in<<"No\n";
    else
        in<<"Yes, "<<opt->nbData<<" data sets\n";
    in<<"  f                                  Compute confidence intervals : ";
    if (opt->ci && opt->bootstraps_file==""){
        in<<"Yes, sampling "<<opt->nbSampling<<" times\n";
        in<<"  q                  Standard deviation of lognormal relaxed clock: "<<opt->q<<" (for computing confidence intervals)\n";
    } else if (opt->bootstraps_file!=""){
        in<<"Use bootstrap trees from "<<opt->bootstraps_file<<"\n";
    } else
        in<<"No\n";
    if (opt->nullblen==-1 || opt->ci){
        in<<"  s                                               Sequence length : "<<opt->seqLength<<"\n";
    }
    in<<"  e                                          Exclude outlier tips : ";
    if (opt->e>0){
        in<<"Yes, detect and exclude outliers from the analysis\n";
        in<<"  m                   Number of sampling nodes to detect outliers : "<<opt->m<<"\n";
        in<<"  e                       The Zscore threshold to detect outliers : "<<opt->e<<"\n";
    }
    else
        in<<"No\n";
    in<<"  u                Min internal branch length of time scaled tree : ";
    if (opt->minblen==-1){
        in<<"To estimate\n";
        in<<"  R     Rounding number for min branch length of time scaled tree : ";
        if (opt->round_time>0) in<<opt->round_time<<"\n";
        else in<<"To guess\n";
    } else {
        in<<opt->minblen<<"\n";
    }
    in<<"  U                Min external branch length of time scaled tree : ";
    if (opt->minblenL==-1){
        if (opt->minblen>0) {
            opt->minblenL = opt->minblen;
            in<<opt->minblenL<<"\n";
        } else {
            in<<"To estimate\n";
        }
    } else {
        in<<opt->minblenL<<"\n";
    }
    in<<"  l                        Collapsed internal branch length limit : ";
    if (opt->bootstraps_file!=""){
        in<<"Don't collapse\n";
    } else if (opt->nullblen==-1){
        in<<0.5/opt->seqLength<<"\n";
    } else {
        in<<opt->nullblen<<"\n";
    }
    in<<"  D                                            Output date format :";
    if (opt->outDateFormat==0){
        in<<" Based on input date format\n";
    }
    if (opt->outDateFormat==1){
        in<<" Real number\n";
    }
    if (opt->outDateFormat==2){
        in<<" Year-Month-Day\n";
    }
    in<<"  E  Estimate rates for external and internal branches separately : ";
    if (opt->splitExternal){
        in<<"Yes\n";
    } else{
        in<<"No\n";
    }
    in<<"\n  h to print Help ";
    in<<"\n  y to accept or type a letter to change an option (x = Exit) ";
}

void printHelp( void )
{
    const string BOLD = "\033[00;01m";
    const string LINE = "\033[00;04m";
    const string FLAT = "\033[00;00m";
    const string VERSION = "v.2.3";
    
    cout<<BOLD<<"LSD: LEAST-SQUARES METHODS TO ESTIMATE RATES AND DATES - "<<VERSION<<"\n\n";
    cout<<BOLD<<"DESCRIPTION\n"
    <<FLAT<<"\tThis program estimates the rate and the dates of the input phylogenies given some temporal constraints.\n"
    <<FLAT<<"\tIt minimizes the square errors of the branch lengths under normal distribution model.\n\n";
    cout<<BOLD<<"SYNOPSIS\n"
    <<FLAT<<"\t" <<BOLD<<"./lsd " <<FLAT<<"[" <<BOLD<<"-i " <<LINE<<"inputFile" <<FLAT<<"] "
    <<FLAT<<"[" <<BOLD<<"-d " <<LINE<<"inputDateFile" <<FLAT<<"] "
    <<FLAT<<"[" <<BOLD<<"-o " <<LINE<<"outputFile" <<FLAT<<"] "
    <<FLAT<<"[" <<BOLD<<"-s " <<LINE<<"sequenceLength" <<FLAT<<"] "
    <<FLAT<<"[" <<BOLD<<"-g " <<LINE<<"outgroupFile" <<FLAT<<"] "
    <<FLAT<<"[" <<BOLD<<"-f " <<LINE<<"nbSamplings" <<FLAT<<"] "
    <<FLAT<<"\n";
    
    cout<<BOLD<<"OPTIONS\n"
    <<FLAT<<"\t" <<BOLD<<"-a " <<LINE<<"rootDate\n"
    <<FLAT<<"\t   To specify the root date if there's any. If the root date is not a number, but a string (ex: 2020-01-10, or b(2019,2020)) then it should\n"
    <<FLAT<<"\t   be put between the quotes.\n"
    <<FLAT<<"\t" <<BOLD<<"-b " <<LINE<<"varianceParameter\n"
    <<FLAT<<"\t   The parameter (between 0 and 1) to compute the variances in option -v. It is the pseudo positive constant to add to the branch lengths\n"
    <<FLAT<<"\t   when calculating variances, to adjust the dependency of variances to branch lengths. By default b is the maximum between median branch length\n"
    <<FLAT<<"\t   and 10/seqlength; but it should be adjusted  based on how/whether the input tree is relaxed or strict. The smaller it is the more variances\n"
    <<FLAT<<"\t   would be linear to branch lengths, which is relevant for strict clock. The bigger it is the less effect of branch lengths on variances, \n"
    <<FLAT<<"\t   which might be better for relaxed clock.\n"
    <<FLAT<<"\t" <<BOLD<<"-d " <<LINE<<"inputDateFile\n"
    <<FLAT<<"\t   This options is used to read the name of the input date file which contains temporal constraints of internal nodes\n"
    <<FLAT<<"\t   or tips. An internal node can be defined either by its label (given in the input tree) or by a subset of tips that have it as \n"
    <<FLAT<<"\t   the most recent common ancestor (mrca). A date could be a real or a string or format year-month-day.\n"
    <<FLAT<<"\t   The first line of this file is the number of temporal constraints. A temporal constraint can be fixed date, or a \n"
    <<FLAT<<"\t   lower bound " <<FLAT<<LINE<<"l(value)" <<FLAT<<", or an upper bound " <<FLAT<<LINE<<"u(value)" <<FLAT<<", or an interval " <<LINE<<"b(v1,v2)\n"
    <<FLAT<<"\t   For example, if the input tree has 4 taxa a,b,c,d, and an internal node named n, then following is a possible date file:\n"
    <<FLAT<<"\t    6\n"
    <<FLAT<<"\t    a l(2003.12)\n"
    <<FLAT<<"\t    b u(2007.07)\n"
    <<FLAT<<"\t    c 2005\n"
    <<FLAT<<"\t    d b(2001.2,2007.11)\n"
    <<FLAT<<"\t    mrca(a,b,c,d) b(2000,2001)\n"
    <<FLAT<<"\t    n l(2004.3)\n"
    <<FLAT<<"\t   If this option is omitted, and option -a, -z are also omitted, the program will estimate relative dates by giving T[root]=0 and T[tips]=1.\n"
    <<FLAT<<"\t" <<BOLD<<"-D " <<LINE<<"outDateFormat\n"
    <<FLAT<<"\t    Specify output date format: 1 for real, 2 for year-month-day. By default the program will guess the format of input dates and uses it for\n"
    <<FLAT<<"\t    output dates.\n"
    <<FLAT<<"\t" <<BOLD<<"-e " <<LINE<<"ZscoreOutlier\n"
    <<FLAT<<"\t   This option is used to estimate and exclude outlier nodes before dating process.\n"
    <<FLAT<<"\t   LSD2 normalize the branch residus and decide a node is outlier if its related residus is great than the " <<LINE<<"ZscoreOutlier.\n"
    <<FLAT<<"\t   A normal value of " <<LINE<<"ZscoreOutlier" <<FLAT<<"could be 3, but you can adjust it bigger/smaller depending if you want to have\n"
    <<FLAT<<"\t   less/more outliers. Note that for now, some functionalities could not be combined with outliers estimation, for example \n"
    <<FLAT<<"\t   estimating multiple rates, imprecise date constraints.\n"
    <<FLAT<<"\t" <<BOLD<<"-E\n"
    <<FLAT<<"\t   Use this option if you want to estimate the rate of internal and external branches separately.\n"
    <<FLAT<<"\t" <<BOLD<<"-f " <<LINE<<"samplingNumberCI or bootstrapTreeFile\n"
    <<FLAT<<"\t   This option calculates the confidence intervals of the estimated rate and dates. If the bootstrap trees file is specified, then dating is \n"
    <<FLAT<<"\t   processed on each of these trees, and the 2.5th and 97.5th percentile of calculated rates and dates are reported as confidence intervals.\n"
    <<FLAT<<"\t   If one of the bootstrap trees does not have the same topology as the original tree, then only confidence intervals of root dates are referred.\n"
    <<FLAT<<"\t   If there's no bootstrap trees, then we simulate a set of trees. In this case, the number of simulated trees should be specified. Those \n"
    <<FLAT<<"\t   trees have the same topology as the original one, and their branch lengths are generated as follow:\n"
    <<FLAT<<"\t          b_i = Poisson(B_i*seqLen)*lognormal(1,q)\n"
    <<FLAT<<"\t   where B_i is the branch length estimated by lsd2 and q is the standard deviation of relaxed clock settable via option -q.\n"
    <<FLAT<<"\t   So the first term represents a simulation for the number of substitutions on the branch i; and the second term adds the effect of relaxed\n"
    <<FLAT<<"\t   clock on that branch length.\n"
    <<FLAT<<"\t   Confidence intervals are written in the nexus tree with label CI_date, and can be visualzed with Figtree under Node bar feature with CI_height.\n"
    <<FLAT<<"\t" <<BOLD<<"-g " <<LINE<<"outgroupFile\n"
    <<FLAT<<"\t   If your data contain outgroups, then specify the name of the outgroup file here. The program will use the outgroups to root the trees.\n"
    <<FLAT<<"\t   If you use this combined with options -G, then the outgroups will be removed. The format of this file should be:\n"
    <<FLAT<<"\t        n\n"
    <<FLAT<<"\t        OUTGROUP1\n"
    <<FLAT<<"\t        OUTGROUP2\n"
    <<FLAT<<"\t        ...\n"
    <<FLAT<<"\t        OUTGROUPn\n"
    <<FLAT<<"\t" <<BOLD<<"-F \n"
    <<FLAT<<"\t   By default without this option, we impose the constraints that the date of every node is equal or smaller then the\n"
    <<FLAT<<"\t   dates of its descendants, so the running time is quasi-linear. Using this option we ignore this temporal constraints, and\n"
    <<FLAT<<"\t   the the running time becomes linear, much faster.\n"
    <<FLAT<<"\t" <<BOLD<<"-h " <<LINE<<"help\n"
    <<FLAT<<"\t   Print this message.\n"
    <<FLAT<<"\t" <<BOLD<<"-i " <<LINE<<"inputTreesFile\n"
    <<FLAT<<"\t   The name of the input trees file. It contains tree(s) in newick format, each tree on one line. Note that the taxa sets of all\n"
    <<FLAT<<"\t   trees must be the same.\n"
    <<FLAT<<"\t" <<BOLD<<"-j\n"
    <<FLAT<<"\t   Verbose mode for output messages.\n"
    <<FLAT<<"\t" <<BOLD<<"-G\n"
    <<FLAT<<"\t   Use this option to remove the outgroups (given in option -g) in the estimated tree. If this option is not used, the outgroups \n"
    <<FLAT<<"\t   will be kept and the root position in estimated on the branch defined by the outgroups.\n"
    <<FLAT<<"\t" <<BOLD<<"-l " <<LINE<<"nullBlen\n"
    <<FLAT<<"\t   A branch in the input tree is considered informative if its length is greater this value. By default it is 0.5/seq_length. Only \n"
    <<FLAT<<"\t   informative branches are forced to be bigger than a minimum branch length (see option -u for more information about this).\n"
    <<FLAT<<"\t" <<BOLD<<"-m " <<LINE<<"samplingNumberOutlier\n"
    <<FLAT<<"\t   The number of dated nodes to be sampled when detecting outlier nodes. This should be smaller than the number of dated nodes,\n"
    <<FLAT<<"\t   and is 10 by default.\n"
    <<FLAT<<"\t" <<BOLD<<"-n " <<LINE<<"datasetNumber\n"
    <<FLAT<<"\t   The number of trees that you want to read and analyse.\n"
    <<FLAT<<"\t" <<BOLD<<"-o " <<LINE<<"outputFile\n"
    <<FLAT<<"\t   The base name of the output files to write the results and the time-scale trees.\n"
    <<FLAT<<"\t" <<BOLD<<"-p " <<LINE<<"partitionFile\n"
    <<FLAT<<"\t   The file that defines the partition of branches into multiple subsets in the case that you know each subset has a different rate.\n"
    <<FLAT<<"\t   In the partition file, each line contains the name of the group, the prior proportion of the group rate compared to the main rate\n"
    <<FLAT<<"\t   (selecting an appropriate value for this helps to converge faster), and a list of subtrees whose branches are supposed to have the \n"
    <<FLAT<<"\t   same substitution rate. All branches that are not assigned to any subtree form a group having another rate.\n"
    <<FLAT<<"\t   A subtree is defined between {}: its first node corresponds to the root of the subtree, and the following nodes (if there any) \n"
    <<FLAT<<"\t   correspond to the tips of the subtree. If the first node is a tip label then it takes the mrca of all tips as the root of the subtree.\n"
    <<FLAT<<"\t   If the tips of the subtree are not defined (so there's only the defined root), then by \n"
    <<FLAT<<"\t   default this subtree is extended down to the tips of the full tree. For example the input tree is \n"
    <<FLAT<<"\t   ((A:0.12,D:0.12)n1:0.3,((B:0.3,C:0.5)n2:0.4,(E:0.5,(F:0.2,G:0.3)n3:0.33)n4:0.22)n5:0.2)root;\n"
    <<FLAT<<"\t   and you have the following partition file:\n"
    <<FLAT<<"\t         group1 1 {n1} {n5 n4}\n"
    <<FLAT<<"\t         group2 1 {n3}\n"
    <<FLAT<<"\t   then there are 3 rates: the first one includes the branches (n1,A), (n1,D), (n5,n4), (n5,n2), (n2,B), (n2,C); the second one \n"
    <<FLAT<<"\t   includes the branches (n3,F), (n3,G), and the last one includes all the remaining branches. If the internal nodes don't have labels,\n"
    <<FLAT<<"\t   then they can be defined by mrca of at least two tips, for example n1 is mrca(A,D)\n"
    <<FLAT<<"\t" <<BOLD<<"-q " <<LINE<<"standardDeviationRelaxedClock\n"
    <<FLAT<<"\t   This value is involved in calculating confidence intervals to simulate a lognormal relaxed clock. We multiply the simulated branch lengths\n"
    <<FLAT<<"\t   with a lognormal distribution with mean 1, and standard deviation q. By default q is 0.2. The bigger q is, the more your tree is relaxed\n"
    <<FLAT<<"\t   and give you bigger confidence intervals.\n"
    <<FLAT<<"\t" <<BOLD<<"-r " <<LINE<<"rootingMethod\n"
    <<FLAT<<"\t   This option is used to specify the rooting method to estimate the position of the root for unrooted trees, or\n"
    <<FLAT<<"\t   re-estimate the root for rooted trees. The principle is to search for the position of the root that minimizes\n"
    <<FLAT<<"\t   the objective function.\n"
    <<FLAT<<"\t   Use " <<BOLD<<"-r l" <<FLAT<<" if your tree is rooted, and you want to re-estimate the root locally around the given root.\n"
    <<FLAT<<"\t   Use " <<BOLD<<"-r a" <<FLAT<<" if you want to estimate the root on all branches (ignoring the given root if the tree is rooted).\n"
    <<FLAT<<"\t       In this case, if the constrained mode is chosen (option -c), method \"a\" first estimates the root without using the constraints.\n"
    <<FLAT<<"\t       After that, it uses the constrained mode to improve locally the position of the root around this pre-estimated root.\n"
    <<FLAT<<"\t   Use " <<BOLD<<"-r as" <<FLAT<<" if you want to estimate to root using constrained mode on all branches.\n"
    <<FLAT<<"\t   Use " <<BOLD<<"-r k" <<FLAT<<" if you want to re-estimate the root position on the same branche of the given root.\n"
    <<FLAT<<"\t       If combined with option -g, the root will be estimated on the branche defined by the outgroups.\n"
    <<FLAT<<"\t" <<BOLD<<"-R " <<LINE<<"round_time\n"
    <<FLAT<<"\t   This value is used to round the minimum branch length of the time scaled tree. The purpose of this is to make the minimum branch length\n"
    <<FLAT<<"\t   a meaningful time unit, such as day, week, year ... By default this value is 365, so if the input dates are year, the minimum branch\n"
    <<FLAT<<"\t   length is rounded to day. The rounding formula is round(R*minblen)/R.\n"
    <<FLAT<<"\t" <<BOLD<<"-s " <<LINE<<"sequenceLength\n"
    <<FLAT<<"\t   This option is used to specify the sequence length when estimating confidence intervals (option -f). It is used to generate \n"
    <<FLAT<<"\t   integer branch lengths (number of substitutions) by multiplying this with the estimated branch lengths. By default it is 1000.\n"
    <<FLAT<<"\t" <<BOLD<<"-S " <<LINE<<"minSupport\n"
    <<FLAT<<"\t   Together with collapsing internal short branches (see option -l), users can also collapse internal branches having weak support values (if\n"
    <<FLAT<<"\t   provided in the input tree) by using this option. The program will collapse all internal branches having support <= the specifed value.\n"
    <<FLAT<<"\t" <<BOLD<<"-t " <<LINE<<"rateLowerBound\n"
    <<FLAT<<"\t   This option corresponds to the lower bound for the estimating rate. It is 1e-10 by default.\n"
    <<FLAT<<"\t" <<BOLD<<"-u " <<LINE<<"minBlen\n"
    <<FLAT<<"\t   Using this option, lsd2 forces every branch of the time scaled tree to be >= minBlen. minBlen is either a positive real or letter e. \n"
    <<FLAT<<"\t   For the later case lsd2 estimates a minimum branch length, which is m/rate (m is the minimum branch length of the input tree, and \n"
    <<FLAT<<"\t   rate is an pre-estimated rate). This estimated value is then rounded to the number of days or weeks or years, using the rounding\n"
    <<FLAT<<"\t   parameter -R. Without this option minBlen = 0 by default.\n"
    <<FLAT<<"\t" <<BOLD<<"-U " <<LINE<<"minExBlen\n"
    <<FLAT<<"\t   Similar to option -u but applies for external branches if specified. If it's not specified then the minimum branch length of external\n"
    <<FLAT<<"\t   branches is set the same as the one of internal branch.\n"
    <<FLAT<<"\t" <<BOLD<<"-v " <<LINE<<"variance\n"
    <<FLAT<<"\t   Use this option to specify the way you want to apply variances for the branch lengths. Variances are used to recompense big errors on\n"
    <<FLAT<<"\t   long estimated branch lengths. The variance of the branch Bi is Vi = (Bi+b) where b is specified by option -b.\n"
    <<FLAT<<"\t   If " <<FLAT<<LINE<<"variance" <<FLAT<<"=0, then we don't use variance. If " <<FLAT<<LINE<<"variance" <<FLAT<<"=1, then LSD uses the input branch lengths to calculate variances.\n"
    <<FLAT<<"\t   If " <<FLAT<<LINE<<"variance" <<FLAT<<"=2, then LSD runs twice where the second time it calculates the variances based on the estimated branch\n"
    <<FLAT<<"\t   lengths of the first run. By default " <<FLAT<<LINE<<"variance" <<FLAT<<"=1.\n"
    <<FLAT<<"\t" <<BOLD<<"-V \n"
    <<FLAT<<"\t   Get the actual version.\n"
    <<FLAT<<"\t" <<BOLD<<"-w " <<LINE<<"givenRte\n"
    <<FLAT<<"\t   This option is used to specify the name of the file containing the substitution rates.\n"
    <<FLAT<<"\t   In this case, the program will use the given rates to estimate the dates of the nodes.\n"
    <<FLAT<<"\t   This file should have the following format\n"
    <<FLAT<<"\t        RATE1\n"
    <<FLAT<<"\t        RATE2\n"
    <<FLAT<<"\t        ...\n"
    <<FLAT<<"\t  where RATEi is the rate of the tree i in the inputTreesFile.\n"
    <<FLAT<<"\t" <<BOLD<<"-z " <<LINE<<"tipsDate\n"
    <<FLAT<<"\t   To specify the tips date if they are all equal. If the tips date is not a number, but a string (ex: 2020-01-10, or b(2019,2020))\n"
    <<FLAT<<"\t   then it should be put between the quotes.\n";
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
        cout<<"The file "<<outfile.c_str()<<" does not exist.\n";
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
        cout<<"The file "<<outfile.c_str()<<" does not exist.\n";
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

double getInputDate( string msg ,int& type)
{
    string word = getInputString( msg );
    double date;
    bool validDate = true;
    if( isReal(word.c_str()) ){
        date = atof( word.c_str() );
        type = 1;
        validDate = true;
    } else{
        validDate = readDateFromString(word.c_str(), date);
        type = 2;
    }
    if (!validDate) myErrorMsg("Your word is not recognized as a real or a valid date format year-month-day.\n");
    return date;
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

double getInputNonNegativeReal( string msg )
{
    string word;
    do
    {
        word = getInputString( msg );
        if( isReal(word.c_str()) && atof( word.c_str())>=0)
            break;
        myErrorMsg("Your word is not recognized as a non negative real.\n");
    } while( true );
    return atof( word.c_str() );
}

double getInputNonNegativeRealOrE( string msg )
{
    string word;
    do
    {
        word = getInputString( msg );
        if( isReal(word.c_str()) && atof( word.c_str())>=0)
            break;
        if (word.compare("e")==0)
            return -1;
        myErrorMsg("Your word is not recognized as a non negative real or letter e.\n");
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
        case 'd':
        case 'o':
        case 'p':
        case 's':
        case 'c':
        case 'v':
        case 'b':
        case 'r':
        case 'R':
        case 'g':
        case 'G':
        case 't':
        case 'w':
        case 'f':
        case 'n':
        case 'y':
        case 'q':
        case 'e':
        case 'm':
        case 'j':
        case 'u':
        case 'U':
        case 'l':
        case 'S':
        case 'D':
        case 'x':
        case 'a':
        case 'z':
        case 'h':
        case 'E':
            return true;
    }
    return false;
}

void setOptionsWithLetter( Pr* opt, char letter)
{
    //char* fnOut;
    switch( letter )
    {
        case 'x':
            exit( EXIT_SUCCESS );
        case 'i':
            opt->inFile = getInputFileName("Enter your Input File name> ");
            checkRooted(opt);
            break;
        case 'd':
            cout<<"Do you have a date file? y/n "<<endl;
            char letter[3];
            do{
                fgets(letter,3,stdin);
                if (*letter=='n' || *letter=='N') {
                    cout<<"There is no date file, so the program will estimate relative dates with root\n date = 0 and tips date = 1.\n Type 'y' to continue or 'n' to modify the root\n date and the tips date"<<endl;
                    char letter1[3];
                    do {
                        fgets( letter1, 3, stdin );
                        if (*letter1=='n' || *letter1=='N'){
                            opt->MRCA = getInputString("Enter the root date (default=0)> ");
                            opt->LEAVES = getInputString("Enter the tips date (default=1)> ");
                        }
                        else if (*letter1=='y' || *letter1=='Y'){
                            opt->MRCA="0";
                            opt->LEAVES="1";
                            opt->mrca=0;
                            opt->leaves=1;
                        }
                        else {
                            cout<<"Type 'y' to continue or 'n' to modify the root date and tips date"<<endl;
                        }
                    } while (*letter1!='n' && *letter1!='N' && *letter1!='y' && *letter1!='Y');
                }
                else if (*letter=='y' || *letter=='Y'){
                    opt->inDateFile = getInputFileName("Enter you input date file name>");
                }
            } while (*letter!='n' && *letter!='N' && *letter!='y' && *letter!='Y');
            break;
        case 'a':
            opt->MRCA = getInputString("Enter the root date > ");
            break;
        case 'z':
            opt->LEAVES = getInputString("Enter the tips date > ");
            break;
        case 'D':
            if (opt->outDateFormat==0) opt->outDateFormat=1;
            else if (opt->outDateFormat==1) opt->outDateFormat=2;
            else if (opt->outDateFormat==2) opt->outDateFormat=0;
            break;
        case 'p':
            opt->partitionFile = getInputFileName("Enter your Partition File name> ");
            break;
        case 's':
            opt->seqLength = getPositiveInputInteger("Enter your sequence length> ");
            break;
        case 'n':
            opt->nbData = getPositiveInputInteger("Enter your number of dataset> ");
            break;
        case 'o':
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
            opt->constraint=!opt->constraint;
            break;
        case 'v':
            if (opt->variance==0 || opt->variance==1) opt->variance+=1;
            else if (opt->variance==2) opt->variance=0;
            break;
        case 'b':
            if (opt->variance) opt->c = getInputPositiveReal("Enter the parameter for the variances> ");
            break;
        case 'r':
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
                    cout<<"The trees are not rooted, you must use either option -g to specify the outgroups\n file or -r to estimate the root"<<endl;
                }
                else if (opt->estimate_root.compare("as")==0){
                    opt->estimate_root="a";
                }
            }
            break;
        case 'g':
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
        case 'G':
            if (opt->fnOutgroup!=""){
                opt->removeOutgroup=!opt->removeOutgroup;
                if (!opt->removeOutgroup) {
                    opt->estimate_root="k";
                }
            }
            break;
        case 'w':
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
            opt->rho_min = getInputPositiveReal("Enter the lower bound for the rate> ");
            break;
        case 'f':
            if (!opt->ci) {
                opt->nbSampling = getInputInteger("Enter the number of sampling for calculating confidence intervals> ");
                opt->bootstraps_file = "";
                opt->ci=true;
            }
            else{
                if (opt->bootstraps_file==""){
                   opt->bootstraps_file = getInputFileName("Enter your bootstrap trees file> ");
                   opt->ci = true;
                } else{
                    opt->ci = false;
                }
            }
            break;
        case 'j':
            opt->verbose = !opt->verbose;
            break;
        case 'm':
            if (opt->e>0){
                opt->m = getInputInteger("Enter the number of sampling dated nodes for outliers detection> ");
            }
            break;
        case 'q':
            opt->q = getInputPositiveReal("Enter the standard deviation for lognormal relaxed clock to calculate confidence intervals> ");
            break;
        case 'e':
            opt->e = getInputReal("Enter the Zscore threshold for outliers detection> ");
            break;
        case 'u':
            opt->minblen = getInputNonNegativeRealOrE("Enter the minimum branch length of the time scaled tree> ");
            break;
        case 'U':
            opt->minblenL = getInputNonNegativeRealOrE("Enter the minimum external branch length of the time scaled tree> ");
            break;
        case 'l':
            opt->nullblen = getInputReal("Collapse all internal branches that are shorter than this value (default is 1/(2*seq_length))> ");
            break;
        case 'R':
            opt->round_time = getInputPositiveReal("Rounding number for minimum branch length of time scaled tree>");
            break;
        case 'S':
            opt->support = getInputNonNegativeReal("Threshold of support value to collapse>");
            break;
        case 'E':
            opt->splitExternal = !opt->splitExternal;
        case 'h':
            printHelp();
            break;
    }
}




