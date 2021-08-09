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
#include "utils.h"
#include "lsd.h"
using namespace lsd;

InputOutputStream::InputOutputStream () {
    inTree = nullptr;
    inOutgroup = nullptr;
    inDate = nullptr;
    inRate = nullptr;
    inPartition = nullptr;
    inBootstrapTree = nullptr;
    outResult = nullptr;
    outTree1 = nullptr;
    outTree2 = nullptr;
    outTree3 = nullptr;
}

InputOutputStream::InputOutputStream(string tree, string outgroup, string date, string rate, string partition, string bootstrap) {
    inTree = nullptr;
    inOutgroup = nullptr;
    inDate = nullptr;
    inRate = nullptr;
    inPartition = nullptr;
    inBootstrapTree = nullptr;
    setTree(tree);
    setOutgroup(outgroup);
    setDate(date);
    setRate(rate);
    setPartition(partition);
    setBootstrapTree(bootstrap);
    outResult = new ostringstream;
    outTree1 = new ostringstream;
    outTree2 = new ostringstream;
    outTree3 = new ostringstream;
}

InputOutputStream::~InputOutputStream() {
    if (inTree) {
        delete inTree;
        inTree = nullptr;
    }
    if (inOutgroup) {
        delete inOutgroup;
        inOutgroup = nullptr;
    }
    if (inDate) {
        delete inDate;
        inDate = nullptr;
    }
    if (inPartition) {
        delete inPartition;
        inPartition = nullptr;
    }
    if (inRate) {
        delete inRate;
        inRate = nullptr;
    }
    if (inBootstrapTree) {
        delete inBootstrapTree;
        inBootstrapTree = nullptr;
    }
    if (outResult) {
        delete outResult;
        outResult = nullptr;
    }
    if (outTree1) {
        delete outTree1;
        outTree1 = nullptr;
    }
    if (outTree2) {
        delete outTree2;
        outTree2 = nullptr;
    }
    if (outTree3) {
        delete outTree3;
        outTree3 = nullptr;
    }
}

void InputOutputStream::setTree(string str) {
    if (inTree)
        delete inTree;
    inTree = new istringstream(str);
}

void InputOutputStream::setOutgroup(string str) {
    if (str.empty())
        return;
    if (inOutgroup)
        delete inOutgroup;
    inOutgroup = new istringstream(str);
}

void InputOutputStream::setDate(string str) {
    if (str.empty())
        return;
    if (inDate)
        delete inDate;
    inDate = new istringstream(str);
}

void InputOutputStream::setPartition(string str) {
    if (str.empty())
        return;
    if (inPartition)
        delete inPartition;
    inPartition = new istringstream(str);
}

void InputOutputStream::setBootstrapTree(string str) {
    if (str.empty())
        return;
    if (inBootstrapTree)
        delete inBootstrapTree;
    inBootstrapTree = new istringstream(str);
}
void InputOutputStream::setRate(string str) {
    if (str.empty())
        return;
    if (inRate)
        delete inRate;
    inRate = new istringstream(str);
}

InputOutputFile::InputOutputFile(Pr *opt) : InputOutputStream() {
    treeIsFile = true;
    bootstrapTreeIsFile = true;
    // open the tree file
    ifstream *tree_file = new ifstream(opt->inFile);
    inTree = tree_file;
    if (!tree_file->is_open()){
        cerr << "Error: cannot open the input tree file " << opt->inFile << endl;
        exit(EXIT_FAILURE);
    }
    
    // open outgroup file
    if (opt->fnOutgroup!=""){
        ifstream *outgroup_file = new ifstream(opt->fnOutgroup);
        inOutgroup = outgroup_file;
        if (!outgroup_file->is_open()) {
            cerr << "Error: cannot open outgroup file " << opt->fnOutgroup << endl;
            exit(EXIT_FAILURE);
        }
    }
    
    // open date file
    if (opt->inDateFile != "") {
        ifstream *date_file = new ifstream(opt->inDateFile);
        inDate = date_file;
        if (!date_file->is_open()) {
            cerr << "Error: cannot open date file " << opt->inDateFile << endl;
            exit(EXIT_FAILURE);
        }
    }
    
    // open partition file
    if (opt->partitionFile != "") {
        ifstream *part_file = new ifstream(opt->partitionFile);
        inPartition = part_file;
        if (!part_file->is_open()) {
            cerr << "Error: cannot open partition file " << opt->partitionFile << endl;
            exit(EXIT_FAILURE);
        }
    }
    
    // open bootstrap tree file
    if (opt->bootstraps_file != "") {
        ifstream *bootstrap_file = new ifstream(opt->bootstraps_file);
        inBootstrapTree = bootstrap_file;
        if (!bootstrap_file->is_open()) {
            cerr << "Error: cannot open bootstrap file " << opt->bootstraps_file << endl;
            exit(EXIT_FAILURE);
        }
    }
    // open given rate file
    if (opt->rate != "") {
        ifstream *rate_file = new ifstream(opt->rate);
        inRate = rate_file;
        if (!rate_file->is_open()) {
            cerr << "Error: cannot open rate file " << opt->rate << endl;
            exit(EXIT_FAILURE);
        }
    }
    
    
    // open the result file
    ofstream *result_file = new ofstream(opt->outFile);
    outResult = result_file;
    if (!result_file->is_open()) {
        cerr << "Error: cannot create the output file "<< opt->outFile << endl;
        exit(EXIT_FAILURE);
    }
    
    ofstream *tree2_file = new ofstream(opt->treeFile2);
    outTree2 = tree2_file;
    if (!tree2_file->is_open()) {
        cerr << "Error: can not create the output tree file " << opt->treeFile2 << endl;
        exit(EXIT_FAILURE);
    }
    
    ofstream *tree3_file = new ofstream(opt->treeFile3);
    outTree3 = tree3_file;
    if (!tree3_file->is_open()) {
        cerr << "Error: can not create the output tree file " << opt->treeFile3 << endl;
        exit(EXIT_FAILURE);
    }
}

InputOutputFile::~InputOutputFile() {
    if (inTree && treeIsFile) {
        ((ifstream*)inTree)->close();
    }
    if (inOutgroup) {
        ((ifstream*)inOutgroup)->close();
    }
    if (inDate) {
        ((ifstream*)inDate)->close();
    }
    if (inPartition) {
        ((ifstream*)inPartition)->close();
    }
    if (inBootstrapTree && bootstrapTreeIsFile) {
        ((ifstream*)inBootstrapTree)->close();
    }
    if (inRate) {
        ((ifstream*)inRate)->close();
    }
    if (outResult) {
        ((ofstream*)outResult)->close();
    }
    if (outTree1) {
        ((ofstream*)outTree1)->close();
    }
    if (outTree2) {
        ((ofstream*)outTree2)->close();
    }
    if (outTree3) {
        ((ofstream*)outTree3)->close();
    }
}

void InputOutputFile::setTree(string str) {
    if (inTree) {
        if (treeIsFile) {
            ((ifstream*)inTree)->close();
        }
        delete inTree;
    }
    // change tree to stringstream
    treeIsFile = false;
    inTree = new istringstream(str);
}

void InputOutputFile::setBootstrapTree(string str) {
    if (inBootstrapTree) {
        if (bootstrapTreeIsFile) {
            ((ifstream*)inBootstrapTree)->close();
        }
        delete inBootstrapTree;
    }
    // change tree to stringstream
    bootstrapTreeIsFile = false;
    inBootstrapTree = new istringstream(str);
}

string readWord(istream& f,string fn){
    string s;
    char c=readChar(f,fn);
    int i=0;
    while (c<33 || c>126) c=readChar(f,fn);
    s=c;
    while (f.get(c) && c>=33 && c<=126 && c!='(' && c!=')' && c!=','){
        s=s+c;
        i++;
    }
    return s;
}

string readWord(string line,int& pos){
    string s;
    char c=line.at(pos);
    while (pos < line.length() && (c<33 || c>126)){ c=line.at(pos);pos++;}
    while (pos < line.length() && c>=33 && c<=126){
        c=line.at(pos);
        s=s+c;
        pos++;
    }
    return s;
}

int getLineNumber(istream &myfile){
    streampos pos = myfile.tellg();
    string line;
    int nb=0;
    while (getline(myfile, line)){
        nb++;
    }
    // go back to the current position in file
    myfile.clear();
    myfile.seekg(pos);
    return nb;
}

char readChar(istream& f,string fn){
    char r;
    if (f.get(r)) return r;
    else {
        cerr<<"Error in "<<fn<<endl;
        exit(EXIT_FAILURE);
    }
}

double readdouble(istream& f,string fn){
    double r;
    if (f >>r) return r;
    else {
        cerr<<"Error in "<<fn<<" : real expected"<<endl;
        exit(EXIT_FAILURE);
    }
}
string realToYearMonth(double year){
    ostringstream oss;
    double days = abs(year) - floor(abs(year));
    int y;
    if (year>0) {
        y = floor(year);
    } else {
        y = ceil(year);
    }
    int d = round(365*days);
    int m = 0;
    if (d>=1 && d<=31){
        m = 1;
    }
    if (d>=32 && d<=59){
        m = 2;
    }
    if (d>=60 && d<=90){
        m = 3;
    }
    if (d>=91 && d<=120){
        m = 4;
    }
    if (d>=121 && d<=151){
        m = 5;
    }
    if (d>=152 && d<=181){
        m = 6;
    }
    if (d>=182 && d<=212){
        m = 7;
    }
    if (d>=213 && d<=243){
        m = 8;
    }
    if (d>=244 && d<=273){
        m = 9;
    }
    if (d>=274 && d<=304){
        m = 10;
    }
    if (d>=305 && d<=334){
        m = 11;
    }
    if (d>=335){
        m = 12;
    }
    oss<<y;
    if (m!=0) {
        if (m<=9) oss<<"-0"<<m;
        else oss<<"-"<<m;
    }
    return oss.str();
}
string realToYearMonthDay(double year){
    ostringstream oss;
    double days = abs(year) - floor(abs(year));
    int y;
    if (year>0) {
        y = floor(year);
    } else {
        y = ceil(year);
    }
    int d = round(365*days);
    int m = 0;
    if (d>=1 && d<=31){
        m = 1;
    }
    if (d>=32 && d<=59){
        m = 2;
        d = d -31;
    }
    if (d>=60 && d<=90){
        m = 3;
        d = d - 59;
    }
    if (d>=91 && d<=120){
        m = 4;
        d = d - 90;
    }
    if (d>=121 && d<=151){
        m = 5;
        d = d - 120;
    }
    if (d>=152 && d<=181){
        m = 6;
        d = d - 151;
    }
    if (d>=182 && d<=212){
        m = 7;
        d = d - 181;
    }
    if (d>=213 && d<=243){
        m = 8;
        d = d -212;
    }
    if (d>=244 && d<=273){
        m = 9;
        d = d - 243;
    }
    if (d>=274 && d<=304){
        m = 10;
        d = d - 273;
    }
    if (d>=305 && d<=334){
        m = 11;
        d = d - 304;
    }
    if (d>=335){
        m = 12;
        d = d - 334;
    }
    oss<<y;
    if (d!=0 && m!=0) {
        if (m<=9) oss<<"-0"<<m;
        else oss<<"-"<<m;
        if (d<=9) oss<<"-0"<<d;
        else oss<<"-"<<d;
    }
    return oss.str();
}

double monthToReal(int m){
    return monthDayToReal(m,15);
}

int maxDate(int m){
    switch (m)
    {
        case 1:
            return 31;
            break;
        case 2:
            return 28;
            break;
        case 3:
            return 31;
            break;
        case 4:
            return 30;
            break;
        case 5:
            return 31;
            break;
        case 6:
            return 30;
            break;
        case 7:
            return 31;
            break;
        case 8:
            return 31;
            break;
        case 9:
            return 30;
            break;
        case 10:
            return 31;
            break;
        case 11:
            return 30;
            break;
        case 12:
            return 31;
            break;
    }
    cerr<<"Invalid month "<<m<<endl;
    exit(EXIT_FAILURE);
}
double monthDayToReal(int m,int d){
    switch (m)
    {
        case 1:
            if (d>=1 && d<=31) return (double)d/365.;
            break;
        case 2:
            if (d>=1 && d<=29) return (d+31.)/365.;
            break;
        case 3:
            if (d>=1 && d<=31) return (d+59.)/365.;
            break;
        case 4:
            if (d>=1 && d<=30) return (d+90.)/365.;
            break;
        case 5:
            if (d>=1 && d<=31) return (d+120.)/365.;
            break;
        case 6:
            if (d>=1 && d<=30) return (d+151.)/365.;
            break;
        case 7:
            if (d>=1 && d<=31) return (d+181.)/365.;
            break;
        case 8:
            if (d>=1 && d<=31) return (d+212.)/365.;
            break;
        case 9:
            if (d>=1 && d<=30) return (d+243.)/365.;
            break;
        case 10:
            if (d>=1 && d<=31) return (d+273.)/365.;
            break;
        case 11:
            if (d>=1 && d<=30) return (d+304.)/365.;
            break;
        case 12:
            if (d>=1 && d<=31) return (d+334.)/365.;
            break;
    }
    cerr<<"Invalid month-day "<<m<<"-"<<d<<endl;
    exit(EXIT_FAILURE);
}

double readDate(istream& f,string fn,Pr* pr,double& month,double& day){
    double y;
    month=-1;
    day=-1;
    int sign = 1;
    if (f >> y) {
        if (y<0) {
            sign = -1;
            y = -y;
        }
        char c = readChar(f,fn);
        if (c==')' || c==','){
            if (pr->inDateFormat!=2){
                if (y>=9 && y<=9999) pr->inDateFormat=1;
                else if (pr->inDateFormat!=1) pr->inDateFormat=0;
            }
            return y*sign;
        }
        else if (c=='-' && y==round(y)){
            int m;
            if (f >> m){
                month = m;
                if (pr->inDateFormat != 2) pr->inDateFormat=3;
                c = readChar(f,fn);
                if (c=='-'){
                    int d;
                    if (f >> d){
                        day = d;
                        c = readChar(f,fn);
                        if (c==')' || c==',') {
                            pr->inDateFormat=2;
                            return (y+monthDayToReal(m,d))*sign;
                        }
                    }
                } else if (c==')' || c==','){
                    return (y+monthToReal(m))*sign;
                }
            }
        }
    }
    cerr<<"Error reading input date : real or date format year-month-date or\n year-month expected"<<endl;
    exit(EXIT_FAILURE);
}

double readDate1(istream& f,string fn,char c,Pr* pr,double& month,double& day){
    month=-1;
    day=-1;
    string wd="";
    wd+=c;
    double y;
    int sign = 1;
    while (f.get(c) && c>=33 && c<=126 && c!='-') {
        wd+=c;
    }
    try {
        y=stod(wd.c_str());
    } catch (const std::invalid_argument&) {
        cerr<<"Error reading input date : real or date format year-month-date or\n year-month expected"<<wd<<endl;
        exit(EXIT_FAILURE);
    }
    if (c=='-' && y==round(y)){
        if (y<0){
            sign = -1;
            y = -y;
        }
            int m;
            if (f >> m){
                month = m;
                if (pr->inDateFormat != 2) pr->inDateFormat=3;
                c = readChar(f,fn);
                if (c=='-'){
                    int d;
                    if (f >> d){
                        day = d;
                        pr->inDateFormat=2;
                        return (y+monthDayToReal(m,d))*sign;
                    }
                } else{
                    f.unget();
                    return (y+monthToReal(m))*sign;
                }
            } 
    }
    else {
        if (pr->inDateFormat!=2){
            if (y>=9 && y<=9999) pr->inDateFormat=1;
            else if (pr->inDateFormat!=1) pr->inDateFormat=0;
        }
        return y;
    }
    cerr<<"Error reading input date : real or date format year-month-date or\n year-month expected"<<endl;
    exit(EXIT_FAILURE);
}

bool readDateFromString(const char* str,double& f){
    string y;
    while( *str!='\0' && *str>='0' && *str<='9'){
        y = y+(*str);
        str++;
    }
    if (*str=='-'){
        str++;
        string m;
        while( *str!='\0' && *str>='0' && *str<='9'){
            m = m+(*str);
            str++;
        }
        if(*str=='-'){
            str++;
            string d;
            while( *str!='\0' && *str>='0' && *str<='9'){
                d = d+(*str);
                str++;
            }
            f = (atoi(y.c_str())+monthDayToReal(atoi(m.c_str()),atoi(d.c_str())));
            return true;
        }
    }
    return false;
}

void readWholeDate(istream &dateFile,Pr* pr,int& type,double& v1,double& v2, double& m1,double& m2,double& d1,double& d2,int& dateFormat){
    char c = readChar(dateFile,"the  input date");
    while (c<33 || c>126) c=readChar(dateFile,"the input date");
    if (c=='l' || c=='L' || c=='u' || c=='U' || c=='b' || c=='B'){//interval value
        char p = readChar(dateFile,"the input date");
        if (p=='('){
            if (c=='l' || c=='L'){
                type='l';
                v1=readDate(dateFile,"the input date",pr,m1,d1);
                if (v1 == (int) v1){
                    if (m1<0 && dateFormat!=3) dateFormat = 1;
                    else if (d1<0) dateFormat = 3;
                }
            }
            else if (c=='u' || c=='U'){
                type='u';
                v1=readDate(dateFile,"the input date",pr,m1,d1);
                if (v1 == (int) v1){
                    if (m1<0 && dateFormat!=3) dateFormat = 1;
                    else if (d1<0) dateFormat = 3;
                }
            }
            else if (c=='b' || c=='B'){
                type='b';
                v1=readDate(dateFile,"the input date",pr,m1,d1);
                v2=readDate(dateFile,"the input date",pr,m2,d2);
                if (v1 == (int) v1){
                    if (m1<0 && dateFormat!=3) dateFormat = 1;
                    else if (d1<0) dateFormat = 3;
                }
                if (v2 == (int) v2){
                    if (m1<0 && dateFormat!=3) dateFormat = 1;
                    else if (d1<0) dateFormat = 3;
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
            while (c<33 || c>126) c=readChar(dateFile,"the input date");
        } else{
            cerr<<"Error reading inpute date: flexible temporal constraints must be defined\n as either 'l(lower_bound)' or 'u(upper_bound)' or 'b(lower_bound,upper_bound)'"<<endl;
            exit(EXIT_FAILURE);
        }
    }
    else {
        v1 = readDate1(dateFile,"the input date",c,pr,m1,d1);
        if (m1<0 && dateFormat!=3) dateFormat = 1;
        else if (d1<0) dateFormat = 3;
        type='p';
    }
}


void adjustNodeDateToYMD(Node*& node,int m1,int d1,int m2,int d2){
    if (node->type == 'p'){
        double year = node->D;
        int sign = 1;
        if (year < 0 ) {
            year = -year;
            sign = -1;
        }
        year = floor(year);
        if (m1<0){//only year
            node->type = 'b';
            node->lower = sign*(year + monthDayToReal(1,1));
            node->upper = sign*(year + monthDayToReal(12,31));
        } else if (d1<0){//only year-month
            node->type = 'b';
            node->lower = sign*(year + monthDayToReal(m1,1));
            node->upper = sign*(year + monthDayToReal(m1,maxDate(m1)));
        }
    }
    else if (node->type == 'l'){
        double year = node->lower;
        int sign = 1;
        if (year < 0 ) {
            year = -year;
            sign = -1;
        }
        year = floor(year);
        if (m1<0){
            node->lower = sign*(year + monthDayToReal(1,1));
        } else if (d1<0){
            node->lower = sign*(year + monthDayToReal(m1,1));
        }
    }
    else if (node->type == 'u'){
        double year = node->upper;
        int sign = 1;
        if (year < 0 ) {
            year = -year;
            sign = -1;
        }
        year = floor(year);
        if (m1<0){
            node->upper = sign*(year + monthDayToReal(12,31));
        } else if (d1<0){
            node->upper = sign*(year + monthDayToReal(m1,maxDate(m1)));
        }
    }
    else if (node->type == 'b'){
        double year = node->lower;
        int sign = 1;
        if (year < 0 ) {
            year = -year;
            sign = -1;
        }
        year = floor(year);
        if (m1<0){
            node->lower = sign*(year + monthDayToReal(1,1));
        } else if (d1<0){
            node->lower = sign*(year + monthDayToReal(m1,1));
        }
        
        year = node->upper;
        sign = 1;
        if (year < 0 ) {
            year = -year;
            sign = -1;
        }
        year = floor(year);
        if (m2<0){
            node->upper = sign*(year + monthDayToReal(12,31));
        } else if (d2<0){
            node->upper = sign*(year + monthDayToReal(m2,maxDate(m2)));
        }
    }
}

void adjustDateToYMD(Date*& date,int m1,int d1,int m2,int d2){
    if (date->type == 'p'){
        double year = date->date;
        int sign = 1;
        if (year < 0 ) {
            year = -year;
            sign = -1;
        }
        year = floor(year);
        if (m1<0){//only year
            date->type = 'b';
            date->lower = sign*(year + monthDayToReal(1,1));
            date->upper = sign*(year + monthDayToReal(12,31));
        } else if (d1<0){//only year-month
            date->type = 'b';
            date->lower = sign*(year + monthDayToReal(m1,1));
            date->upper = sign*(year + monthDayToReal(m1,maxDate(m1)));
        }
    }
    else if (date->type == 'l'){
        double year = date->lower;
        int sign = 1;
        if (year < 0 ) {
            year = -year;
            sign = -1;
        }
        year = floor(year);
        if (m1<0){
            date->lower = sign*(year + monthDayToReal(1,1));
        } else if (d1<0){
            date->lower = sign*(year + monthDayToReal(m1,1));
        }
    }
    else if (date->type == 'u'){
        double year = date->upper;
        int sign = 1;
        if (year < 0 ) {
            year = -year;
            sign = -1;
        }
        year = floor(year);
        if (m1<0){
            date->upper = sign*(year + monthDayToReal(12,31));
        } else if (d1<0){
            date->upper = sign*(year + monthDayToReal(m1,maxDate(m1)));
        }
    }
    else if (date->type == 'b'){
        double year = date->lower;
        int sign = 1;
        if (year < 0 ) {
            year = -year;
            sign = -1;
        }
        year = floor(year);
        if (m1<0){
            date->lower = sign*(year + monthDayToReal(1,1));
        } else if (d1<0){
            date->lower = sign*(year + monthDayToReal(m1,1));
        }
        
        year = date->upper;
        sign = 1;
        if (year < 0 ) {
            year = -year;
            sign = -1;
        }
        year = floor(year);
        if (m2<0){
            date->upper = sign*(year + monthDayToReal(12,31));
        } else if (d2<0){
            date->upper = sign*(year + monthDayToReal(m2,maxDate(m2)));
        }
    }
}

void adjustNodeDateToYM(Node*& node,int m1,int d1,int m2,int d2){
    if (node->type == 'p'){
        double year = node->D;
        int sign = 1;
        if (year < 0 ) {
            year = -year;
            sign = -1;
        }
        year = floor(year);
        if (m1<0){//only year
            node->type = 'b';
            node->lower = sign*(year + monthToReal(1));
            node->upper = sign*(year + monthToReal(12));
        }
    }
    else if (node->type == 'l'){
        double year = node->lower;
        int sign = 1;
        if (year < 0 ) {
            year = -year;
            sign = -1;
        }
        year = floor(year);
        if (m1<0){
            node->lower = sign*(year + monthToReal(1));
        }
    }
    else if (node->type == 'u'){
        double year = node->upper;
        int sign = 1;
        if (year < 0 ) {
            year = -year;
            sign = -1;
        }
        year = floor(year);
        if (m1<0){
            node->upper = sign*(year + monthToReal(12));
        }
    }
    else if (node->type == 'b'){
        double year = node->lower;
        int sign = 1;
        if (year < 0 ) {
            year = -year;
            sign = -1;
        }
        year = floor(year);
        if (m1<0){
            node->lower = sign*(year + monthToReal(1));
        }
        
        year = node->upper;
        sign = 1;
        if (year < 0 ) {
            year = -year;
            sign = -1;
        }
        year = floor(year);
        if (m2<0){
            node->upper = sign*(year + monthToReal(12));
        }
    }
}

void adjustDateToYM(Date*& date,int m1,int d1,int m2,int d2){
    if (date->type == 'p'){
        double year = date->date;
        int sign = 1;
        if (year < 0 ) {
            year = -year;
            sign = -1;
        }
        year = floor(year);
        if (m1<0){//only year
            date->type = 'b';
            date->lower = sign*(year + monthToReal(1));
            date->upper = sign*(year + monthToReal(12));
        }
    }
    else if (date->type == 'l'){
        double year = date->lower;
        int sign = 1;
        if (year < 0 ) {
            year = -year;
            sign = -1;
        }
        year = floor(year);
        if (m1<0){
            date->lower = sign*(year + monthToReal(1));
        }
    }
    else if (date->type == 'u'){
        double year = date->upper;
        int sign = 1;
        if (year < 0 ) {
            year = -year;
            sign = -1;
        }
        year = floor(year);
        if (m1<0){
            date->upper = sign*(year + monthToReal(12));
        }
    }
    else if (date->type == 'b'){
        double year = date->lower;
        int sign = 1;
        if (year < 0 ) {
            year = -year;
            sign = -1;
        }
        year = floor(year);
        if (m1<0){
            date->lower = sign*(year + monthToReal(1));
        }
        
        year = date->upper;
        sign = 1;
        if (year < 0 ) {
            year = -year;
            sign = -1;
        }
        year = floor(year);
        if (m2<0){
            date->upper = sign*(year + monthToReal(12));
        }
    }
}

vector<double> read_double_from_line(string line){
    stringstream ss(line);
    vector<double> results;
    double d;
    while (ss >> d){
        results.push_back(d);
    }
    return results;
}

double readDouble(string line,int& pos){
    string s;
    char c=line.at(pos);
    while (pos < line.length() && (c<33 || c>126)){c=line.at(pos);pos++;}
    while (pos < line.length() && c>=33 && c<=126){
        c=line.at(pos);
        s=s+c;
        pos++;
    }
    double result;
    stringstream ss(s);
    ss >> result;
    return result;
}

int readInt(istream& f,string msg){
    int c;
    if (f >> c) {
        return (int)(c);
    } else {
        cerr<<msg<<endl;
        exit(EXIT_FAILURE);
    }
}

string readSupport(istream& f,string fn){
    string s="";
    char c=readChar(f,fn);
    while (c!=':' && c!=';') {
        s=s+c;
        c=readChar(f,fn);
    }
    return s;
}

void concat(list<int> & l1,list<int> l2){
    for (list<int>::iterator iter = l2.begin();iter != l2.end();iter++){
        l1.push_back(*iter);
    }
}

void concatPos(list<int> l1,list<int> &l2){
    for (list<int>::iterator iter = l1.begin();iter != l1.end();iter++){
        l2.push_back(*iter);
    }
}

void concat(stack<int> & l1,list<int> l2){
    for (list<int>::iterator iter = l2.begin();iter != l2.end();iter++){
        l1.push(*iter);
    }
}

void concat(list<int> & l1,stack<int> l2){
    while (!l2.empty()) {
        int i=l2.top();
        l2.pop();
        l1.push_back(i);
    }
}

int getPosition(Node** nodes,string s,int n,int m){
    int i = n;
    int count=0;
    int k=0;
    while (i<m && count<2){
        if (nodes[i]->L.compare(s)==0){
            count++;
            k=i;
        }
        i++;
    }
    if (count==0) return -1;
    else if (count>1){
        cerr<<"There are at least two nodes that have the same label "<<s<<endl;
        exit(EXIT_FAILURE);
    }
    else return k;
}

bool isAncestor(Node** nodes,int i,int j){
    int x=j;
    while (x!=-1){
        if (x==i) return true;
        else x=nodes[x]->P;
    }
    return false;
}

int mrca(Node** nodes,int i,int j){
    int c = i;
    while (!isAncestor(nodes,c,j)){
        c=nodes[c]->P;
    }
    return c;
}

int index(list<int> & L, int e){
    int ind=0;
    for (list<int>::iterator i=L.begin();i!=L.end();i++){
        if (*i==e) return ind;
        ind++;
    }
    return -1;
}

int index(string s,string* & L,int n){
    for (int i=0;i<n;i++){
        if (s.compare(L[i])==0){
            return i;
        }
    }
    return -1;
}
bool contain(int s,list<int> l){
    for (list<int>::iterator iter = l.begin();iter!=l.end();iter++){
        if (s==*iter) return true;
    }
    return false;
}

bool contain(int s,vector<int> l){
    for (int i=0;i<l.size();i++){
        if (s==l[i]) return true;
    }
    return false;
}

bool contain(string s,list<string> l){
    for (list<string>::iterator iter = l.begin();iter!=l.end();iter++){
        if (s.compare(*iter)==0) return true;
    }
    return false;
}

list<int> intersect(list<int> l1,list<int> l2){
    list<int>::iterator it1=l1.begin();
    list<int>::iterator it2=l2.begin();
    list<int> common;
    while (it1!=l1.end() && it2!=l2.end() && *it1==*it2) {
        common.push_back(*it1);
        it1++;
        it2++;
    }
    return common;
}

string readLabel(char ch,istream& f,int& a){
    string s="";
    s=s+ch;
    char c=readChar(f,"input tree");
    while (c!=':' && c!=';' && c!=')'){
        s=s+c;
        f.get(c);
    }
    if (c==')') {
        while (c==')') {
            f.get(c);
            a--;
        }
    }
    return s.c_str();
}


char readBracket(istream& f,string fn){
    char c=readChar(f,fn);
    while (c!='('){
        c=readChar(f,fn);
    }
    return c;
}

char readCommaBracket(istream& f,string fn,string& s){
    char c=readChar(f,fn);
    s="";
    while (c==' ' || c=='	'){ c=readChar(f,fn);}
    while (c!=',' && c!=')'){
        s=s+c;
        c=readChar(f,fn);
    }
    return c;
}
char read2P(istream f,string fn){
    char c= readChar(f,fn);
    if (c==';') return c;
    while (c!=':' && c!=';') c=readChar(f,fn);
    return c;
}

void unrooted2rooted(Pr* &pr,Node** nodes){
    nodes[0] = new Node();
    nodes[0]->P=-1;
    int s = nodes[1]->suc[0];
    double br=nodes[s]->B;
    nodes[s]->B=br/2;
    nodes[1]->B=br/2;
    nodes[s]->P=0;
    nodes[1]->P=0;
    nodes[1]->suc.erase(nodes[1]->suc.begin());
    nodes[0]->suc.push_back(1);
    nodes[0]->suc.push_back(s);
    pr->rooted = true;
}

void rooted2unrooted(Pr* &pr,Node** nodes){
    int s1 = nodes[0]->suc[0];
    int s2 = nodes[0]->suc[1];
    if (s1 != 1){
        s2 = s1;
        s1 = 1;
    }
    nodes[s1]->P = -1;
    nodes[s2]->P = s1;
    nodes[s1]->suc.push_back(s2);
    nodes[s2]->B = nodes[s2]->B + nodes[s1]->B;
    pr->rooted = false;
}

/*Node** unrooted2rooted(Pr* & pr,Node** nodes){
    Node** nodes_new = new Node*[pr->nbBranches+1];
    for (int i=pr->nbINodes; i<=pr->nbBranches; i++) {
        nodes_new[i]=new Node();
        nodes_new[i]->P=nodes[i]->P;
        nodes_new[i]->B=nodes[i]->B;
        nodes_new[i]->L=nodes[i]->L;
        nodes_new[i]->D=nodes[i]->D;
        nodes_new[i]->type=nodes[i]->type;
        nodes_new[i]->lower=nodes[i]->lower;
        nodes_new[i]->upper=nodes[i]->upper;
        nodes_new[i]->status=nodes[i]->status;
    }
    nodes_new[0]=new Node();
    nodes_new[0]->P=-1;
    for (int i=1; i<pr->nbINodes; i++) {
        nodes_new[i] = new Node();
        nodes_new[i]->P=nodes[i]->P;
        nodes_new[i]->B=nodes[i]->B;
        nodes_new[i]->L=nodes[i]->L;
        nodes_new[i]->suc=nodes[i]->suc;
    }
    int s = nodes[1]->suc[0];
    double br=nodes[s]->B;
    nodes_new[s]->B=br/2;
    nodes_new[1]->B=br/2;
    nodes_new[s]->P=0;
    nodes_new[1]->P=0;
    nodes_new[1]->suc.erase(nodes_new[1]->suc.begin());
    nodes_new[0]->suc.push_back(1);
    nodes_new[0]->suc.push_back(s);
    for (int i=1;i<=pr->nbBranches;i++) delete nodes[i];
    delete[] nodes;
    pr->rooted = true;
    return nodes_new;
}*/

Node** unrooted2rootedS(Pr* &pr,Node** nodes,int s){//simplier version, use only for remove outgroup
    Node** nodes_new = new Node*[pr->nbBranches+1];
    for (int i=pr->nbINodes; i<=pr->nbBranches; i++) {
        nodes_new[i]=new Node();
        nodes_new[i]->P=nodes[i]->P;
        nodes_new[i]->B=nodes[i]->B;
        nodes_new[i]->L=nodes[i]->L;
    }
    for (int i=0; i<pr->nbINodes; i++) {
        nodes_new[i] = new Node();
        nodes_new[i]->P=nodes[i]->P;
        nodes_new[i]->B=nodes[i]->B;
        nodes_new[i]->L=nodes[i]->L;
    }
    nodes_new[0]=new Node();
    nodes_new[0]->P=-1;
    double br=nodes[s]->B;
    nodes_new[s]->B=br/2;
    nodes_new[1]->B=br/2;
    nodes_new[s]->P=0;
    nodes_new[1]->P=0;
    for (int i=0;i<=pr->nbBranches;i++) delete nodes[i];
    delete[] nodes;
    pr->rooted = true;
    computeSuc_polytomy(pr,nodes_new);
    return nodes_new;
}

void computeVariance(Pr* pr,Node** nodes){
    if (pr->variance==1 || pr->variance==2){
        for (int i=1;i<=pr->nbBranches;i++){
            nodes[i]->V=(nodes[i]->B+pr->b);
        }
    }
    else{
        for (int i=1;i<=pr->nbBranches;i++){
            nodes[i]->V=1./(double)(pr->seqLength);
        }
    }
}

void computeVarianceEstimateRoot(Pr* pr,Node** nodes,double br){
    if (pr->variance==1 || pr->variance==2){
        for (int i=1;i<=pr->nbBranches;i++){
            if (nodes[i]->P!=0) {
                nodes[i]->V=(nodes[i]->B+pr->b);
            }
            else nodes[i]->V=variance(pr,br);
        }
    }
    else{
        for (int i=1;i<=pr->nbBranches;i++){
            nodes[i]->V=1./(double)(pr->seqLength);
        }
    }
}

double variance(Pr* pr,double b){
    if (pr->variance==1 || pr->variance==2) return (b+pr->b);
    else return 1./(double)(pr->seqLength);
}

void computeNewVariance(Pr* pr,Node** nodes){
    if (pr->variance){
        for (int i=1;i<=pr->nbBranches;i++){
            if (nodes[i]->D>=nodes[nodes[i]->P]->D){
                nodes[i]->V=(pr->rho*nodes[i]->D-pr->rho*nodes[nodes[i]->P]->D+pr->b);
            }
            else{
                nodes[i]->V=(pr->b);
            }
        }
    } else {
        for (int i=1;i<=pr->nbBranches;i++){
            nodes[i]->V=1./(double)(pr->seqLength);
        }
    }
}

void computeNewVarianceEstimateRoot(Pr* pr,Node** nodes){
    double br=0;
    for (vector<int>::iterator iter=nodes[0]->suc.begin(); iter!=nodes[0]->suc.end(); iter++) {
        br+=(nodes[*iter]->D-nodes[0]->D)*pr->rho;
    }
    for (int i=1;i<=pr->nbBranches;i++){
        if (nodes[i]->P==0) {
            if (br>=0){
                nodes[i]->V=variance(pr,br);
            }
            else{
                nodes[i]->V=pr->b;
            }
        }
        else{
            if (nodes[i]->D>=nodes[nodes[i]->P]->D){
                nodes[i]->V=pr->rho*nodes[i]->D - pr->rho*nodes[nodes[i]->P]->D + pr->b;
            }
            else{
                nodes[i]->V=pr->b;
            }
        }
    }
}

void myExit( string msg, ... )
{
    va_list ptr;
    fprintf( stderr, "Error: " );
    va_start( ptr, msg );
    vfprintf( stderr, msg.c_str(), ptr );
    va_end( ptr );
    exit( EXIT_FAILURE );
}


void myErrorMsg( string msg, ... )
{
    va_list ptr;
    fprintf( stderr, "Error: " );
    va_start( ptr, msg );
    vfprintf( stderr, msg.c_str(), ptr );
    va_end( ptr );
    
}

bool isReal( const char* str )
{
    if( *str=='-' )
        str++;
    while( *str!='\0' && *str!='e' && *str!= 'E' && *str!='.')
    {
        if( !( ('0'<=*str && *str<='9')))
            return false;
        str++;
    }
    if (*str=='.'){
        str++;
        while( *str!='\0' && *str!='e' && *str!= 'E'){
            if( !( ('0'<=*str && *str<='9')))
               return false;
            str++;
        }
        if (*str=='e' || *str=='E'){
            str++;
            if (*str=='-' || *str=='+'){
                str++;
                while( *str!='\0'){
                    if( !( ('0'<=*str && *str<='9'))) return false;
                    str++;
                }
            } else {
                return false;
            }
        } else {
            return true;
        }
    } else if (*str=='e' || *str=='E'){
        str++;
        if (*str=='-' || *str=='+'){
            str++;
            while( *str!='\0'){
                if( !( ('0'<=*str && *str<='9'))) return false;
                str++;
            }
        } else {
            return false;
        }
    }
    return true;
}

bool isInteger( const char* str )
{
    if( *str=='-' )
        str++;
    while( *str!='\0' )
    {
        if( *str<'0' || '9'<*str )
            return false;
        str++;
    }
    return true;
}

void sort(int* & tab,int size){
    for (int i=0;i<size;i++){
        for (int j=i;j<size;j++){
            if (tab[i]>tab[j]){
                int temp = tab[i];
                tab[i]=tab[j];
                tab[j]=temp;
            }
        }
    }
}
void sort(double* & tab,int size){
    for (int i=0;i<size;i++){
        for (int j=i;j<size;j++){
            if (tab[i]>tab[j]){
                double temp = tab[i];
                tab[i]=tab[j];
                tab[j]=temp;
            }
        }
    }
}

double* sortTab(double* tab,int size){
    double* sortedTab = new double[size];
    for (int i=0;i<size;i++){
        sortedTab[i] = tab[i];
    }
    for (int i=0;i<size;i++){
        for (int j=i;j<size;j++){
            if (sortedTab[i]>sortedTab[j]){
                double temp = sortedTab[i];
                sortedTab[i]=sortedTab[j];
                sortedTab[j]=temp;
            }
        }
    }
    return sortedTab;
}

int index(int* & tab,int value,int size){
    for (int i=0;i<size;i++){
        if (value==tab[i]) return i;
    }
    return -1;
}

int mrca(Node** nodes,list<int> taxa){
    
    int t=taxa.front();
    
    taxa.pop_front();
    
    bool flag = false;
    
    while (!flag && nodes[t]->P!=-1){
        
        t = nodes[t]->P;
        
        flag=true;
        
        for (list<int>::iterator ia=taxa.begin();ia!=taxa.end();ia++){
            
            int j=*ia;
            
            if (!isAncestor(nodes,t,j)) {flag=false;break;}
            
        }
        
    }
    
    return t;
    
}
int mrca(Node** nodes,vector<int> taxa){
    int first=taxa.front();
    int t=first;
    bool flag = false;
    while (!flag && nodes[t]->P!=-1){
        t = nodes[t]->P;
        flag=true;
        for (vector<int>::iterator ia=taxa.begin();ia!=taxa.end();ia++){
            int j=*ia;
            if (j!=first && !isAncestor(nodes,t,j)) {flag=false;break;}
        }
    }
    return t;
}

void computeSuc(int* & Pre,int* & Suc1,int* & Suc2,int size,int n){
    for (int i=0;i<n;i++){
        Suc1[i]=-1;
    }
    for (int i=0;i<size;i++){
        if (Pre[i]!=-1){
            if (Suc1[Pre[i]]==-1) Suc1[Pre[i]]=i;
            else if (Suc1[Pre[i]]>i){
                Suc2[Pre[i]]=Suc1[Pre[i]];
                Suc1[Pre[i]]=i;
            }
            else Suc2[Pre[i]]=i;
        }
    }
}

bool markLeaf(Node* no){
    return no->status>=16;
}

bool limit(Node* no){
    return (no->status % 4 !=0);
}

void activeTC(Node* no){
    no->status+=4;
}

bool lower(Node* no){
    return (no->status %2 == 1);
}

bool upper(Node* no){
    return (no->status /2) %2 ==1 ;
}

bool tc(Node* no){
    return (no->status/4)%2==1;
}
bool leaf(Node* no){
    return ((no->status/8) % 2 == 1);
}

void activeMarkLeaf(Node* no){
    if (!markLeaf(no)) no->status+=16;
}

void desactiveMarkLeaf(Node* no){
    if (markLeaf(no)) no->status-=16;
}

void activeUpper(Node* no){
    no->status+=2+8;
}
void activeLower(Node* no){
    no->status+=1+8;
}
void desactive(Node* no){
    if (no->type=='p') {
        no->status=8;
    }
    else no->status=0;
}
void desactiveTC(Node* no){
    if (tc(no)) no->status-=4;
}
void desactiveLimit(Node* no){
    if (lower(no)) no->status-=1+8;
    if (upper(no)) no->status-=2+8;
}

bool initConstraintReRooted(Pr* pr,Node** nodes,int r,int p_r){
    bool constraintConsistent=true;
    for (int i=0;i<pr->nbINodes;i++){
        nodes[i]->type = 'n';
        nodes[i]->status = 0;
    }
    for (vector<Date*>::iterator iter=pr->internalConstraints.begin();iter!=pr->internalConstraints.end();iter++){
        Date* no = (*iter);
        int k=-1;
        if (no->mrca.size()==0){
            k=no->id;
        }
        else{
            k=mrca(nodes,(*iter)->mrca);
            no->id = k;
        }
        bool bl = (nodes[k]->addConstraint(*iter));
        constraintConsistent = constraintConsistent && bl;
    }
    if (constraintConsistent){
        for (int i=0; i<=pr->nbBranches; i++) {
            if (nodes[i]->type=='l' || nodes[i]->type=='b') {
                nodes[i]->D = nodes[i]->lower;
                activeLower(nodes[i]);
            }
            else if (nodes[i]->type=='u'){
                nodes[i]->D = nodes[i]->upper;
                activeUpper(nodes[i]);
            }
        }
        return checkAllConstraintConsistent(pr,nodes);
    } else {
        return false;
    }
}

bool checkAllConstraintConsistent(Pr* pr,Node** nodes){
    double* lowerX = new double[pr->nbBranches+1];
    bool* bl = new bool[pr->nbBranches+1];
    double* dates = new double[pr->nbBranches+1];
    double* lower = new double[pr->nbBranches+1];
    for (int i=0;i<=pr->nbBranches;i++){
        dates[i] = nodes[i]->D;
        lower[i] = nodes[i]->lower;
        if (nodes[i]->type=='l' || nodes[i]->type=='b'){
            bl[i]=true;
            lowerX[i]=nodes[i]->lower;
        }
        else if (nodes[i]->type=='p'){
            bl[i]=true;
            lowerX[i]=dates[i];
        }
        else {
            bl[i]=false;
        }
    }
    vector<int> pre = preorder_polytomy(pr,nodes);
    for (vector<int>::iterator iter = pre.begin();iter!=pre.end();iter++){
        int i=*iter;
        if (bl[i]){
            for (vector<int>::iterator iter=nodes[i]->suc.begin(); iter!=nodes[i]->suc.end(); iter++) {
                int s=*iter;
                if (!bl[s] || (nodes[s]->type=='l' && (nodes[s]->lower<lowerX[i]))){
                    lowerX[s]=lowerX[i];
                    lower[s]=lowerX[i];
                    bl[s]=true;
                }
                else if (nodes[s]->type=='u' && nodes[s]->upper<lowerX[i]){
                    delete[] lowerX;
                    delete[] bl;
                    delete[] dates;
                    delete[] lower;
                    return false;
                }
                else if (nodes[s]->type=='b'){
                    if (nodes[s]->upper>=lowerX[i]){
                        if (nodes[s]->lower<lowerX[i]){
                            lowerX[s]=lowerX[i];
                            lower[s]=lowerX[i];
                            bl[s]=true;
                        }
                    }
                    else {
                        delete[] lowerX;
                        delete[] bl;
                        delete[] dates;
                        delete[] lower;
                        return false;
                    }
                }
                else if ((nodes[s]->type=='p') && nodes[s]->D<lowerX[i]){
                    delete[] lowerX;
                    delete[] bl;
                    delete[] dates;
                    delete[] lower;
                    return false;
                }
                if (nodes[s]->D<lowerX[s]){
                    dates[s]=lowerX[s];
                }
            }
        }
    }
    delete[] lower;
    delete[] dates;
    delete[] lowerX;
    delete[] bl;
    return true;
}
bool initConstraint(Pr* pr,Node** nodes){
    bool constraintConsistent=true;
    for (int i=0;i<pr->nbINodes;i++){
        nodes[i]->type = 'n';
        nodes[i]->status = 0;
    }
    for (vector<Date*>::iterator iter=pr->internalConstraints.begin();iter!=pr->internalConstraints.end();iter++){
        Date* no = (*iter);
        int k=-1;
        if (no->mrca.size()==0){
            k=no->id;
        }
        else{
            k=mrca(nodes,no->mrca);
            no->id = k;
        }
        bool bl = (nodes[k]->addConstraint(*iter));
        constraintConsistent = constraintConsistent && bl;
    }
    if (constraintConsistent){
        for (int i=0; i<=pr->nbBranches; i++) {
            if (nodes[i]->type=='l' || nodes[i]->type=='b') {
                nodes[i]->D = nodes[i]->lower;
                activeLower(nodes[i]);
            }
            else if (nodes[i]->type=='u'){
                nodes[i]->D = nodes[i]->upper;
                activeUpper(nodes[i]);
            }
        }
        return checkAllConstraintConsistent(pr,nodes);
    } else {
        return false;
    }
}


Node** cloneLeaves(Pr* pr,Node** nodes,int f){
    Node** nodes_new =  new Node*[pr->nbBranches+1+f];
    for (int i=0; i<pr->nbINodes; i++) {
        nodes_new[i+f] = new Node();
    }
    for (int i=pr->nbINodes; i<=pr->nbBranches; i++) {
        nodes_new[i+f]=new Node();
        nodes_new[i+f]->P=nodes[i]->P+f;
        nodes_new[i+f]->B=nodes[i]->B;
        nodes_new[i+f]->L=nodes[i]->L;
        nodes_new[i+f]->V=nodes[i]->V;
        nodes_new[i+f]->type=nodes[i]->type;
        nodes_new[i+f]->lower=nodes[i]->lower;
        nodes_new[i+f]->upper=nodes[i]->upper;
        nodes_new[i+f]->D=nodes[i]->D;
        nodes_new[i+f]->status=nodes[i]->status;
        nodes_new[i+f]->minblen=nodes[i]->minblen;
    }
    return nodes_new;
}

void cloneInternalNodes(Pr* pr,Node** nodes,Node** &nodes_new,int f){
    for (int i=0; i<pr->nbINodes; i++) {
        nodes_new[i+f]->P=nodes[i]->P+f;
        nodes_new[i+f]->B=nodes[i]->B;
        nodes_new[i+f]->L=nodes[i]->L;
        nodes_new[i+f]->minblen=nodes[i]->minblen;
        nodes_new[i+f]->suc.clear();
        for (vector<int>::iterator iter=nodes[i]->suc.begin(); iter!=nodes[i]->suc.end(); iter++) {
            nodes_new[i+f]->suc.push_back((*iter)+f);
        }
    }
    for (int i=pr->nbINodes; i<=pr->nbBranches; i++) {
        nodes_new[i+f]->P=nodes[i]->P+f;
        nodes_new[i+f]->B=nodes[i]->B;
    }
}

/*bool outlierCheck(Pr* pr,Node** nodes){
    bool* out = new bool[pr->nbBranches+1];
    for (int i=0;i<=pr->nbBranches;i++) out[i]=false;
    for (int i=0;i<pr->outlier.size();i++){
        out[pr->outlier[i]] = true;
    }
    list<int> pos = postorder_polytomy(pr,nodes);
    for (list<int>::iterator iter=pos.begin();iter!=pos.end();iter++){
        int i = *iter;
        bool bl = true;
        for (vector<int>::iterator it = nodes[i]->suc.begin();it!=nodes[i]->suc.end();it++){
            bl = bl && out[*it];
        }
        out[i] = bl;
    }
    for (vector<int>::iterator it = nodes[0]->suc.begin();it!=nodes[0]->suc.end();it++){
        if (out[*it]==true){
            delete[] out;
            return false;
        }
    }
    delete[] out;
    return true;
}*/

bool reroot_rootedtree(double& br,int r,int s10,int s20,Pr* pr,Node** nodes,Node** &nodes_new){
    cloneInternalNodes(pr,nodes,nodes_new,0);
    if (r==s10 || r==s20){
        br = nodes[s10]->B+nodes[s20]->B;
        nodes_new[s10]->B=br;
        nodes_new[s20]->B=br;
        computeVarianceEstimateRoot(pr,nodes_new,br);
        return initConstraint(pr,nodes_new);
    }
    else {
        nodes_new[0]->L="";
        nodes_new[0]->P=-1;
        nodes_new[r]->P=0;
        nodes_new[nodes[r]->P]->P=0;
        nodes_new[0]->suc.clear();
        nodes_new[0]->suc.push_back(r);
        nodes_new[0]->suc.push_back(nodes[r]->P);
        int ii=r;
        int i=nodes[r]->P;
        int j=nodes[i]->P;
        while (j!=0){
            nodes_new[i]->suc.clear();
            nodes_new[i]->suc.push_back(j);
            for (vector<int>::iterator iter=nodes[i]->suc.begin(); iter!=nodes[i]->suc.end(); iter++) {
                if (*iter!=ii) {
                    nodes_new[i]->suc.push_back(*iter);
                }
            }
            nodes_new[j]->P=i;
            nodes_new[j]->B=nodes[i]->B;
            ii=i;
            i=j;
            j=nodes[i]->P;
        }
        int k=s10;
        if (k==i) k=s20;
        nodes_new[k]->P=i;
        nodes_new[i]->suc.clear();
        nodes_new[i]->suc.push_back(k);
        for (vector<int>::iterator iter=nodes[i]->suc.begin(); iter!=nodes[i]->suc.end(); iter++) {
            if (*iter!=ii) {
                nodes_new[i]->suc.push_back(*iter);
            }
        }
        br=nodes[r]->B;
        nodes_new[k]->B=nodes[i]->B+nodes[k]->B;
        nodes_new[r]->B=br;
        nodes_new[nodes[r]->P]->B=br;
        computeVarianceEstimateRoot(pr,nodes_new,br);
        return initConstraintReRooted(pr, nodes_new,k,i);
    }
}

bool reroot_rootedtree(double& br,int r,Pr* pr,Node** nodes){
    Node** nodes_new = cloneLeaves(pr,nodes,0);
    vector<int>::iterator iter=nodes[0]->suc.begin();
    int s10=(*iter);
    iter++;
    int s20=(*iter);
    for (int i=pr->nbINodes; i<=pr->nbBranches; i++) {
        nodes_new[i]->status=nodes[i]->status;
    }
    cloneInternalNodes(pr,nodes,nodes_new,0);
    if (r==s10 || r==s20){
        br = nodes[s10]->B+nodes[s20]->B;
        nodes_new[s10]->B=br;
        nodes_new[s20]->B=br;
        computeVarianceEstimateRoot(pr,nodes_new,br);
        nodes = nodes_new;
        return initConstraint(pr,nodes_new);
    }
    else {
        nodes_new[0]->L="";
        nodes_new[0]->P=-1;
        nodes_new[r]->P=0;
        nodes_new[nodes[r]->P]->P=0;
        nodes_new[0]->suc.clear();
        nodes_new[0]->suc.push_back(r);
        nodes_new[0]->suc.push_back(nodes[r]->P);
        int ii=r;
        int i=nodes[r]->P;
        int j=nodes[i]->P;
        while (j!=0){
            nodes_new[i]->suc.clear();
            nodes_new[i]->suc.push_back(j);
            for (vector<int>::iterator iter=nodes[i]->suc.begin(); iter!=nodes[i]->suc.end(); iter++) {
                if (*iter!=ii) {
                    nodes_new[i]->suc.push_back(*iter);
                }
            }
            nodes_new[j]->P=i;
            nodes_new[j]->B=nodes[i]->B;
            ii=i;
            i=j;
            j=nodes[i]->P;
        }
        int k=s10;
        if (k==i) k=s20;
        nodes_new[k]->P=i;
        nodes_new[i]->suc.clear();
        nodes_new[i]->suc.push_back(k);
        for (vector<int>::iterator iter=nodes[i]->suc.begin(); iter!=nodes[i]->suc.end(); iter++) {
            if (*iter!=ii) {
                nodes_new[i]->suc.push_back(*iter);
            }
        }
        br=nodes[r]->B;
        nodes_new[k]->B=nodes[i]->B+nodes[k]->B;
        nodes_new[r]->B=br;
        nodes_new[nodes[r]->P]->B=br;
        computeVarianceEstimateRoot(pr,nodes_new,br);
        nodes = nodes_new;
        return initConstraintReRooted(pr, nodes_new,k,i);
    }
}

bool reroot_rootedtree(double& br,int r,int s10,int s20,Pr* pr,Node** nodes,Node** &nodes_new,int* & P_ref,int* & tab){
    cloneInternalNodes(pr,nodes,nodes_new,0);
    for (int i=0; i<=pr->nbBranches; i++) {
        tab[i]=i;
        P_ref[i]=nodes[i]->P;
    }
    if (r==s10 || r==s20){
        for (int i=0;i<=pr->nbBranches;i++){
            nodes_new[i]->P=nodes[i]->P;
            P_ref[i]=nodes[i]->P;
            nodes_new[i]->B=nodes[i]->B;
        }
        br = nodes[s10]->B+nodes[s20]->B;
        nodes_new[s10]->B=br;
        nodes_new[s20]->B=br;
        computeVarianceEstimateRoot(pr,nodes_new,br);
        return initConstraint(pr, nodes_new);
    }
    else {
        nodes_new[0]->P=-1;
        P_ref[0]=-1;
        nodes_new[r]->P=0;
        P_ref[r]=0;
        nodes_new[nodes[r]->P]->P=0;
        P_ref[nodes[r]->P]=0;
        nodes_new[0]->suc.clear();
        nodes_new[0]->suc.push_back(r);
        nodes_new[0]->suc.push_back(nodes[r]->P);
        int ii=r;
        int i=nodes[r]->P;
        int j=nodes[i]->P;
        tab[i]=r;
        while (j!=0){
            nodes_new[i]->suc.clear();
            nodes_new[i]->suc.push_back(j);
            for (vector<int>::iterator iter=nodes[i]->suc.begin(); iter!=nodes[i]->suc.end(); iter++) {
                if (*iter!=ii) {
                    nodes_new[i]->suc.push_back(*iter);
                }
            }
            tab[j]=i;
            nodes_new[j]->P=i;
            P_ref[j]=i;
            nodes_new[j]->B=nodes[i]->B;
            ii=i;
            i=j;
            j=nodes[i]->P;
        }
        int k=s10;
        if (k==i) k=s20;
        nodes_new[k]->P=i;
        P_ref[k]=i;
        nodes_new[i]->suc.clear();
        nodes_new[i]->suc.push_back(k);
        for (vector<int>::iterator iter=nodes[i]->suc.begin(); iter!=nodes[i]->suc.end(); iter++) {
            if (*iter!=ii) {
                nodes_new[i]->suc.push_back(*iter);
            }
        }
        br=nodes[r]->B;
        nodes_new[k]->B=nodes[i]->B+nodes[k]->B;
        nodes_new[r]->B=br;
        nodes_new[nodes[r]->P]->B=br;
        computeVarianceEstimateRoot(pr,nodes_new,br);
        return initConstraintReRooted(pr, nodes_new,k,i);
    }
}

int reroot_rootedtree(int r,Pr* pr,Node** nodes,Node** & nodes_new){//used in extrait outgroups
    for (int i=0; i<pr->nbINodes; i++) {
        //nodes_new[i] = new Node();
        nodes_new[i]->P=nodes[i]->P;
        nodes_new[i]->B=nodes[i]->B;
        nodes_new[i]->L=nodes[i]->L;
    }
    vector<int>::iterator iter=nodes[0]->suc.begin();
    int s10=*iter;
    iter++;
    int s20=*iter;
    if (r==s10 || r==s20){
        for (int i=0;i<=pr->nbBranches;i++){
            nodes_new[i]->P=nodes[i]->P;
            nodes_new[i]->B=nodes[i]->B;
        }
        double br = nodes[s10]->B+nodes[s20]->B;
        nodes_new[s10]->B=br/2.;
        nodes_new[s20]->B=br/2.;
        if (r==s10) {
            return s20;
        }
        else return s10;
    }
    else {
        nodes_new[0]->L="";
        nodes_new[0]->P=-1;
        nodes_new[r]->P=0;
        nodes_new[nodes[r]->P]->P=0;
        int i=nodes[r]->P;
        int j=nodes[i]->P;
        while (j!=0){
            nodes_new[j]->P=i;
            nodes_new[j]->B=nodes[i]->B;
            i=j;
            j=nodes[i]->P;
        }
        int k=s10;
        if (k==i) k=s20;
        nodes_new[k]->P=i;
        double br=nodes[r]->B;
        nodes_new[k]->B=nodes[i]->B+nodes[k]->B;
        nodes_new[r]->B=br/2.;
        nodes_new[nodes[r]->P]->B=br/2.;
        return nodes[r]->P;
    }
}

void computeObjective(Pr* pr,Node** nodes){
    pr->objective = 0;
    for (int i=1;i<=pr->nbBranches;i++){
        //p+=(B[i]-rho*T[i]+rho*T[P[i]])*(B[i]-rho*T[i]+rho*T[P[i]])/(2*V[i]) +log(2*M_PI*V[i])/2;
        pr->objective+=(nodes[i]->B-pr->rho*nodes[i]->D+pr->rho*nodes[nodes[i]->P]->D)*(nodes[i]->B-pr->rho*nodes[i]->D+pr->rho*nodes[nodes[i]->P]->D)/(nodes[i]->V);
    }//p = -(log likelihood)
}

void computeObjectiveMultiRates(Pr* pr,Node** nodes){
    pr->objective = 0;
    for (int i=1;i<=pr->nbBranches;i++){
        double rate = pr->rho*pr->multiplierRate[nodes[i]->rateGroup];
        pr->objective+=(nodes[i]->B-rate*nodes[i]->D+rate*nodes[nodes[i]->P]->D)*(nodes[i]->B-rate*nodes[i]->D+rate*nodes[nodes[i]->P]->D)/(nodes[i]->V);
    }//p = -(log likelihood)
}

void computeObjectiveMultiRates(Pr* pr,Node** nodes,double* B,double* V){
    pr->objective = 0;
    double test = 0;
    for (int i=1;i<=pr->nbBranches;i++){
        double rate = pr->rho*pr->multiplierRate[nodes[i]->rateGroup];
        pr->objective+=(B[i]-rate*nodes[i]->D+rate*nodes[nodes[i]->P]->D)*(B[i]-rate*nodes[i]->D+rate*nodes[nodes[i]->P]->D)/(V[i]);
    }//p = -(log likelihood)
}

void computeObjectiveEstimateRoot(int r,int p_r,double br,Pr* pr,Node** nodes){
    pr->objective = (br-pr->rho*nodes[r]->D-pr->rho*nodes[p_r]->D+2*pr->rho*nodes[0]->D)*(br-pr->rho*nodes[r]->D-pr->rho*nodes[p_r]->D+2*pr->rho*nodes[0]->D)/nodes[r]->V;
    for (int i=1;i<=pr->nbBranches;i++){
        //p+=(B[i]-rho*T[i]+rho*T[P[i]])*(B[i]-rho*T[i]+rho*T[P[i]])/(2*V[i]) +log(2*M_PI*V[i])/2;
        if (i!=r && i!=p_r) {
            pr->objective+=(nodes[i]->B-pr->rho*nodes[i]->D+pr->rho*nodes[nodes[i]->P]->D)*(nodes[i]->B-pr->rho*nodes[i]->D+pr->rho*nodes[nodes[i]->P]->D)/(nodes[i]->V);
        }
    }//p = -(log likelihood)
}


string newick(int i,int terminate,Pr* pr,Node** nodes,int& nbTips){;
    ostringstream b;
    if (i>0){
        b<<nodes[i]->B;
    }
    if (i>=pr->nbINodes){
        nbTips++;
        return nodes[i]->L+":"+b.str();
    }
    else{
        string newLabel="(";
        for (vector<int>::iterator iter=nodes[i]->suc.begin(); iter!=nodes[i]->suc.end(); iter++) {
            int s = *iter;
            string l=newick(s,terminate,pr,nodes,nbTips);
            if (iter==nodes[i]->suc.begin()) newLabel+=l;
            else newLabel+=","+l;
        }
        if (i!=terminate) {
            return newLabel+")"+nodes[i]->L+":"+b.str();
        }
        else{
            return newLabel+")"+nodes[i]->L+";\n";
        }
    }
}

string nexus(int i,Pr* pr,Node** nodes){
    ostringstream b,date;
    if (i>0){
        b<<nodes[i]->B;
    }
    if (pr->outDateFormat==2){
        date<<realToYearMonthDay(nodes[i]->D);
    } else if (pr->outDateFormat==3){
        date<<realToYearMonth(nodes[i]->D);
    } else{
        date<<nodes[i]->D;
    }
    if (i>=pr->nbINodes) return nodes[i]->L+"[&date=\""+date.str()+"\"]:"+b.str();
    else{
        string newLabel="(";
        for (vector<int>::iterator iter=nodes[i]->suc.begin(); iter!=nodes[i]->suc.end(); iter++) {
            int s = *iter;
            string l=nexus(s,pr,nodes);
            if (iter==nodes[i]->suc.begin()) newLabel+=l;
            else newLabel+=","+l;
        }
        if (i>0) {
            return newLabel+")"+nodes[i]->L+"[&date=\""+date.str()+"\"]:"+b.str();
            //if (abs(nodes[i]->B)>0) return newLabel+")"+nodes[i]->L+"[&date=\""+date.str()+"\"]:"+b.str();
            //else return newLabel+")"+nodes[i]->L+":"+b.str();
        }
        else{
            return newLabel+")"+nodes[i]->L+"[&date=\""+date.str()+"\"];\n";
        }
    }
}

string nexusDate(int i,Pr* pr,Node** nodes){
    ostringstream b,date;
    if (i>0){
        b<< (nodes[i]->D - nodes[nodes[i]->P]->D);
    }
    if (pr->outDateFormat==2){
        date<<realToYearMonthDay(nodes[i]->D);
    } else if (pr->outDateFormat==3){
        date<<realToYearMonth(nodes[i]->D);
    } else{
        date<<nodes[i]->D;
    }
    if (i>=pr->nbINodes) return nodes[i]->L+"[&date=\""+date.str()+"\"]:"+b.str();
    else{
        string newLabel="(";
        for (vector<int>::iterator iter=nodes[i]->suc.begin(); iter!=nodes[i]->suc.end(); iter++) {
            int s = *iter;
            string l=nexusDate(s,pr,nodes);
            if (iter==nodes[i]->suc.begin()) newLabel+=l;
            else newLabel+=","+l;
        }
        if (i>0) {
            return newLabel+")"+nodes[i]->L+"[&date=\""+date.str()+"\"]:"+b.str();
            //if (abs(nodes[i]->B)>0) return newLabel+")"+nodes[i]->L+"[&date=\""+date.str()+"\"]:"+b.str();
            //else return newLabel+")"+nodes[i]->L+":"+b.str();
        }
        else{
            return newLabel+")"+nodes[i]->L+"[&date=\""+date.str()+"\"];\n";
        }
    }
}

string nexusIC(int i,Pr* pr,Node** nodes,double* D_min,double* D_max,double* H_min,double* H_max){
    ostringstream b,date;
    if (i>0){
        b<<nodes[i]->B;
    }
    if (pr->outDateFormat==2){
        date<<realToYearMonthDay(nodes[i]->D);
    } else if (pr->outDateFormat==3){
        date<<realToYearMonth(nodes[i]->D);
    } else{
        date<<nodes[i]->D;
    }
    if (i>=pr->nbINodes && nodes[i]->type=='p') return nodes[i]->L+"[&date="+date.str()+"]:"+b.str();
    else{
        ostringstream hmin,hmax,dmin,dmax;
        if (pr->outDateFormat==2){
            dmin<< realToYearMonthDay(D_min[i]);
            dmax<< realToYearMonthDay(D_max[i]);
        } else if (pr->outDateFormat==3){
            dmin<< realToYearMonth(D_min[i]);
            dmax<< realToYearMonth(D_max[i]);
        } else {
            dmin<< D_min[i];
            dmax<< D_max[i];
        }
        hmin<< H_min[i];
        hmax<< H_max[i];
        if (i>=pr->nbINodes) {
            return nodes[i]->L+"[&date=\""+date.str()+"\",CI_height={"+hmin.str()+","+hmax.str()+"},CI_date=\"{"+dmin.str()+","+dmax.str()+"}\"]:"+b.str();
        }
        else{
            string newLabel="(";
            for (vector<int>::iterator iter=nodes[i]->suc.begin(); iter!=nodes[i]->suc.end(); iter++) {
                int s = *iter;
                string l=nexusIC(s,pr,nodes,D_min,D_max,H_min,H_max);
                if (iter==nodes[i]->suc.begin()) newLabel+=l;
                else newLabel+=","+l;
            }
            if (i>0) {
                return newLabel+")"+nodes[i]->L+"[&date=\""+date.str()+"\",CI_height={"+hmin.str()+","+hmax.str()+"},CI_date=\"{"+dmin.str()+","+dmax.str()+"}\"]:"+b.str();
                //if (abs(nodes[i]->B)>0) return newLabel+")"+nodes[i]->L+"[&date=\""+date.str()+"\",CI_height=\"{"+hmin.str()+","+hmax.str()+"}\",CI_date=\"{"+dmin.str()+","+dmax.str()+"}\"]:"+b.str();
                //else return newLabel+")"+nodes[i]->L+":"+b.str();
            }
            else{
                return newLabel+")"+nodes[i]->L+"[&date=\""+date.str()+"\",CI_height={"+hmin.str()+","+hmax.str()+"},CI_date=\"{"+dmin.str()+","+dmax.str()+"}\"];\n";
            }
            
        }
    }
}

string nexusICDate(int i,Pr* pr,Node** nodes,double* D_min,double* D_max,double* H_min,double* H_max){
    ostringstream b,date;
    if (i>0){
        b<< (nodes[i]->D - nodes[nodes[i]->P]->D) ;
    }
    if (pr->outDateFormat==2){
        date<<realToYearMonthDay(nodes[i]->D);
    } else if (pr->outDateFormat==3){
        date<<realToYearMonth(nodes[i]->D);
    } else{
        date<<nodes[i]->D;
    }
    if (i>=pr->nbINodes && nodes[i]->type=='p') return nodes[i]->L+"[&date="+date.str()+"]:"+b.str();
    else{
        ostringstream hmin,hmax,dmin,dmax;
        if (pr->outDateFormat==2){
            dmin<< realToYearMonthDay(D_min[i]);
            dmax<< realToYearMonthDay(D_max[i]);
        } else if (pr->outDateFormat==3){
            dmin<< realToYearMonth(D_min[i]);
            dmax<< realToYearMonth(D_max[i]);
        } else{
            dmin<< D_min[i];
            dmax<< D_max[i];
        }
        hmin<< H_min[i];
        hmax<< H_max[i];
        if (i>=pr->nbINodes) {
            return nodes[i]->L+"[&date=\""+date.str()+"\",CI_height={"+hmin.str()+","+hmax.str()+"},CI_date={"+dmin.str()+","+dmax.str()+"}]:"+b.str();
        }
        else{
            string newLabel="(";
            for (vector<int>::iterator iter=nodes[i]->suc.begin(); iter!=nodes[i]->suc.end(); iter++) {
                int s = *iter;
                string l=nexusICDate(s,pr,nodes,D_min,D_max,H_min,H_max);
                if (iter==nodes[i]->suc.begin()) newLabel+=l;
                else newLabel+=","+l;
            }
            if (i>0) {
                return newLabel+")"+nodes[i]->L+"[&date=\""+date.str()+"\",CI_height={"+hmin.str()+","+hmax.str()+"},CI_date=\"{"+dmin.str()+","+dmax.str()+"}\"]:"+b.str();
                /*if (abs(nodes[i]->B)>0) return newLabel+")"+nodes[i]->L+"[&date=\""+date.str()+"\",CI_height={"+hmin.str()+","+hmax.str()+"},CI_date=\"{"+dmin.str()+","+dmax.str()+"}\"]:"+b.str();
                 else return newLabel+")"+nodes[i]->L+":"+b.str();*/
            }
            else{
                return newLabel+")"+nodes[i]->L+"[&date=\""+date.str()+"\",CI_height={"+hmin.str()+","+hmax.str()+"},CI_date=\"{"+dmin.str()+","+dmax.str()+"}\"];\n";
            }
            
        }
    }
}

double* variance(int w,int m,double* B,int c,int s){
    double* V = new double[m+1];
    if (w==1 || w==2){
        for (int i=1;i<=m;i++){
            //V[i]=(B[i]+c/(c+s))/(c+s);
            V[i]=(B[i]+(double)c/s)/s;
        }
    }
    else{
        for (int i=1;i<=m;i++){
            V[i]=1./(double)(s);
        }
    }
    return V;
}

double* newvariance(int m,double rho,int* P,double* T,int c,int s){
    double* V = new double[m+1];
    for (int i=1;i<=m;i++){
        if (T[i]>=T[P[i]]) V[i]=(rho*T[i]-rho*T[P[i]]+(double)c/s)/s;
        else V[i]=((double)c/s)/s;
    }
    return V;
}

list<string> getOutgroup(istream &o, string fn){
    list<string> result;
    int lineNb=getLineNumber(o);
    int ino=readInt(o,"Error in the outgroup file, the file should begin with an integer (the number of outgroups)");
    if (lineNb-1<ino) {
        cout<<"The number of given outgroups is small than the number of outgroups to read.\n Please change the number of outgroups to read at the first line of the outgroup\n file."<<endl;
        exit(EXIT_FAILURE);
    }
    for (int i=0;i<ino;i++) result.push_back(readWord(o,fn));
    return result;
}

void initialize_status(Pr* &pr,Node** &nodes){
    for (int i=0; i<=pr->nbBranches; i++) {
        if (nodes[i]->type=='p') nodes[i]->status=8;
        else nodes[i]->status=0;
    }
    /*for (int i=0;i<pr->outlier.size();i++){
        nodes[pr->outlier[i]]->status=0;
    }*/
}

list<int> getActiveSet(Pr* pr,Node** nodes){
    list<int> active_set;
    for (int i=0;i<=pr->nbBranches;i++){
        if (tc(nodes[i])) active_set.push_back(i);
        if (limit(nodes[i])) active_set.push_back(-i);
    }
    return active_set;
}

/*void computeSuc_polytomy(Pr* pr,Node** nodes){
    for (int i=0;i<pr->nbINodes;i++){
        nodes[i]->suc.clear();
    }
    for (int i=1;i<=pr->nbBranches;i++){
        nodes[nodes[i]->P]->suc.push_back(i);
    }
}*/

void computeSuc_polytomy(Pr* pr,Node** nodes){
    int root = (int)(!pr->rooted);
    for (int i=root;i<pr->nbINodes;i++){
        nodes[i]->suc.clear();
    }
    bool* visited = new bool[pr->nbINodes];
    for (int i=root; i<pr->nbINodes;i++) visited[i] = false;
    int v = pr->nbBranches;
    while (v!= (pr->nbINodes-1)){
        int vt = v;
        int p = nodes[v]->P;
        nodes[p]->suc.push_back(v);
        while (p!=root && !visited[p]){
            visited[p] = true;
            v = p;
            p = nodes[v]->P;
            nodes[p]->suc.push_back(v);
        }
        v = vt-1; 
    }
    delete[] visited;
}

list<int> pos_polytomy(int i,Pr* pr,Node** nodes){
    list<int> l;
    if (i>=pr->nbINodes) return l;
    else{
        for (vector<int>::iterator iter = nodes[i]->suc.begin();iter!=nodes[i]->suc.end();iter++){
            if (*iter<pr->nbINodes) {
                list<int> l1 = pos_polytomy(*iter,pr,nodes);
                for (list<int>::iterator iter1=l1.begin();iter1!=l1.end();iter1++){
                    l.push_back(*iter1);
                }
            }
        }
        l.push_back(i);
        return l;
    }
}

list<int> postorder_polytomy(Pr* pr,Node** nodes){
    int root=0;
    for (root=0;root<pr->nbINodes;root++){
        if (nodes[root]->P==-1) break;
    }
    return pos_polytomy(root,pr,nodes);
}

vector<int> pre_polytomy(int i,Pr* pr,Node** nodes){
    vector<int> l;
    if (i>=pr->nbINodes) return l;
    else{
        l.push_back(i);
        for (vector<int>::iterator iter = nodes[i]->suc.begin();iter!=nodes[i]->suc.end();iter++){
            if (*iter<pr->nbINodes) {
                vector<int> l1 = pre_polytomy(*iter,pr,nodes);
                for (vector<int>::iterator iter1=l1.begin();iter1!=l1.end();iter1++){
                    l.push_back(*iter1);
                }
            }
        }
        return l;
    }
}

vector<int> preorder_polytomy(Pr* pr,Node** nodes){
    int root=0;
    for (root=0;root<pr->nbINodes;root++){
        if (nodes[root]->P==-1) break;
    }
    return pre_polytomy(root,pr,nodes);
}

vector<int> pre_polytomy_withTips(int i,Pr* pr,Node** nodes){
    vector<int> l;
    l.push_back(i);
    if (i>=pr->nbINodes){
        return l;
    }
    else{
        for (vector<int>::iterator iter = nodes[i]->suc.begin();iter!=nodes[i]->suc.end();iter++){
            vector<int> l1 = pre_polytomy_withTips(*iter,pr,nodes);
            for (vector<int>::iterator iter1=l1.begin();iter1!=l1.end();iter1++){
                l.push_back(*iter1);
            }
        }
        return l;
    }
}

vector<int> preorder_polytomy_withTips(Pr* pr,Node** nodes){
    int root=0;
    for (root=0;root<=pr->nbBranches;root++){
        if (nodes[root]->P==-1) break;
    }
    return pre_polytomy_withTips(root,pr,nodes);
}

list<int> down_polytomy(int i,Pr* pr,Node** nodes){
    list<int> result;
    result.push_back(i);
    activeMarkLeaf(nodes[i]);
    nodes[i]->D = nodes[nodes[i]->P]->D + nodes[i]->minblen;
    if (i<pr->nbINodes){
        for (vector<int>::iterator iter=nodes[i]->suc.begin();iter!=nodes[i]->suc.end();iter++){
            if (tc(nodes[*iter])){
                list<int> l1=down_polytomy(*iter,pr,nodes);
                concat(result,l1);
            }
        }
    }
    return result;
}

stack<int>* computeFeuilles_polytomy(list<int> ls,Pr* pr,Node** nodes){
    stack<int>* feuilles = new stack<int>[ls.size()];
    int count=0;
    for (list<int>::iterator iter=ls.begin();iter!=ls.end();iter++){
        int i=*iter;
        bool alone=true;
        list<int> ai;
        for (vector<int>::iterator iter=nodes[i]->suc.begin();iter!=nodes[i]->suc.end();iter++){
            if (tc(nodes[*iter])){
                alone=false;
                ai.push_back(*iter);
            }
        }
        if ((!tc(nodes[i])) && (i>=pr->nbINodes || alone)){//nodes[i] is alone
            feuilles[count].push(i);
        }
        else {
            if (!tc(nodes[i])){
                feuilles[count].push(i);
            }
            else{
                int j=i;
                while (j!=-1 && tc(nodes[j])){
                    feuilles[count].push(j);
                    int k=j;
                    j=nodes[j]->P;
                    if (j!=-1 && tc(nodes[k])){
                        activeMarkLeaf(nodes[j]);//leaf
                        nodes[j]->D = nodes[k]->D - nodes[k]->minblen;
                        for (vector<int>::iterator iter=nodes[j]->suc.begin(); iter!=nodes[j]->suc.end(); iter++) {
                            if ((*iter)!=k && tc(nodes[*iter])) {
                                ai.push_back(*iter);
                            }
                        }
                    }
                }
            }
            for (list<int>::iterator it=ai.begin();it!=ai.end();it++){
                concat(feuilles[count],down_polytomy(*it,pr,nodes));
            }
        }
        count++;
    }
    return feuilles;
}


list<int> suc_polytomy(int i,int j,Pr* pr,Node** nodes,int* & Pre,double* & add,list<int> &suc){
    list<int> result;
    if (leaf(nodes[i]) && i!=j) {
        nodes[j]->D = nodes[nodes[j]->P]->D + nodes[j]->minblen;
        activeMarkLeaf(nodes[j]);
    }
    if (j >= pr->nbINodes) {
        result.push_back(j);
        Pre[j]=i;
        if (markLeaf(nodes[j])) suc.push_back(j);
    }
    else{
        for (vector<int>::iterator iter=nodes[j]->suc.begin(); iter!=nodes[j]->suc.end(); iter++) {
            int s=*iter;
            if (!tc(nodes[s])) {
                if (i!=j) add[s] = add[j];
                else add[s] = 0;
                if (markLeaf(nodes[s]) || s<pr->nbINodes) {
                    suc.push_back(s);
                }
                Pre[s]=i;
            }
            else{
                if (i!=j) add[s] = nodes[s]->minblen + add[j];
                else add[s] = nodes[s]->minblen;
                list<int> l1=suc_polytomy(i,s,pr,nodes,Pre,add,suc);
                concatPos(l1,result);
            }
        }
        if (j!=i) result.push_back(j);
    }
    return result;
}

void reduceTree_polytomy(Pr* pr,Node** nodes,int* &Pre,list<int>* & Suc,double* &add,list<int>* &internal){
    int count=0;
    Pre[0]=-1;
    for (int i=0;i<pr->nbINodes;i++){
        if ((!tc(nodes[i])) && (!markLeaf(nodes[i]) || nodes[i]->type=='p')){
            bool bl=false;
            for (vector<int>::iterator iter=nodes[i]->suc.begin(); iter!=nodes[i]->suc.end(); iter++) {
                if (tc(nodes[*iter])) {
                    bl=true;
                    break;
                }
            }
            if (bl){
                internal[count]=suc_polytomy(i,i,pr,nodes,Pre,add,Suc[i]);
                count++;
            }
            else {
                for (vector<int>::iterator iter=nodes[i]->suc.begin(); iter!=nodes[i]->suc.end(); iter++) {
                    Pre[*iter]=i;
                    if (markLeaf(nodes[*iter]) || *iter<pr->nbINodes) Suc[i].push_back(*iter);
                }
            }
        }
        else {
            for (vector<int>::iterator iter=nodes[i]->suc.begin(); iter!=nodes[i]->suc.end(); iter++) {
                if (Pre[*iter]==-1) Pre[*iter]=i;
                if (markLeaf(nodes[*iter]) || *iter<pr->nbINodes) Suc[i].push_back(*iter);
            }
        }
    }
}

void checkRooted(Pr* &opt){
    ifstream tree;
    tree.open(opt->inFile.c_str());
    if (!tree.is_open()){
        cerr<<"Error: can not open the input file"<<endl;
        exit(EXIT_FAILURE);
    }
    stack<int> pileNode;
    char c = readBracket(tree,"input tree");
    int a=1;
    int s=0;
    int nbChild=0;
    int n=0;
    int m=0;
    do{
        c = readChar(tree,"input tree");
        if (c==')'){
            a--;
            s=0;
            nbChild=0;
            while (!pileNode.empty() && s!=-1) {
                s=pileNode.top();pileNode.pop();
                if (s!=-1){
                    nbChild++;
                }
            }
            string sp=readSupport(tree,"input tree");
            if (a>0) {
                double d = readdouble(tree,"input tree");
            }
            pileNode.push(1);
            n++;
            m++;
        }
        else if (c!='(' && c!=',' && c!=-1 && c!=';' && c!='\n'){
            string s = readLabel(c,tree,a);
            pileNode.push(1);
            double d =readdouble(tree,"input tree");
            m++;
        }
        else if (c=='(') {a++;pileNode.push(-1);}
        else if (c=='\n') {
            c=readChar(tree,"input tree");
        }
    }  while (a>0);
    tree.close();
    m--;
    if (m==2*n) {
        opt->rooted=true;
    }
    else if (m==2*n+1 && nbChild==3){
        opt->rooted=false;
    }
    else if (nbChild==2) {
        opt->rooted=true;
    }
    else{
        opt->rooted=false;
    }
    if (opt->rooted){
        opt->nbINodes=n;
        opt->nbBranches=m;
    }
    else{
        opt->nbINodes=n+1;
        opt->nbBranches=m+1;
    }
}

int firstCharacter(string s,int p){
    while (p<s.length() && (s.at(p)<33 || s.at(p)>126)) {
        p++;
    }
    return p;
}

int lastCharacter(string s,int p){
    while (p<s.length() && (s.at(p)>=33 && s.at(p)<=126)) {
        p++;
    }
    return p;
}

int getInternalNodeId(Pr* pr,Node** nodes,string s){
    int k = getPosition(nodes,s,0,pr->nbBranches+1);
    if (k==-1 && (s.substr(0,4).compare("mrca")==0)){
        vector<int> mr;
        char c='(';
        int p = 5;
        while (c!=')'){
            int newp = s.find_first_of(",)",p);
            if (newp==-1){
                cerr<<s<<": wrong format"<<endl;
                exit(EXIT_FAILURE);
            }
            c=s.at(newp);
            string s1 = s.substr(p,newp-p);
            p = newp+1;
            int k1=getPosition(nodes,s1,0,pr->nbBranches+1);
            if (k1!=-1){
                mr.push_back(k1);
            }
            else{
                cerr<<"taxa "<<s1<<" not found"<<endl;
                exit(EXIT_FAILURE);
            }
        }
        if (mr.size()>0){ k=mrca(nodes,mr);}
    }
    return k;
}

int assignRecursive(int r,Node** nodes,int g){
    vector<int> children = nodes[r]->suc;
    int s = 0;
    for (int i=0; i<children.size(); i++) {
        if (nodes[children[i]]->rateGroup!=-1 && nodes[children[i]]->rateGroup!=g) {
            nodes[children[i]]->rateGroup = g;
            s+=1+assignRecursive(children[i],nodes,g);
        }
    }
    return s;
}

int assignRateGroupToSubTree(Subtree* subtree,Pr* pr,Node** nodes,int g){
    Pair* root = subtree->root;
    int r = getInternalNodeId(pr,nodes,root->name);
    bool toReroot = false;
    if (r >= pr->nbINodes) {
        subtree->tips.push_back(root);
        toReroot = true;
        root->include = false;
    }
    int s = 0;
    vector<int> tipsId;
    for (int i=0; i<subtree->tips.size(); i++) {
        Pair* tip = subtree->tips[i];
        int t = getInternalNodeId(pr,nodes,tip->name);
        nodes[t]->rateGroup = g;
        s++;
        tipsId.push_back(t);
    }
    if (toReroot) r = mrca(nodes,tipsId);
    if (root->include) {
        s += assignRecursive(r,nodes,g);
    }
    else{
        vector<int> children = nodes[r]->suc;
        for (int k=0; k<children.size(); k++) {
            if (nodes[children[k]]->rateGroup!=g) {
                s += assignRecursive(children[k],nodes,g);
            }
        }
    }
    return s;
}

void assignRateGroupToTree(Pr* pr,Node** nodes){
    vector<int> subroot;
    
    for (int i=0; i<pr->ratePartition.size(); i++) {
        Part* group = pr->ratePartition[i];
        for (int j=0; j<group->subtrees.size(); j++) {
            Pair* root = group->subtrees[j]->root;
            int r = getInternalNodeId(pr,nodes,root->name);
            if (contain(r,subroot)) {
                cout<<"Warning: "<<group->name<<" there are overlapped subtrees in the partition file"<<endl;
            }
            else{
                subroot.push_back(r);
            }
            if (root->include) {
                nodes[r]->rateGroup = -1;
            }
            else{
                vector<int> children = nodes[r]->suc;
                for (int k=0; k<children.size(); k++) {
                    nodes[children[k]]->rateGroup = -1;
                }
            }
        }
    }
    double nbBranchesPartition = 0;
    for (int i=0; i<pr->ratePartition.size(); i++) {
        Part* group = pr->ratePartition[i];
        for (int j=0; j<group->subtrees.size(); j++) {
            Subtree* subtree = group->subtrees[j];
            nbBranchesPartition += assignRateGroupToSubTree(subtree,pr,nodes,i+1);
        }
    }
    for (int i=0; i<pr->ratePartition.size(); i++) {
        Part* group = pr->ratePartition[i];
        for (int j=0; j<group->subtrees.size(); j++) {
            Pair* root = group->subtrees[j]->root;
            int r = getInternalNodeId(pr,nodes,root->name);
            if (root->include) {
                nodes[r]->rateGroup = i+1;
                nbBranchesPartition ++;
            }
            else{
                vector<int> children = nodes[r]->suc;
                for (int k=0; k<children.size(); k++) {
                    if (nodes[children[k]]->rateGroup != (i+1)) {
                        nodes[children[k]]->rateGroup = i+1;
                        nbBranchesPartition ++;
                    }
                }
            }
        }
    }
    
    if (nbBranchesPartition==pr->nbBranches) {
        pr->multiplierRate[0]=-1;//Full partition
    }
}

void calculate_tree_height(Pr* pr,Node** & nodes){
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
    vector<int> pre = preorder_polytomy(pr,nodes);
    for (vector<int>::iterator iter = pre.begin();iter!=pre.end();iter++){
        int i = *iter;
        if (i==0){
            nodes[i]->H = 0;
        } else{
            int P = nodes[i]->P;
            nodes[i]->H = nodes[P]->H + nodes[i]->B;
        }
    }
    double maxH =0;
    double maxHD = 0;
    for (int i = pr->nbINodes; i <= pr->nbBranches; i++){
        int P = nodes[i]->P;
        nodes[i]->H = nodes[P]->H + nodes[i]->B;
        if (nodes[i]->H > maxH){
            maxH = nodes[i]->H;
        }
        if (nodes[i]->D > maxHD){
            maxHD = nodes[i]->D;
        }
    }
    for (int i=0; i<=pr->nbBranches;i++){
        nodes[i]->H = maxH - nodes[i]->H;
        nodes[i]->HD = maxHD - nodes[i]->D;
    }
}

void splitExternalBranches(Pr* pr,Node** nodes){
    pr->ratePartition.clear();
    Part* part = new Part("externalBranches");
    for (int i=pr->nbINodes;i<=pr->nbBranches;i++){
        Pair* node = new Pair(false, nodes[i]->L);
        Subtree* subtree = new Subtree(node);
        part->subtrees.push_back(subtree);
    }
    pr->ratePartition.push_back(part);
    pr->multiplierRate.clear();
    pr->multiplierRate.push_back(1);
    pr->multiplierRate.push_back(1);
}

void splitLongBranches(Pr* pr,Node** nodes,double th){
    pr->ratePartition.clear();
    Part* part = new Part("longBranches");
    for (int i=pr->nbINodes;i<=pr->nbBranches;i++){
        if (nodes[i]->B > th){
            Pair* node = new Pair(false, nodes[i]->L);
            Subtree* subtree = new Subtree(node);
            part->subtrees.push_back(subtree);
        }
    }
    for (int i=1;i<pr->nbINodes;i++){
        if (nodes[i]->B > th){
            int p = nodes[i]->P;
            ostringstream np,ni;
            np<<p;
            if (nodes[p]->L == "") nodes[p]->L = "longBranch"+np.str();
            ni<<i;
            if (nodes[i]->L == "") nodes[i]->L = "longBranch"+ni.str();
            Pair* root = new Pair(false,nodes[p]->L);
            Pair* tip = new Pair(false,nodes[i]->L);
            vector<Pair*> tips;
            tips.push_back(tip);
            Subtree* subtree = new Subtree(root,tips);
            part->subtrees.push_back(subtree);
        }
    }
    pr->ratePartition.push_back(part);
    pr->multiplierRate.clear();
    pr->multiplierRate.push_back(1);
    pr->multiplierRate.push_back(1);
}

bool isIn(int i,vector<int> v){
    for (int j=0;j<v.size();j++){
        if (v[j]==i) return true;
    }
    return false;
}

vector<int> intersect(vector<int> v1, vector<int> v2){
    vector<int> intersect;
    for (int i=0;i<v2.size();i++){
        if (isIn(v2[i],v1)){
            intersect.push_back(v2[i]);
        }
    }
    return intersect;
}

double median(vector<double> array){
    sort(array.begin(), array.begin()+array.size());
    if (array.size() % 2 == 0){
        return (array[(array.size()/2)-1] + array[array.size()/2])/2;
    } else {
        return array[(array.size()-1)/2];
    }
}


/*void imposeMinBlen(ostream& file,Pr* pr, Node** nodes, double median_rate,bool medianRateOK){
    double minblen = pr->minblen;
    double round_time = pr->round_time;
    double m = 1./(pr->seqLength*median_rate);
    if (round_time <0){
        if (pr->inDateFormat == 2 || pr->inDateFormat == 1){
            round_time = 365;
        } else {
            if (m>=1) round_time = 100;
            else {
                round_time = 10;
                double mm = m;
                while (mm<1){
                    mm = mm*10;
                    round_time = round_time*10;
                }
            }
        }
    }
    string unit="";
    if (round_time==365) unit=" days";
    if (round_time==52) unit=" weeks";
    if (medianRateOK && pr->minblen<0){
        if (pr->inDateFile!="" || pr->inDateFormat==2 || pr->round_time!=-1){
            minblen = round(round_time*m)/(double)round_time;
            if (pr->minblenL < 0){
                cout<<"Minimum branch length of time scaled tree (settable via option -u and -U): "<<m<<",\n rounded to "<<minblen<<" ("<<round(round_time*m)<<unit<<"/"<<round_time<<") using factor "<<round_time<<" (settable via option -R)"<<endl;
                file<<"Minimum branch length of time scaled tree (settable via option -u and -U): "<<m<<",\n rounded to "<<minblen<<" ("<<round(round_time*m)<<"/"<<round_time<<") using factor "<<round_time<<" (settable via option -R)\n";
            } else {
                cout<<"Minimum internal branches lengths of time scaled tree (settable via option -u):\n "<<m<<", rounded to "<<minblen<<" ("<<round(round_time*m)<<"/"<<round_time<<") using factor "<<round_time<<" (settable via option -R)"<<endl;
                cout<<"Minimum external branches lengths of time scaled tree was set to "<<pr->minblenL<<"\n (settable via option -U)"<<endl;
                file<<"Minimum internal branches lengths of time scaled tree (settable via option -u):\n "<<m<<", rounded to "<<minblen<<" ("<<round(round_time*m)<<"/"<<round_time<<") using factor "<<round_time<<" (settable via option -R)\n";
                file<<"Minimum external branches lengths of time scaled tree was set to "<<pr->minblenL<<"\n (settable via option -U)"<<endl;
            }
        }
        else {
            minblen = m;
            if (pr->minblenL < 0){
                cout<<"Minimum branch length of time scaled tree (settable via option -u and -U): "<<m<<endl;
                file<<"Minimum branch length of time scaled tree (settable via option -u and -U): "<<m<<"\n";
            } else {
                cout<<"Minimum internal branches lengths of time scaled tree (settable via option -u): "<<m<<endl;
                cout<<"Minimum external branches lengths of time scaled tree was set to "<<pr->minblenL<<" (settable via option -U)"<<endl;
                file<<"Minimum internal branches lengths of time scaled tree (settable via option -u): "<<m<<"\n";
                file<<"Minimum external branches lengths of time scaled tree was set to "<<pr->minblenL<<" (settable via option -U)"<<endl;
            }
        }
    } else if (!medianRateOK && pr->minblen<0){
        cout<<"Can not estimate minimum branch lengths for time scaled tree: the temporal\n constraints provided are not enough, or conflict."<<endl;
        cout<<"Minimum branch length is then set to 0 (settable via option -u and -U)."<<endl;
        std::ostringstream oss;
        oss<<"- Can not estimate minimum branch lengths for time scaled tree: the temporal\n constraints provided are not enough, or conflict. Minimum branch length is then set to 0.\n";
        pr->warningMessage.push_back(oss.str());
        minblen = 0;
        
    } else{
        if (pr->minblen < 0) minblen = 0;
        if (pr->minblenL <0){
            cout<<"Minimum branch length of time scaled tree was set to "<<minblen<<" (settable\n via option -u and -U)"<<endl;
        } else {
            cout<<"Minimum internal branches lengths of time scaled tree was set to "<<minblen<<"\n (settable via option -u)"<<endl;
            cout<<"Minimum external branches lengths of time scaled tree was set to "<<pr->minblenL<<"\n (settable via option -U)"<<endl;
        }
    }
    nodes[0]->minblen = minblen;
    double minblenL = minblen;
    if (pr->minblenL >= 0) minblenL = pr->minblenL;
    for (int i=1;i<=pr->nbBranches;i++){
        if (i<pr->nbINodes) {
            nodes[i]->minblen = minblen;
        } else{
            nodes[i]->minblen = minblenL;
        }
    }
}*/

double median_branch_lengths(Pr* pr,Node** nodes){
    vector<double> bl;
    for (int i=1;i<=pr->nbBranches;i++){
        if (nodes[i]->B >= pr->nullblen){
            bl.push_back(nodes[i]->B);
        }
    }
    if (bl.size()==0){
        cerr<<"Not any branch length >= "<<pr->nullblen<<" (informative branch length threshold set via option -l)"<<endl;
        exit(EXIT_FAILURE);
    }
    return median(bl);
}

void collapse(int i,int j,Pr* pr,Node** nodes,Node** nodes_new,int &cc,int* &tab, double toCollapse, bool useSupport, double* support){
    for (vector<int>::iterator iter=nodes[j]->suc.begin(); iter!=nodes[j]->suc.end(); iter++) {
        int s= *iter;
        if (s<pr->nbINodes && (abs(nodes[s]->B) <= toCollapse  || (useSupport && support[s]<= pr->support))) {
            tab[s]=-1;
            collapse(i,s, pr, nodes, nodes_new, cc,tab,toCollapse,useSupport,support);
        }
        else{
            nodes_new[s]->P=i;
            nodes_new[s]->B = nodes[s]->B + nodes_new[j]->B;
            if (s<pr->nbINodes) {
                tab[s]=cc;
                cc++;
            }
        }
    }
}

int collapseTree(Pr* pr,Node** nodes,Node** nodes_new,int* &tab, double toCollapse, bool& useSupport){
    double* support = new double[pr->nbINodes];
    int root = (int)(!pr->rooted);
    if (useSupport){
        for (int i=root; i<pr->nbINodes;i++){
                try {
                    support[i] = std::stod(nodes[i]->L.c_str());
                } catch (const std::invalid_argument&) {
                    useSupport = false;
                }
        }
        if (!useSupport){
                std::ostringstream oss;
                oss<<"- Can not read support values, invalid arguments\n";
                pr->warningMessage.push_back(oss.str());
        }
    }
    for (int i=root;i<=pr->nbBranches;i++){
        nodes_new[i]= new Node();
        nodes_new[i]->P=nodes[i]->P;
        nodes_new[i]->type=nodes[i]->type;
        nodes_new[i]->lower=nodes[i]->lower;
        nodes_new[i]->upper=nodes[i]->upper;
        nodes_new[i]->D=nodes[i]->D;
        nodes_new[i]->B=nodes[i]->B;
        nodes_new[i]->L=nodes[i]->L;
    }
    if (pr->ratePartition.size()>0) {
        for (int i=0;i<=pr->nbBranches;i++){
            nodes_new[i]->rateGroup = nodes[i]->rateGroup;
        }
    }
    int cc = root+1;//number of internal nodes reduced
    if (!pr->rooted){
        tab[0]=0;
    }
    tab[root]=root;
    if (pr->removeOutgroup == false && pr->fnOutgroup!=""){
        for (vector<int>::iterator iter = nodes[root]->suc.begin(); iter != nodes[root]->suc.end();iter++){
            if (*iter < pr->nbINodes){
                tab[*iter]=cc;
                cc++;
            }
            nodes_new[*iter]->P = root;
        }
        for (int i=root+1;i<pr->nbINodes;i++){
            if ((abs(nodes_new[i]->B) > toCollapse && (!useSupport || support[i] > pr->support)) || nodes[i]->P==root){
                for (vector<int>::iterator iter=nodes[i]->suc.begin(); iter!=nodes[i]->suc.end(); iter++) {
                    int s=*iter;
                    if (s<pr->nbINodes && (abs(nodes[s]->B) <= toCollapse  || (useSupport && support[s]<= pr->support))) {
                        tab[s]=-1;
                        collapse(i, s,  pr, nodes, nodes_new,cc,tab, toCollapse, useSupport,  support);
                    }
                    else {
                        nodes_new[s]->P=i;
                        if (s<pr->nbINodes) {
                            tab[s]=cc;
                            cc++;
                        }
                    }
                }
            }
        }
    } else {
        for (int i=root;i<pr->nbINodes;i++){
            if ((abs(nodes_new[i]->B) > toCollapse && (!useSupport || support[i] > pr->support)) || i==root){
                for (vector<int>::iterator iter=nodes[i]->suc.begin(); iter!=nodes[i]->suc.end(); iter++) {
                    int s=*iter;
                    if (s<pr->nbINodes && (abs(nodes[s]->B) <= toCollapse  || (useSupport && support[s]<= pr->support))) {
                        tab[s]=-1;
                        collapse(i, s,  pr, nodes, nodes_new,cc,tab, toCollapse, useSupport,  support);
                    }
                    else {
                        nodes_new[s]->P=i;
                        if (s<pr->nbINodes) {
                            tab[s]=cc;
                            cc++;
                        }
                    }
                }
            }
        }
    }
    for (int i=pr->nbINodes;i<=pr->nbBranches;i++){
        tab[i]=cc+i-pr->nbINodes;
    }
    delete[] support;
    return cc;
}

void collapseTreeReOrder(Pr* pr,Node** nodes,Pr* prReduced,Node** nodesReduced,int* &tab){
    int root = (int)(!pr->rooted);
    nodesReduced[root]=new Node();
    nodesReduced[root]->P=-1;
    nodesReduced[root]->type=nodes[root]->type;
    nodesReduced[root]->lower=nodes[root]->lower;
    nodesReduced[root]->upper=nodes[root]->upper;
    nodesReduced[root]->D=nodes[root]->D;
    for (int i=(root+1); i<=pr->nbBranches; i++) {
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
    for (vector<Date*>::iterator iter=pr->internalConstraints.begin();iter!=pr->internalConstraints.end();iter++){
        Date* no = (*iter);
        if (no->mrca.size()==0){
            no->id = tab[no->id];
        }
        else{
            vector<int> new_mrca;
            for (vector<int>::iterator iter = no->mrca.begin(); iter != no->mrca.end(); iter++){
                new_mrca.push_back(tab[*iter]);
            }
            no->mrca = new_mrca;
        }
    }
    if (pr->ratePartition.size()>0) {
        for (int i=root;i<=pr->nbBranches;i++){
            if (tab[i]!=-1) nodesReduced[tab[i]]->rateGroup = nodes[i]->rateGroup;
        }
    }
}

void collapseUnInformativeBranches(Pr* &pr,Node** &nodes,bool verbose){
    Node** nodes_new = new Node*[pr->nbBranches+1];
    int* tab = new int[pr->nbBranches+1];
    bool useSupport = (pr->support>=0);
    int nbC = collapseTree(pr, nodes, nodes_new,tab, pr->nullblen,useSupport);//nbC is the number of internal nodes reduced
    if (verbose){
        if (!useSupport) {
            cout<<"Collapse "<<(pr->nbINodes - nbC)<<" (over "<<(pr->nbINodes - (!pr->rooted) -1 )<<") internal branches having branch length <= "<<pr->nullblen<<"\n (settable via option -l)"<<endl;
        } else {
            cout<<"Collapse "<<(pr->nbINodes - nbC)<<" (over "<<(pr->nbINodes - (!pr->rooted) -1 )<<") internal branches having branch length <= "<<pr->nullblen<<"\n (settable via option -l) or support value <= "<<pr->support<<"\n (settable via option -S)"<<endl;
        }
        if (  (double)(pr->nbINodes - nbC)/(pr->nbINodes - (!pr->rooted) -1) > 0.1){
            ostringstream oss;
            oss<<"- "<<(pr->nbINodes - nbC)*100/(double)(pr->nbINodes - (!pr->rooted) -1)<<"% internal branches were collapsed.\n";
            pr->warningMessage.push_back(oss.str());
        }
    }
    Node** nodesReduced = new Node*[nbC+pr->nbBranches-pr->nbINodes+1];
    Pr* prReduced = new Pr(nbC,nbC+pr->nbBranches-pr->nbINodes);
    prReduced->copy(pr);
    collapseTreeReOrder( pr, nodes_new, prReduced, nodesReduced,tab);
    for (int i=(!pr->rooted);i<pr->nbBranches+1;i++){
        delete nodes_new[i];
    }
    delete[] nodes_new;
    delete[] tab;
    computeSuc_polytomy(prReduced, nodesReduced);
    pr = prReduced;
    nodes = nodesReduced;
}

bool checkTopology(Pr* pr,Node** nodes1, Node** nodes2){
    for (int i=pr->nbINodes;i<=pr->nbBranches;i++){
        if (nodes1[i]->L != nodes2[i]->L) return false;
    }
    for (int i=1;i<=pr->nbBranches;i++){
        if (nodes1[i]->P != nodes2[i]->P) return false;
    }
    return true;
}

/*bool checkTopology(Pr* pr,Node** nodes1, Node** nodes2){
    int* tab = new tab[pr->nbBranches+1];
    for (int i=pr->nbINodes;i<=pr->nbBranches;i++){
        if (nodes1[i]->L != nodes2[i]->L) return false;
    }
    for (int i=1;i<=pr->nbBranches;i++){
        if (nodes1[i]->P != nodes2[i]->P) return false;
    }
    return true;
}*/

double* rtt(Pr* pr,Node** nodes){
    double* r2t = new double[pr->nbBranches+1];
    for (int i=0;i<=pr->nbBranches;i++){
        double r=0;
        int j=i;
        while (j!=0){
            r+=nodes[j]->B;
            j=nodes[j]->P;
        }
        r2t[i]=r;
    }
    return r2t;
}

void starting_pointLower(Pr* pr,Node** nodes,list<int> & active_set){
    for (int i =0;i<=pr->nbBranches;i++) {
        //if (nodes[i]->type!='p') {
        if (nodes[i]->type=='l' || nodes[i]->type=='b') {
            activeLower(nodes[i]);
            nodes[i]->D=nodes[i]->lower;
            active_set.push_back(-i);
        }
        /*   else if (nodes[i]->type=='u') {
         activeUpper(nodes[i]);
         nodes[i]->D=nodes[i]->upper;
         active_set.push_back(-i);
         }
         }*/
    }
}

void starting_pointUpper(Pr* pr,Node** nodes,list<int> & active_set){
    for (int i =0;i<=pr->nbBranches;i++) {
        //if (nodes[i]->type!='p') {
        if (nodes[i]->type=='u' || nodes[i]->type=='b') {
            activeUpper(nodes[i]);
            nodes[i]->D=nodes[i]->upper;
            active_set.push_back(-i);
        }
        /* else if (nodes[i]->type=='l') {
         activeLower(nodes[i]);
         nodes[i]->D=nodes[i]->lower;
         active_set.push_back(-i);
         }
         }*/
    }
}
