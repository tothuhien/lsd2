#include <unistd.h>
#include "utils.h"

Pr* getOptions( int, char** );
Pr* getCommandLine( int, char**);
Pr* getInterface( );
void     printHelp( void );
string    getInputFileName( string );
string    getOutgroupFileName( string );
void     printInterface(ostream& result, Pr*);
void     setOptionsWithLetter( Pr* , char);
double   getInputReal( string );
double getInputRealOrE( string msg );
double getInputDate( string msg, int& type );
double   getInputPositiveReal( string );
int      getInputInteger( string );
int      getPositiveInputInteger( string );
string    getInputString( string );
bool     isOptionActivate( Pr*, char );
FILE*    openOutputFile( char** );
FILE*    myFopen( char*, char* );
double getInputNonNegativeRealOrE( string msg );
