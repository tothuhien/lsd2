#include <unistd.h>
#include "utils.h"

#ifndef VERSION
#define VERSION "v1.4.12"
#endif
#ifndef BOLD
#define BOLD      "\033[00;01m"
#endif
#ifndef FLAT
#define FLAT      "\033[00;00m"
#endif
#ifndef LINE
#define LINE      "\033[00;04m"
#endif
//max length of a file name
#ifndef MAX_FILE_NAME
#define MAX_FILE_NAME 256
#endif
#ifndef MAX_NAME_LENGTH
#define MAX_NAME_LENGTH 256
#endif


Pr* getOptions( int, char** );
Pr* getCommandLine( int, char**);
Pr* getInterface( );
void     printHelp( void );
string    getInputFileName( string );
string    getOutgroupFileName( string );
void     printInterface(ostream& result, Pr*);
void     setOptionsWithLetter( Pr* , char);
double   getInputReal( string );
double   getInputPositiveReal( string );
int      getInputInteger( string );
int      getPositiveInputInteger( string );
string    getInputString( string );
bool     isOptionActivate( Pr*, char );
FILE*    openOutputFile( char** );
FILE*    myFopen( char*, char* );

