/*
 Interface for LSD as a library
 */

#ifndef LSD_H_
#define LSD_H_

/** define a new namespace for LSD */
namespace lsd {
    
    /**
     abstract input/output stream for interface with LSD
     added by BQM 2020-04-09
     */
    class InputOutputStream {
    public:
        /** input tree stream */
        istream *inTree;
        
        /** input outgroup stream */
        istream *inOutgroup;
        
        /** input date stream */
        istream *inDate;
        
        /** input partition stream */
        istream *inPartition;
        
        /** input bootstrap tree stream */
        istream *inBootstrapTree;
        
        /** input rate stream */
        istream *inRate;
        
        /** output result stream */
        ostream *outResult;
        
        /** output tree 1 stream */
        ostream *outTree1;
        
        /** output tree 2 stream */
        ostream *outTree2;
        
        /** output tree 3 stream */
        ostream *outTree3;
        
        /** constructor */
        InputOutputStream();
        
        /**
         constructor that initialise object as (i/o)stringstream
         @param tree tree string
         @param outgroup outgroup string
         @param date date string
         */
        InputOutputStream(string tree, string outgroup, string date,string rate,string partition, string bootstrap);
        
        /** destructor */
        virtual ~InputOutputStream();
        
        /**
         set the content of the tree stream
         @param str a tree string
         */
        virtual void setTree(string str);
        
        /**
         set the content of the outgroup stream
         @param str a string in the following format: first line is the number
         of outgroup, each of the following line is a name of the outgroup
         */
        virtual void setOutgroup(string str);
        
        /**
         set the content of the date stream
         @param str a string
         */
        virtual void setDate(string str);
        
        /**
         set the content of the partition stream
         @param str a string
         */
        virtual void setPartition(string str);
        
        /**
         set the content of the bootstrap tree stream
         @param str a string
         */
        virtual void setBootstrapTree(string str);
        
        /**
         set the content of the rate stream
         @param str a string
         */
        virtual void setRate(string str);
    };
    
    /**
     main function to call LSD
     @param argc input for main()
     @param argv input for main()
     @param inputOutput stream for input/output, nullptr for reading/writing files
     */
    int buildTimeTree( int argc, char** argv, InputOutputStream *inputOutput = nullptr);
    
} // namespace

#endif
