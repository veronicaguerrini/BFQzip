#ifndef BCRdecode_hpp
#define BCRdecode_hpp

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdio.h>

#include "parameters.h" // Defines ulong and uchar.
#include "../external/rankbv/rankbv.hpp"
#include <string.h>

#include <deque>


#define BUFFERSIZE 1024

#define SIZE_ALPHA 256

#define DIMBLOCK  1024     // 1048576

using std::string;
using std::vector;

class BCRdecode
{
public:
    BCRdecode(int,bool);
    ~BCRdecode();

    string  fileOutBwt;
    string  fileOutCyc;
    string  ext;
    
    struct sortElement {
        dataTypeNSeq seqN;
        dataTypeNChar posN;
    };
    
    char TERMINATE_CHAR;
    char TERMINATE_CHAR_LEN;
 
    /*** For fq_compression_ext ***/
    dataTypedimAlpha sizeAlpha; //number of symbols
    vector<dataTypeNChar> freq;
    
    //alphabet map
    dataTypedimAlpha alpha[SIZE_ALPHA];
    vector<dataTypedimAlpha> alphaInverse;
    
    dataTypeNChar n; //BWT length
    
    dataTypeNSeq numSeq; //number of input string
    
    vector <string> BWT_MOD;//to store modified bases
    vector <rankbv_t*> rbv; //to track modified bases
    
    bool ignore_headers; //ignore headers
    
        /***/
    
    dataTypeNChar SIZEBUFFERcycFiles;
    
    dataTypelenSeq lengthRead;    //Length of each text
    dataTypeNChar lengthTot;   //Total length of all texts without $-symbols
    dataTypeNChar lengthTot_plus_eof;   //Total length of all texts with $-symbols

    dataTypeNSeq nText;   //number total of texts in filename1
    int numthreads;
    
    dataTypeNChar** tableOcc; //contains the number of occurrences of each symbol
    
    std::vector< std::deque<sortElement> > vectPair;  //Is is used both encoding, decoding, searching.
         
    std::vector< std::deque<sortElement> > vectPairTmp;
    std::vector< std::vector< std::deque<sortElement> > > vectPairParallel;
    std::vector< std::vector< std::deque<sortElement> > > vectPairParallelTmp;
    
    std::vector <bool> vectInsTexts;
    dataTypeNSeq textToInsert;
    vector <dataTypeNChar> vectSizeCurrentPile;

    vector< vector< vector<dataTypeNChar> > > vectorOcc;
    vector <dataTypeNChar> numBlocksInPartialBWT;

    void BCRInverse(string, string, dataTypelenSeq);
    
    void computeVectorUnbuildBCR();
    
    int concatenateFastaOrFastq( string );
  
    int decodeBCRmultipleReverse();
    #if OMP
        int RecoverNsymbolsReverseByVector(int t_id, uchar * newSymb);
    #else
        int RecoverNsymbolsReverseByVector(uchar * newSymb);
    #endif
    
    dataTypeNChar rankManySymbolsByVector(FILE & , dataTypeNChar *, dataTypeNChar, uchar *, uchar *, FILE &);
    
	void transposeCycFile(string , dataTypelenSeq, dataTypeNSeq );
	ulong read_block(uchar *, ulong , ulong , dataTypeNSeq, FILE *);
	ulong write_block(uchar *, ulong , ulong , dataTypelenSeq, FILE *);
	
    int convertFromCycFileToFastaOrFastq(string, string);
	
	
    int findBlockToRead(dataTypeNChar *, dataTypedimAlpha , dataTypeNChar *, dataTypeNChar *);
    
private:

    dataTypedimAlpha illumina_8_level_binning(dataTypedimAlpha newqs);
    
    FILE * openFilePartialIn(dataTypedimAlpha currentPile);
    dataTypeNChar readOnFilePartial(uchar *buffer, dataTypeNChar toRead, FILE * InFileBWT);
    dataTypeNChar writeOnFilePartial(uchar *buffer, dataTypeNChar numchar, FILE * OutFileBWT);
    int closeFilePartial(FILE * pFile);


    
};

#endif
