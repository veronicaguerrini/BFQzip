
 /*
 * Setting
 */

 
#ifndef _PARAMETERS_H_
#define _PARAMETERS_H_

#include <fstream>

#define DEBUG 0

#define OMP 1

#if B==0
    #define binningQS 0
#else
    #define binningQS 1
#endif

#define SIZEBUFFER 1024     //Size of the buffer for I/O partial ebwt/LCP/DA/SA

#define SIZE_ALPHA 256  


using namespace std;

typedef unsigned char uchar;
typedef unsigned short ushort;
typedef unsigned int uint;
typedef unsigned long ulong;

#define dataTypedimAlpha uchar  //size of the alphabet (in biologic case 6 ($,A,C,G,N,T))

/* USE: 0 - below 255 (unsigned char)
 *		1 - between 256 and 65.536 (unsigned short)
 *		2 - between 65.536 and 4.294.967.296 (unsigned int)
 *		3 - otherwise (unsigned long)
 */

// Type size for Sequences Length (in biologic case 100)
#define dataTypeLengthSequences 0		

// Type size for Number of Sequences in the input file
#define dataTypeNumSeq 2		

// Type size for Number of Character in the input file (length of the BWT)
#define dataTypeNumChar 3		


//Set the types
#if dataTypeLengthSequences == 0
	#define dataTypelenSeq uchar
#elif dataTypeLengthSequences == 1
	#define dataTypelenSeq ushort
#elif dataTypeLengthSequences == 2
	#define dataTypelenSeq uint
#elif dataTypeLengthSequences == 3
	#define dataTypelenSeq ulong	
#endif

#if dataTypeNumSeq == 0
	#define dataTypeNSeq uchar
#elif dataTypeNumSeq == 1
	#define dataTypeNSeq ushort
#elif dataTypeNumSeq == 2
	#define dataTypeNSeq uint
#elif dataTypeNumSeq == 3
	#define dataTypeNSeq ulong	
#endif


#if dataTypeNumChar == 0
	#define dataTypeNChar uchar
#elif dataTypeNumChar == 1
	#define dataTypeNChar ushort
#elif dataTypeNumChar == 2
	#define dataTypeNChar uint
#elif dataTypeNumChar == 3
	#define dataTypeNChar ulong	
#endif

////////////

//Print of the output (BWT/DA/SA/LCP)
//Store a txt file containing (BWT/DA/SA/LCP)
#define printFinalOutput 1

//Verbose
#define verboseEncode 0
#define verboseDecode 0

//if you want to delete the partial files, please set it to 1
#define deletePartialBWT 0
#define deletePartialLCP 1 
#define deletePartialGSA 1 
#define deletePartialQS 0
#define deleteCycFiles 0

//Compute the quality score permutation associated to BWT   permutation
#define USE_QS 1


//if you want to compute the LCP array, please set it to 1
#define BUILD_LCP 0

//The pair (da[i], sa[i]) is the gsa[i]
//if you want to compute the SA array, please set it to 1
#define BUILD_SA 0
//if you want to compute the DA array \in [0..nSeq], please set it to 1
#define BUILD_DA 0
//if you want to compute the bit DA array, please set it to 1
//it works in internal memory
#define BUILD_DA_bit 0



//if you want to Store the 'end positions' of the end-markers (one for each sequence), please set it to 1
#define STORE_ENDMARKER_POS 0

//if BUILD_BCR_ALTERNATE=0 then BCR computes the eBWT/SA/DA/LCP in straightforward order of the sequences (lexicographical order)
//if BUILD_BCR_ALTERNATE=1 then BCR computes the eBWT/SA/DA/LCP in alternating order of the sequences (alternating lexicographical order) See paper...
#define BUILD_BCR_ALTERNATE 0 //0 --> else we compute alternate order the BWT of the sequences

//if BCR_SET=1 then BCR computes the EBWT (the input is a set) (one can have strings of different length, so BCR uses the symbol TERMINATE_CHAR_LEN) 
//if BCR_SET=0 then BCR computes the BWT (the input is a single sequence)
#define BCR_SET 1				

//if BCR_INPUT_IN_MEMORY==1, BCR loads the input file in a string and compute the BWT of the string (it computes the BWT of the reverse string),  
//if BCR_INPUT_IN_MEMORY==0, BCR reads from files (cyc files)
#define BCR_INPUT_IN_MEMORY 0	

//if KEEP_eBWT_IN_EXT_MEMORY==1, BCR uses files for partials ebwts
//if KEEP_eBWT_IN_EXT_MEMORY==0, BCR uses strings for partials ebwts
//In both cases, SA, DA, LCP are stored in files.
#define KEEP_eBWT_IN_EXT_MEMORY  1

//Save lengths of each sequence in a file .len
#define STORE_LENGTH_IN_FILE  0

//if OUTPUT_FORMAT == 0, the output format of BCR is at most 5 files - built one after the other
//if OUTPUT_FORMAT == 1, the output format of BCR is as the output of EGSA (.gesa file). BUILD_LCP, BUILD_DA and BUILD_SA must be set to 1. Please, set the types as in eGSA
//if OUTPUT_FORMAT == 2, the output format of BCR is a unique file .egsa. BUILD_LCP must be set to 1 (we do not use a struct), BUILD_DA and BUILD_SA could be set to a either 0 or 1.  Order: ebwt, lcp, da, sa
//if OUTPUT_FORMAT == 3, the output format of BCR is at most 5 files at the same time
//if OUTPUT_FORMAT == 4, the output format of BCR is at most 3 files (ebwt, da), lcp, sa
//if OUTPUT_FORMAT == 5, the output format of BCR is at most 3 files ebwt, (lcp, da), sa
//if OUTPUT_FORMAT == 6, the output format of BCR is at most 3 files ebwt, lcp, (sa, da)
#define OUTPUT_FORMAT 0

//if OUTPUT_linear_SuffixArray == 1, BCR also computes the SA of the concatenated strings 
//if OUTPUT_linear_SuffixArray == 0, BCR does not compute the SA of the concatenated strings 
#define OUTPUT_linear_SuffixArray 0

//Computes the EGSA by starting to another EGSA data structure
//See output EGSA of Felipe Louze
//if BUILD_BCR_FROM_BCRpartials == 1, BCR takes in input the BCR partial files and adds the symbols of new sequences.
//if BUILD_BCR_FROM_BCRpartials == 0, BCR starts from scratch
#define BUILD_BCR_FROM_BCRpartials 0


//BCR reads the cyc files already computed by transpose.cpp in a previous execution.
//if BCR_FROMCYC=0 then BCR builds the cyc files before, otherwise BCR does not build the cyc files.
#define BCR_FROMCYC	0			

#endif
