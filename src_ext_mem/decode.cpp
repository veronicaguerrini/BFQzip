#include "../external/rankbv/rankbv.hpp"
#include "decode.hpp"


#include <assert.h>
#include <stdio.h>
#include <sys/stat.h>
#include <vector>
#include <string.h>     // std::string, std::to_string
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>      // std::stringstream

#include <math.h>

#include <omp.h>

#include <sys/resource.h>

using namespace std;


BCRdecode::BCRdecode (int n_thr, bool titles){
    freq.resize(SIZE_ALPHA,0);
    numSeq = 0;
    numthreads=n_thr; //number of threads set by the user
    TERMINATE_CHAR = '#';     //end-marker symbol, it must be lexicographically smaller than all the alphabet letters
	   TERMINATE_CHAR_LEN = '$';  //it is stored in cyc files, it is ignored by the algorithm, so it must not belong to the alphabet
    ignore_headers = not titles;
    cout << "\tIllumina_8_level_binning = " << binningQS << endl;

    char c_aux[500];
    strcpy (c_aux,"mkdir -p tmp/ ");
    assert (system (c_aux) == 0);
}



//returns a quality score according to Illumina 8-Level Binning 
dataTypedimAlpha BCRdecode::illumina_8_level_binning(dataTypedimAlpha newqs){
    
    if (newqs >= 40) newqs = 40;
    else if (newqs >= 35) newqs = 37;
    else if (newqs >= 30) newqs = 33;
    else if (newqs >= 25) newqs = 27;
    else if (newqs >= 20) newqs = 22;
    else if (newqs >= 10) newqs = 15;
    else if (newqs >= 2) newqs = 6;
    
    return newqs+33;
    
}


void BCRdecode::computeVectorUnbuildBCR()
{
    numBlocksInPartialBWT.resize(sizeAlpha);
    
    //Set number of blocks for each BWT-partial
    for (dataTypedimAlpha x = 0 ; x <= sizeAlpha-1; x++) {
        numBlocksInPartialBWT[x] = (dataTypeNChar)ceil((long double)freq[alphaInverse[x]]/DIMBLOCK);
    #if DEBUG
        std::cerr << "computeVectorUnbuildBCR: " << "freq[" << (unsigned int)x << "]= " << freq[alphaInverse[x]] << " and numBlocksInPartialBWT[" << (unsigned int)x << "]= " << numBlocksInPartialBWT[x] << "\n";
    #endif
    }
    
    
    // Start by allocating an array for array of arrays
    
    vectorOcc.resize(sizeAlpha);    //For each BWT-partial
    
    #if OMP==0
    time_t start,end;
    double dif;
    #endif
    
    
    // alphaInverse[x] is the symbol to which correspond bwt_x
    
    //For each BWT-partial
    //Read BWT-partials in parallel
    
    vector < FILE * > InFileBWT;
    
    #if OMP
        double d_total = omp_get_wtime();
        InFileBWT.resize(numthreads);
    #else
        InFileBWT.resize(1);
    #endif
    
    dataTypedimAlpha x=0;
    int t = 0;
    
    #if OMP
    #pragma omp parallel for default(shared) private(x) firstprivate(t) num_threads(numthreads) schedule(dynamic, 1)
    #endif
    for ( x=0; x < sizeAlpha; x++)
    {
        string filenameIn = "tmp/bwt_" + to_string((int)(x))+".aux";
        #if OMP
        t = omp_get_thread_num(); //t is the thread ID
        #endif
        InFileBWT[t] = fopen(filenameIn.c_str(), "rb");
        if (InFileBWT[t]==NULL) {
            #if OMP
            #pragma omp critical
            #endif
            {
                std::cerr << "openFilePartialIn: file currentPile=" << (unsigned int)x << ": Error opening: " << filenameIn << std::endl;
                exit (EXIT_FAILURE);
            }
        }
        
        fseek(InFileBWT[t], 0, SEEK_SET);
        
        // Allocate an array for each block of BWT-partial
        vectorOcc[x].resize(sizeAlpha);
        
        // Allocate an array of integers for each element of vectorOcc[x]
        for (dataTypedimAlpha y = 0 ; y < sizeAlpha; y++)   {      //For each block
            vectorOcc[x][y].resize(numBlocksInPartialBWT[x],0);
        }
        
        uchar *bufBlock = new uchar[DIMBLOCK];
        
        dataTypeNChar numBlock = 0;
        dataTypeNChar num_read = 0;
        
        
        #if OMP
            double tr_start=omp_get_wtime();
        #else
            time (&start);
        #endif
        
        //Read DIMBLOCK symbols in BWT-partial
        while( ( (num_read =  readOnFilePartial(bufBlock, DIMBLOCK, InFileBWT[t]) ) && (num_read > 0) )  &&  (numBlock < numBlocksInPartialBWT[x]))   //Added check on numBlocks
        {
            for (dataTypeNChar i=0; i<num_read; i++) {
                vectorOcc[x][alpha[(unsigned int)(bufBlock[i])]][numBlock]++;

            }
            numBlock++;
                
        }//end-while
        
        if ( !feof(InFileBWT[t]) && (numBlock > numBlocksInPartialBWT[x])) {
            std::cerr << "computeVectorUnbuildBCR: Error - The file contains more blocks than allocates.\n" ;
            exit(1);
        }
            
        //For FIXED x
        //Compute the sum cumulative for each BWT-partial
        for (dataTypedimAlpha z = 0 ; z < sizeAlpha; z++)  {      //For each symbol z
            for(dataTypeNChar y = 1; y < numBlocksInPartialBWT[x] ; y++)   {      //For each block y>1 of partial-BWT x
                vectorOcc[x][z][y]=vectorOcc[x][z][y-1] + vectorOcc[x][z][y];   //Sum the previous one: ie Blcok y and block y-1
            }
        }
            
        //VectorOcc[x] stores info similar to row tableOcc[x] but keeping track of symbol occurrences into blocks
        //N.B. VectorOcc[x][z][y] stores the total number of z-occurrences in the BWT-partial corresponding to alpha[x] symbol up to the y-th block
            
            
        #if OMP
        #pragma omp critical
            {
                std::cerr << "TIME THREAD " << t << " = " << omp_get_wtime()- tr_start << "(in seconds) and numBlocksInPartialBWT[ " << (unsigned int)x << " ]= " << numBlocksInPartialBWT[x] << "\n\n";
            }
        #else
            time (&end);
            dif = difftime (end,start);
            
            std::cerr << "computeVectorUnbuildBCR: the cycle for " << "symbol = " << (int)x << " tooks " << dif << " seconds" << " and numBlocksInPartialBWT[ " << (unsigned int)x << " ]= " << numBlocksInPartialBWT[x] << "\n";
        #endif
        
            
        delete [] bufBlock;
            
        fclose(InFileBWT[t]);

            
    }//end-for
        
    #if OMP
        std::cout << "Time: " << omp_get_wtime()-d_total << std::endl;
    #endif
        
    #if DEBUG==1
        for (dataTypedimAlpha x = 0 ; x < sizeAlpha; x++) {
            std::cerr << "x = " << (unsigned int)x << " For the " << alphaInverse[x] << "-BWT-partial: the number of symbols is " << freq[alphaInverse[x]] << ".\n";
            std::cerr << "Number of blocks of the symbol " << numBlocksInPartialBWT[x] << "\n";
            for(dataTypedimAlpha z = 0; z < sizeAlpha; ++z) {
                std::cerr << "Symbol " << (unsigned int)z << ":\t";
                for(dataTypeNChar y = 0; y < numBlocksInPartialBWT[x]; ++y) {
                    std::cerr << (int)vectorOcc[x][z][y];
                }
                std::cerr << "\n";
            }
            std::cerr << "\n";
        }
    #endif
}
    

void BCRdecode::BCRInverse (string fileInputHeaders, string fileOutDecode, dataTypelenSeq maxLength)
{
    
    fileOutBwt = "bwt_";
    char c_aux[500];
    strcpy (c_aux,"mkdir -p ");
    strcat (c_aux, "tmp/cycFilesDecode/");
    assert (system (c_aux) == 0);
    fileOutCyc = "tmp/cycFilesDecode/cyc.\0";
    
    ext = ".aux";

    //set nText
    nText = numSeq;
    
    lengthRead=maxLength;

    #if OMP
        if((int)nText < numthreads){
            std::cerr << "The number of sequences is " << nText << ", which is too small to use " << numthreads << " threads.\n Number of threads used: " << nText << std::endl;
            numthreads = nText;
        }
    #endif
    
    time_t startI,endI;
    double difI;
    time (&startI);
    std::cerr << "Invert BWT by Backward direction." << std::endl;
    std::cerr << "\nStart decodeBCRmultipleReverse " << startI << " seconds\n";

    //run decodeBCRmultipleReverse
    textToInsert = nText;
    assert ( decodeBCRmultipleReverse() == 1);
    
    time (&endI);
    difI = difftime (endI,startI);
    std::cerr << "End   decodeBCRmultipleReverse " << endI << " seconds\n";
    std::cerr << "decodeBCRmultipleReverse tooks " << difI << " seconds\n";

    string newFilename;
    #if USE_QS == 1
        newFilename = fileOutDecode + ".fq";
    #else
        newFilename = fileOutDecode + ".fa";
    #endif


    time_t startT,endT;
    time (&startT);
    
    std::cerr << "\nStart convertFromCycFileToFastaOrFastq " << startT << " seconds\n";
    
    //run convertFromCycFileToFastaOrFastq
    assert (convertFromCycFileToFastaOrFastq(newFilename, fileInputHeaders) == 1);
    
    time (&endT);
    std::cerr << "End   convertFromCycFileToFastaOrFastq " << endT << " seconds\n";
    std::cerr << "convertFromCycFileToFastaOrFastq tooks " << difftime (endT,startT) << " seconds\n";
    
    #if OMP
    //concatenate multiple FastaOrFastq files: newFilename_1....newFilename_(numthread-1) in newFilename
    if(numthreads>1){
        time (&startT);
        std::cerr << "\nStart concatenateFastaOrFastq " << startT << " seconds\n";
        assert(concatenateFastaOrFastq(newFilename)==1);
        time (&endT);
        std::cerr << "End   concatenateFastaOrFastq " << endT << " seconds\n";
        std::cerr << "concatenateFastaOrFastq tooks " << difftime (endT,startT) << " seconds\n";
    }
    #endif
    
    //Free the memory
    for (dataTypedimAlpha j = 0 ; j < sizeAlpha; j++) {
        delete [] tableOcc[j];
        tableOcc[j] = NULL;
    }
    delete [] tableOcc;


}

int BCRdecode::concatenateFastaOrFastq (string Filename) {
    
    string toCopyFile;
    int ret;
    
    //Open newFilename to append
    std::ofstream firstFile;
    firstFile.open(Filename.c_str(), std::ios_base::app);
    
    for ( int i = 1; i<numthreads ; i++){
        //Open newFilename_i
        toCopyFile = Filename + "_" + to_string(i);
        std::ifstream nextFile(toCopyFile.c_str()); // Open for reading
        
        string line;
        
        while(getline(nextFile,line))//Append content to newFilename
            firstFile << line << "\n";
        
        nextFile.close();
        
        //Remove it
        ret = remove(toCopyFile.c_str());
        if(ret == 0)
            std::cout << "File " << toCopyFile << " deleted successfully\n";
        else
            std::cout << "Error: unable to delete file " << toCopyFile << std::endl;
    }//end-for
    firstFile.close();
    
    return 1;
}

        
int BCRdecode::findBlockToRead(dataTypeNChar *counters, dataTypedimAlpha currentPile, dataTypeNChar *toRead, dataTypeNChar *numBlock) {
        //Find the block numblock, where the position toRead is
    
        *numBlock = (dataTypeNChar)floor((long double)((*toRead-1)/DIMBLOCK)) ;  //The smallest integral value NOT less than x.
    
        assert(*numBlock < numBlocksInPartialBWT[currentPile]);
    
        if (*numBlock > 0) {
           for (dataTypedimAlpha r=0; r<sizeAlpha; r++)
              counters[r] =  vectorOcc[currentPile][r][(*numBlock)-1];   //vectorOcc is indexed by 0, so we have numBlock-1
           *toRead = *toRead - (*numBlock*DIMBLOCK);  //Number of symbols that we must read yet. it could be = DIMBLOCK
        }
        return 1;
    }
        
    
    //Multiple Decoding the sequences (Build reverse sequence)
    //Reconstruct m sequences backwards by threading through the mapping and reading the characters off of L.

    //fileOutBWT is the suffix of the filename of the partial BWTs
    //fileOutCyc is the prefix of the lengthRead-filename (transpose texts: cyc.i.txt)


ulong BCRdecode::read_block(uchar *block, ulong n, ulong m, dataTypeNSeq max_m, FILE *fileIn) {
    ulong  numEle=0, numEleTot=0;
    if (m == max_m) {
        numEle=fread(block,sizeof(uchar), n*m, fileIn);
        numEleTot += numEle;
    }
    else {
        for (ulong i = 0; i < n; i++) {
            numEle=fread(block+(i*m),sizeof(uchar),m, fileIn);
            numEleTot += numEle;
            fseek(fileIn, max_m-m, SEEK_CUR); //skip max_m-m 
        }
    }
    return numEleTot;
}

ulong  BCRdecode::write_block(uchar *block, ulong n, ulong m, dataTypelenSeq max_m, FILE *fileOut) {
    ulong  numEle=0, numEleTot=0;
    if (m == max_m) {
        numEle=fwrite(block,sizeof(uchar), n*m, fileOut);
        numEleTot += numEle;
    }
    else {
        for (ulong i = 0; i < n; i++) {
            numEle =fwrite(block+(i*m),sizeof(uchar),m, fileOut);           
            numEleTot += numEle;
            fseek(fileOut, max_m-m, SEEK_CUR); //skip max_m-m
        }
    }
    return numEleTot;
}


void BCRdecode::transposeCycFile(string filename, dataTypelenSeq row, dataTypeNSeq col)
{
	const ulong BLOCK_SIZE = 1000;	
	ulong row_block;
	ulong col_block;
	
	string out = filename + "_T";
    
    FILE *cyc_file = NULL;
	cyc_file = fopen(filename.c_str(), "r");
	if (cyc_file == NULL) {
		#if OMP
		#pragma omp critical
		#endif
		{
            std::cerr << "transposeCycFile: could not open file " << filename << " !" << std::endl;
			exit (EXIT_FAILURE);
        }
    }
	
    FILE *cyc_file_T = NULL;
    cyc_file_T = fopen(out.c_str(), "w");
	
	cerr << "transposeCycFile: writing " << out << "." << endl;

	if (BLOCK_SIZE<row)
	    col_block = (col<BLOCK_SIZE)?col:BLOCK_SIZE;
	else
		col_block = (col<BLOCK_SIZE*BLOCK_SIZE/row)?col:BLOCK_SIZE*BLOCK_SIZE/row;
	
	cerr << "col_block " << col_block << endl;
	
    row_block = (row<BLOCK_SIZE*BLOCK_SIZE/col_block)?row:BLOCK_SIZE*BLOCK_SIZE/col_block;

	cerr << "row_block " << row_block << endl;
	
	//Allocate block_in and block_out
    uchar *block_in = new uchar[row_block * col_block];
    uchar *block_out = new uchar[col_block * row_block];

    ulong  numEleRead=0, numEleWrite=0;

    for (ulong i = 0; i < row; i += row_block) 
	{
        ulong row_block_i = (row_block < (row - i))?row_block:row-i;
		
        for (ulong j = 0; j < col; j += col_block) 
		{
            ulong col_block_j = (col_block< (col - j))?col_block:col-j;

			//Read block
			fseek(cyc_file, i*col+j, SEEK_SET);
			
            numEleRead=read_block(block_in, row_block_i, col_block_j, col, cyc_file);

			for (ulong r = 0; r < row_block_i; r++) {
				for (ulong c = 0; c < col_block_j; c++) {
					block_out[c*row_block_i + r] = block_in[r*col_block_j + c];
				}
			}
			
			//Write block
            fseek(cyc_file_T, j*row+i, SEEK_SET);
			
            numEleWrite=write_block(block_out, col_block_j, row_block_i, row, cyc_file_T);

			#if DEBUG
				if (numEleRead!=numEleWrite) {
					#if OMP
					#pragma omp critical
					#endif
					{
						cerr << "*** i=" << (int)i << " j= " << (int)j << ", numEleRead!=numEleWrite: " << numEleRead << "!=" << numEleWrite << ", m= " << m << " n= " << n << endl;
						exit(1);
					}
				}
			#endif
			
        }
    }

	delete [] block_in;
	delete [] block_out;
	
    fclose(cyc_file);
    fclose(cyc_file_T);
	
}


int BCRdecode::decodeBCRmultipleReverse()
    {
        vectInsTexts.resize(nText,0);
        
        sortElement pair;
        
        #if OMP
            int j;
            vectPairParallel.resize(numthreads);
            vectPairParallelTmp.resize(numthreads);
            for (j=0; j<numthreads; j++){
                vectPairParallel[j].resize(sizeAlpha);
                vectPairParallelTmp[j].resize(sizeAlpha);
            }
            //vectPairParallel[j][pile] contains the pairs for the j-th thread
        #else
            vectPair.resize(sizeAlpha);
            vectPairTmp.resize(sizeAlpha);
        #endif

        //Position of $ in F
        dataTypedimAlpha pile=alpha[(unsigned int)(TERMINATE_CHAR)];
        
        #if OMP==0
            for (dataTypeNSeq g = 0 ; g < nText; g++) {
                pair.posN = g + 1;
                pair.seqN = g;
                vectPair[pile].push_back(pair);
            }
            #if DEBUG==1
                std::cerr << "The Initial triples of " << (int)pile << " in first column are!"<< std::endl;
                for (std::list<sortElement>::const_iterator it = vectPair[pile].begin(); it != vectPair[pile].end(); ++it)
                    std::cout <<  "\t" <<  it->posN << "\t" << it->seqN << std::endl;
                std::cerr << "vectSizeCurrentPile[" << (int)pile << "]= " << vectSizeCurrentPile[pile] << std::endl;
            #endif
        #endif

        //newSymb
        uchar *newSymb = NULL;

        #if USE_QS == 1
            newSymb = new uchar[2*nText];
        #else
            newSymb = new uchar[nText];
        #endif
        
        #if OMP==0
            static FILE *InfileOutDecodeCyc;
            string filename;
            time_t start,end;
        #endif
        
        #if OMP
            vector < FILE * > InfileOutDecodeCyc;
            InfileOutDecodeCyc.resize(numthreads);
        
            //numPairsIn is # pairs in vectPairParallel
            dataTypeNSeq numPairsIn=(dataTypeNChar)ceil((long double)nText/numthreads);
            #if DEBUG
                cout << "numPairsIn:" << numPairsIn << endl;
            #endif
        
            #pragma omp parallel default(shared) private(pair) firstprivate(numPairsIn) num_threads(numthreads)
            {
                int tid=omp_get_thread_num();//id_thread
                double start=omp_get_wtime();
                
                dataTypeNSeq s_ind = tid*numPairsIn; //start index
                
                #pragma omp master
                {
                    std::cout << "Number of threads: " << omp_get_num_threads() << std::endl;
                }
				
                if(tid==numthreads-1){//last thread
                    numPairsIn=nText-s_ind; //set numPairsIn for last thread
                }
						
                //push pairs vectPairParallel[tid]
                for (dataTypeNSeq ind=s_ind; ind < s_ind+numPairsIn; ind++) {
                    pair.posN = ind + 1;
                    pair.seqN = ind;
                    vectPairParallel[tid][pile].push_back(pair);
                }
				
                string filename;
				
				//We write each column of length nText (2*nText if USE_QS=1) in filename
				filename = fileOutCyc + to_string(tid) + ".txt";
				
				
				InfileOutDecodeCyc[tid] = fopen(filename.c_str(), "w");
				
                if (InfileOutDecodeCyc[tid]==NULL) {
                    #pragma omp critical
					{
                        std::cerr << "decodeBCRmultipleReverse: could not open file " << filename << " !" << std::endl;
						exit (EXIT_FAILURE);
                    }
                }
					
				
                dataTypelenSeq m;
				
				//As we recover the symbol in reverse order, I store the first found symbol in cyc.(length-1).txt file
                //and the last found symbol in cyc.0.txt file
                for (m = lengthRead ; m > 0 ; m--) {
                    
                    #if DEBUG
                    #pragma omp critical
                    {
                        cout << "***RecoverNsymbolsReverseByVector: t_id: " << tid << " m: " << (int)m << endl;
                    }
                    #endif
                    
                    assert ( RecoverNsymbolsReverseByVector(tid, newSymb) == 1);

                    #if USE_QS == 1
                        dataTypeNChar numcharWrite = fwrite (&newSymb[2*s_ind], sizeof(uchar), 2*numPairsIn , InfileOutDecodeCyc[tid]);
                        assert( numcharWrite == 2*numPairsIn); // we should always read the same number of characters
                    #else
                        dataTypeNChar numcharWrite = fwrite (&newSymb[s_ind], sizeof(uchar), numPairsIn , InfileOutDecodeCyc[tid]);
                        assert( numcharWrite == numPairsIn); // we should always read the same number of characters
                    #endif
                    
                }//end-for
                
				fclose(InfileOutDecodeCyc[tid]);
				
				
				cerr << "Transpose" << endl;
				
				//transpose
				#if USE_QS == 1
					transposeCycFile(filename, lengthRead, 2*numPairsIn);
				#else
					transposeCycFile(filename, lengthRead, numPairsIn);
				#endif
				
                #pragma omp critical
                {
                    std::cerr << "decodeBCRmultipleReverse: THREAD = " << tid << " tooks " << omp_get_wtime()-start << " seconds " << "\n\n";
                }
		
            }//end-pragma
            #else
				filename = fileOutCyc + "txt";
                InfileOutDecodeCyc = fopen(filename.c_str(), "w");
                if (InfileOutDecodeCyc==NULL) {
                    std::cerr << "decodeBCRmultipleReverse: could not open file " << filename << " !" << std::endl;
                    exit (EXIT_FAILURE);
                }
					
                for (dataTypelenSeq m = lengthRead ; m > 0 ; m--) {
                    
                    time (&start);
                    
                    assert ( RecoverNsymbolsReverseByVector(newSymb) == 1);

                    #if USE_QS == 1
                        dataTypeNChar numcharWrite = fwrite (newSymb, sizeof(uchar), 2*nText , InfileOutDecodeCyc);
                        assert( numcharWrite == 2*nText); // we should always read the same number of characters
                    #else
                        dataTypeNChar numcharWrite = fwrite (newSymb, sizeof(uchar), nText , InfileOutDecodeCyc);
                        assert( numcharWrite == nText); // we should always read the same number of characters
                    #endif
                   
                    time (&end);
                    #if DEBUG
                        std::cerr << "decodeBCRmultipleReverse: the cycle for position = " << (int)m << " tooks " << difftime (end,start) << " seconds" << "\n\n";
                    #endif
                    
                }//end-for
				
				fclose(InfileOutDecodeCyc);
				
				#if USE_QS == 1
					transposeCycFile(filename, lengthRead, 2*nText);
				#else
					transposeCycFile(filename, lengthRead, nText);
				#endif
				
            #endif
        
        delete [] newSymb;

        return 1;
    }
    
    
//It is used to reconstruct m sequences backwards by threading through the mapping and reading the characters off of L.
#if OMP
        int BCRdecode::RecoverNsymbolsReverseByVector(int t_id,  uchar * newSymb)    {
#else
        int BCRdecode::RecoverNsymbolsReverseByVector(uchar * newSymb)    {
#endif
    
    string filename;
    FILE *InFileBWT;
    #if USE_QS
        FILE *InFileQS=NULL;
        string filenameQS;
    #endif
    
    dataTypeNChar toRead = 0;
    dataTypeNChar *counters = new dataTypeNChar[sizeAlpha];  //it counts the number of each symbol into the i-Pile-BWT
    uchar *bufferBlock = new uchar[DIMBLOCK];
            
    dataTypedimAlpha currentPile=0, newCurrentPile=0;
    
    dataTypeNChar posInPile=0;
    
    sortElement pair, newPair;
            
    while(currentPile < sizeAlpha){
        
        //Set counters to 0
        for (dataTypedimAlpha i = 0 ; i < sizeAlpha; i++)
            counters[i]=0;
            
        filename = "tmp/" + fileOutBwt + to_string((int)(currentPile)) + ext;
        
        //Open BWT-partial
        InFileBWT = fopen(filename.c_str(), "rb");
        if (InFileBWT==NULL) {
            #if OMP
                #pragma omp critical
            #endif
            {
                std::cerr << "RecoverNsymbolsReverseByVector: " << filename << " file error opening " << std::endl;
                exit (EXIT_FAILURE);
            }
        }
            
        #if USE_QS
            filenameQS = "tmp/" + fileOutBwt + "qs_" + to_string((int)(currentPile)) + ext;
            InFileQS = fopen(filenameQS.c_str(), "rb");
            if (InFileQS==NULL) {
                #if OMP
                    #pragma omp critical
                #endif
                {
                    std::cerr << "RecoverNsymbolsReverseByVector: " << filenameQS << " file error opening " << std::endl;
                    exit (EXIT_FAILURE);
                }
            }
        #endif

        uchar foundSymbol, foundQual='\0';
        
        dataTypeNChar numBlock = 0;
           
        posInPile=0; //posInPile is the position in vectPair[currentPile]
        
        #if OMP
        
            #if DEBUG
                #pragma omp critical
                {
                    cout << "t_id: " << t_id << ", currentPile: " << (int)currentPile << ", vectPairParallel[t_id][currentPile].size: " << vectPairParallel[t_id][currentPile].size() << endl;
                }
            #endif
        
        while (!vectPairParallel[t_id][currentPile].empty() ) {
            //Pick up a pair
            pair = vectPairParallel[t_id][currentPile].front();
            //Remove that pair
            vectPairParallel[t_id][currentPile].pop_front();
        #else
        while (!vectPair[currentPile].empty() ) {
            //Pick up a pair
            pair = vectPair[currentPile].front();
            //Remove that pair
            vectPair[currentPile].pop_front();
        #endif

            if (vectInsTexts[pair.seqN]==0) { //Retrieve the BWT-symbol to insert in newSymb[i]
                    
                #if DEBUG==1 && OMP==0
                    std::cerr << "Inside While Sequence number k= " << posInPile << "\n";
                    std::cerr << "\t" << " P["<< posInPile <<"]=" << (dataTypeNChar)pair.posN << " N["<< posInPile <<"]=" << (dataTypeNSeq)pair.seqN << "\n";
                #endif
                    
                //For any character (of differents sequences) in the same pile
                foundSymbol = '\0';
                #if USE_QS
                    foundQual = '\0';
                #endif
                
                toRead = pair.posN; //skip pair.posN symbols

                //set to zero counters
                for (dataTypedimAlpha i = 0 ; i < sizeAlpha; i++)
                        counters[i]=0;
                    
                if (toRead > 0) {
                        //if ToRead > dimBlock, we can use vectorOcc in order to find the occorrences in the blocks preceding the block where the position toRead is.
                    assert (findBlockToRead(counters, currentPile, &toRead, &numBlock) == 1);
                }

                if (toRead > 0 && toRead <= DIMBLOCK) {   //If toRead == DIMBLOCK, because I can need to known foundSymbol character
                       
                    fseek (InFileBWT, numBlock*DIMBLOCK, 0);
                    #if USE_QS
                        fseek (InFileQS, numBlock * DIMBLOCK, 0 );
                    #endif
                        
                    #if DEBUG==1 && OMP==0
                        std::cerr << "RecoverNsymbolsReverseByVector Move file to the position " << numBlock*DIMBLOCK <<  "\n";
                    #endif
                    
                    //Read BWT-symbol
                    dataTypeNChar numchar;
                    
                    #if ( USE_QS == 1 )
                        dataTypeNChar toReadQS = toRead;
                    #endif
                    
                    numchar = fread(bufferBlock,sizeof(uchar),toRead,InFileBWT);
                    assert(numchar == toRead);  // we should always read/write the same number of characters
                    
                    #if DEBUG == 1 && OMP==0
                        std::cerr << "toRead " << toRead << " is less than DIMBLOCK " << DIMBLOCK << "\n";
                        std::cerr << "numchar (should be equal to toRead): " << numchar << "\n";
                    #endif
                    
                    foundSymbol = bufferBlock[numchar-1];
                    //The symbol is in the last position in the partial BWT that we have read.
                    
                    #if DEBUG == 1 && OMP==0
                        std::cerr << "** counters before:\t";
                        for (dataTypedimAlpha i = 0 ; i < sizeAlpha; i++)
                            std::cerr << " " << counters[i];
                        std::cerr << "\n";
                    #endif
                    
                    //For each symbol in the buffer, it updates the number of occurrences into counters
                    for (dataTypeNChar r=0; r<numchar; r++)
                        counters[alpha[(unsigned int)bufferBlock[r]]]++;//increment the number of letter symbol into counters
                    
                    #if DEBUG == 1 && OMP==0
                        std::cerr << "counters after:\t";
                        for (dataTypedimAlpha i = 0 ; i < sizeAlpha; i++)
                            std::cerr << " " << counters[i];
                        std::cerr << "\n**\n";
                    #endif
                    
                    #if ( USE_QS == 1 )
                        if (toReadQS > 1) {
                            fseek(InFileQS, (toReadQS - 1), SEEK_CUR );
                        }
                        assert (fread( &foundQual, sizeof( uchar ), 1, InFileQS ) == 1);
                        
                        size_t pos1 = ftell( InFileBWT );
                        size_t pos2 = ftell( InFileQS );
                        assert( pos1 == pos2 );
                    #endif

                    #if DEBUG==1 && OMP==0
                        std::cerr << "RecoverNsymbolsReverseByVector foundSymbol " << (unsigned int)foundSymbol << " " << (char)foundSymbol <<  "\n";
                        #if ( USE_QS == 1 )
                            std::cerr << "RecoverNsymbolsReverseByVector foundQual " << (unsigned int)foundQual <<  "\n";
                        #endif
                    #endif
                        
                }  //end-if (toRead <= DIMBLOCK)
                    
                #if DEBUG==1 && OMP==0
                    std::cerr << "counters after FirstVector:\t";
                    for (dataTypedimAlpha i = 0 ; i < sizeAlpha; i++)
                            std::cerr << " " << counters[i];
                    std::cerr << "\n";
                    std::cerr << "pair.seqN = " << pair.seqN << " Symbol = " << foundSymbol << "\n";
                #endif
  
                //write foundSymbol and foundQual
                
                #if USE_QS
                    newSymb[2*pair.seqN] = (rankbv_access(rbv[currentPile],pair.posN-1)==1)?(BWT_MOD[currentPile][rankbv_rank1(rbv[currentPile],pair.posN-1)-1]):(foundSymbol);
                
                    #if binningQS
                        newSymb[2*pair.seqN+1] = illumina_8_level_binning(foundQual-33);
                    #else
                        newSymb[2*pair.seqN+1] = foundQual;
                    #endif
                
                #else
                    newSymb[pair.seqN] = (rankbv_access(rbv, pair.posN-1)==1)?(BWT_MOD[currentPile][rankbv_rank1(rbv,pair.posN-1)-1]):(foundSymbol);
                #endif
                
                //vectInsTexts[pair.seqN]==1 iff the i-th text has been recovered
                if (foundSymbol == TERMINATE_CHAR)
                    vectInsTexts[pair.seqN]=1;
        
                    
                //Set newPair
                    
                //.posN
                //counters[alpha[(unsigned int)foundSymbol]]= # foundSymbol occ. in currentPile
                newPair.posN = counters[alpha[(unsigned int)foundSymbol]];

                for (dataTypedimAlpha g = 0 ; g < currentPile; g++) {
                    //Sum #foundSymbol occ. in each pile g= 0...(currentPile-1)
                    newPair.posN = newPair.posN + tableOcc[g][alpha[(unsigned int)foundSymbol]];
                }
                //.seqN
                newPair.seqN=pair.seqN;

                #if DEBUG==1 && OMP==0
                    std::cerr << "RecoverNsymbolsReverseByVector -Result: vectInsTexts[vectTriple[k].seqN]==0 j: P[q]=" << (dataTypeNChar)newPair.posN <<  " N[q]=" << (dataTypeNSeq)newPair.seqN << std::endl << std::endl;
                #endif
                    
                newCurrentPile=alpha[(unsigned int)foundSymbol];
                
                #if OMP
                    vectPairParallelTmp[t_id][newCurrentPile].push_back (newPair);
                #else
                    vectPairTmp[newCurrentPile].push_back (newPair);
                #endif
                    
            } //end-if (vectInsTexts[vectTriple[k].seqN]==0)
            else {
                //INSERT $ in newSymb[i]
                    
                foundSymbol = TERMINATE_CHAR_LEN;

                #if USE_QS
                    newSymb[2*pair.seqN] = foundSymbol;
                    newSymb[2*pair.seqN+1] = '{';
                #else
                    newSymb[pair.seqN] = foundSymbol;
                #endif
                
                #if DEBUG==1 && OMP==0 && USE_QS==0
                    std::cerr << "RecoverNsymbolsReverseByVector +Result: vectInsTexts[" << pair.seqN << "]!=0 k=" << posInPile << " - : foundSymbol=" << (unsigned int)foundSymbol << " newSymb[pair.seqN]=" << newSymb[pair.seqN] <<  std::endl << std::endl;
                #endif
            }
                
            posInPile++;
        }//end-while ( !vectPairParallel[t_id][currentPile].empty() )
            
        fclose(InFileBWT);
        #if USE_QS
            fclose(InFileQS);
        #endif
                    
        currentPile++;

    }//end-while ( currentPile < sizeAlpha )
    
    delete [] counters;
    delete [] bufferBlock;
            
    #if DEBUG==1 && OMP==0 && USE_QS==0
        std::cerr << "NewSymbols " ;
        for (dataTypeNSeq g = 0 ; g < nText; g++) {
          std::cerr << (char)newSymb[g] << " ";
        }
    #endif
        
    #if OMP
        vectPairParallelTmp[t_id].swap(vectPairParallel[t_id]);
    #else
        vectPairTmp.swap(vectPair);
    #endif
        
    return 1;
}

int BCRdecode::convertFromCycFileToFastaOrFastq( string fileOutDecode, string fileHeaders)
{
    int t_id = 0;
    
    //Open cyc files
    vector <FILE *> inFilesCyc;

    #if OMP
        inFilesCyc.resize(numthreads);
        //for(; t_id<numthreads; t_id++)
    #else
        inFilesCyc.resize(1);
        string fileOutName = fileOutDecode;
    #endif
        //{ inFilesCyc[t_id].resize(lengthRead); }
    
    #if OMP
    dataTypeNChar ind = (dataTypeNChar)ceil((long double)nText/numthreads);
	
    #pragma omp parallel default(shared) firstprivate(t_id) num_threads(numthreads)
    {
        t_id = omp_get_thread_num();//id_thread
        dataTypeNChar s_ind = t_id * ind; //start index
        dataTypeNChar f_ind = (t_id+1)*ind;  //final index
        
        if(t_id==numthreads-1){//last iteration
            f_ind = nText; //final index = nText for last thread
        }
        double start=omp_get_wtime();
        
        string fileOutName = fileOutDecode;
        
        if(t_id!=0)
            fileOutName = fileOutName + "_" + to_string(t_id);
	#endif
        
        //Open cyc files
        #if OMP
            string filenameIn = fileOutCyc + to_string(t_id) + ".txt_T";
        #else
            string filenameIn = fileOutCyc + "txt_T";
        #endif
		
        inFilesCyc[t_id] = fopen(filenameIn.c_str(), "r");
		
        if (inFilesCyc[t_id]==NULL) {
            #if OMP
            #pragma omp critical
            #endif
            {   std::cerr << "convertFromCycFileToFastaOrFastq: error opening " << filenameIn << std::endl;
                exit (EXIT_FAILURE);
            }
		}
            
		fseek(inFilesCyc[t_id], 0, SEEK_SET);
            
        //Open a outFile 
        FILE *outFile = fopen( fileOutName.c_str(), "w");
        
        if (outFile == NULL){
            #if OMP
            #pragma omp critical
            #endif
            {
                std::cerr << "Error opening \"" << fileOutDecode << "\" file" << std::endl;
                exit ( EXIT_FAILURE );
            }
        }

        FILE *fileh;
		
		#if USE_QS
            char header[]="@\n";
        #else
            char header[]=">\n";
        #endif
		
		char plus[]="+\n";
		size_t len = 0;
		char *buf = NULL;
		 
		if(not ignore_headers){
			fileh = fopen(fileHeaders.c_str(), "r");
			if(fileh == NULL){
				#if OMP
				#pragma omp critical
				#endif
				{
					std::cerr << "Error opening \"" << fileh << "\" file" << std::endl;
					exit ( EXIT_FAILURE );
				}
			}
			for(dataTypeNChar i=0; i<s_ind;i++)
				getline(&buf, &len, fileh); // discard those @'s lines in [0, s_ind)
		}
		
		uchar *bufferSeq = new uchar[lengthRead+1];
		
		dataTypeNChar num;
		dataTypelenSeq index_rev, index_buf, index_s;
		uchar symbol;
		
		#if USE_QS == 1
			uchar *bufferQS = new uchar[lengthRead+1];
		#endif
        
        #if OMP
        for ( dataTypeNSeq j = s_ind; j < f_ind; j++ )
        #else
        for ( dataTypeNSeq j = 0; j < nText ; j++ )
        #endif
        {
			//Write headers
            if(not ignore_headers){
				ssize_t size = getline(&buf, &len, fileh); // @'s line
				fwrite(buf, sizeof(char), size, outFile);
			}
			else{
				fwrite(header, sizeof(char), 2, outFile);
			}
			
			num = fread (bufferSeq, sizeof( uchar ),lengthRead, inFilesCyc[t_id] );
			assert(num==lengthRead);

			#if USE_QS == 1
				num = fread (bufferQS, sizeof( uchar ),lengthRead, inFilesCyc[t_id] );
				assert(num==lengthRead);
			#endif
			
			index_rev = lengthRead-1;
			while((bufferSeq[index_rev]==TERMINATE_CHAR_LEN) || (bufferSeq[index_rev]==TERMINATE_CHAR) )
				index_rev--;
				
			index_buf = index_rev +1;
			index_s = index_buf/2;
			
			//Reverse strings
			for (dataTypelenSeq index = 0; index < index_s; index++)
			{
				symbol = bufferSeq[index];
				bufferSeq[index] = bufferSeq[index_rev];
				bufferSeq[index_rev] = symbol;
				#if USE_QS == 1
                    symbol = bufferQS[index];
					bufferQS[index] = bufferQS[index_rev];
					bufferQS[index_rev] = symbol;
                #endif
				index_rev--;
			}
			bufferSeq[index_buf] = '\n';
			
            //Write sequence
            fwrite(bufferSeq, sizeof(uchar),index_buf+1, outFile);

            #if USE_QS == 1
				bufferQS[index_buf] = '\n';
                //Write qualities
                fwrite(plus, sizeof(char), 2, outFile);
				fwrite(bufferQS, sizeof(uchar),index_buf+1, outFile);
            #endif
            
        }  //end-for

        fclose(outFile);
		if(not ignore_headers)
			fclose(fileh);
        free(buf);
		
		delete [] bufferSeq;
		#if USE_QS == 1
            delete [] bufferQS;
        #endif
		
		fclose( inFilesCyc[t_id] );
		
	#if OMP    
        #pragma omp critical
        {
            std::cerr << "convertFromCycFileToFastaOrFastq: TIME THREAD " << t_id << " = " << omp_get_wtime()- start << "\n\n";
        }  
    }//end-pragma
    #endif
    
    return 1;
}

dataTypeNChar BCRdecode::readOnFilePartial(uchar *buffer, dataTypeNChar toRead, FILE * InFileBWT) {
   dataTypeNChar numchar;

   numchar = fread(buffer,sizeof(uchar),toRead,InFileBWT);

   return numchar;
}

dataTypeNChar BCRdecode::writeOnFilePartial(uchar *buffer, dataTypeNChar numchar, FILE * OutFileBWT) {
   dataTypeNChar numcharWrite;
   numcharWrite = fwrite (buffer, sizeof(uchar), numchar , OutFileBWT);
   return numcharWrite;
}

BCRdecode::~BCRdecode()
{

}
