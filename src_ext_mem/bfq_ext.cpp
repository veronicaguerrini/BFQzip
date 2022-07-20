#include <iterator>
#include <map>
#include <iostream>
#include <fstream>
#include <assert.h>
#include <vector>
#include <algorithm>
#include <unistd.h>
#include <math.h>
#include <iomanip>
#include <stack>
#include <sstream>
#include <unordered_map>
#include <limits.h>

#include <omp.h>

#include "../external/rankbv/rankbv.hpp"
#include "../external/malloc_count/malloc_count.h"
#include "decode.hpp"
#include "parameters.h"

#define BUFFERSIZE 1024

#ifndef DEBUG
    #define DEBUG 1
#endif


string input_dna;
string input_qual;
string input_lcp;
string input_titles;
string output;

bool isFileHeaders = false;

/*
 * Debug mode variables: print the BWT and read names/qualities for each base
 */

bool debug = false; //print debug info
bool verbose = false; //verbose output
int max_id_len=20; //in read_info, store at most this number of chars for the IDs
vector<string> read_info;//if debug, store read coordinate for every BWT position
//------------

vector<string> read_ids;//ID of each read

dataTypeNSeq modified = 0;//count how many bases have been modified
dataTypeNSeq qs_smoothed = 0;//count how many qs have been modified
dataTypeNChar bases_inside;//total number of bases inside analyzed clusters

dataTypeNChar clusters_size = 0;//total number of bases inside clusters
dataTypeNChar num_clust_discarded = 0;
dataTypeNChar num_clust_amb_discarded = 0;
dataTypeNChar num_clust_mod = 0;

vector<dataTypeNChar> freqs; //temporary vector used to count frequency of bases inside the cluster
vector<dataTypeNChar> lowQS; //bitvector for qs symbols associated with dna symbols in any cluster (if 0 all occ. have low qs)

//vector < vector<dataTypeNChar> > vectorOcc;

BCRdecode *BCRdec;

vector<dataTypeNChar> fVect; //vector used to split BWT

//minimum LCP required inside clusters
int K_def = 16;
int K = -1;

//do not consider clusters smaller than this threshold
int m_def = 2;
int m = 0;

//default qs
char default_value_def = '>';   //QS = 29
char default_value = '\0';

//quality threshold
int quality_threshold_def = 20; 
int quality_threshold = -1; 

//threshold for frequent bases in clusters
float freq_threshold_def = 40.0;
float freq_threshold = 0.0;

//maximum read length
dataTypelenSeq maxLengthRead = 0;

//threads number
int num_threads_def = 1;
int num_threads = 0;

//terminator character at the end of the reads
char TERM = '#';

vector<bool> LCP_minima;//bitvector that stores the LCP minima
vector<bool> LCP_threshold;//bitvector that stores LCP values that exceed the threshold: >= K

             //A B C D E F G H I J K L M N O P Q R S T U V W X Y
int ORD[25] = {0,0,1,0,0,0,2,0,0,0,0,0,0,4,0,0,0,0,0,3,0,0,0,0,0};
#define ord(c) (ORD[c-65])

void help(){
    
    cout << "bfq_ext [options]" << endl <<
    "Options:" << endl <<
    "-h          Print this help." << endl <<
    "-e <arg>    Input eBWT file (A,C,G,T,#) of DNA (REQUIRED)." << endl <<
    "-q <arg>    Qualities permuted according to the DNA's ebwt (REQUIRED)." << endl <<
    "-a <arg>    LCP array file (REQUIRED)." << endl <<
    "-o <arg>    Output fastq (REQUIRED)." << endl <<
    "-l <arg>    Maximum read length (REQUIRED)." << endl <<
    "-k <arg>    Minimum LCP required in clusters. Default: " << K_def << "." << endl <<
    "-m <arg>    Minimum length of cluster to be processed. Default: " << m_def << "." << endl <<
    "-v <arg>    Quality score for constant replacement (if M=2). Default: " << (int)default_value_def-33 << "." << endl <<
    "-f <arg>    Percentage threshold for frequent bases in clusters. Default: " << freq_threshold_def << "." << endl <<
    "-t <arg>    Quality score threshold for trusted bases. Default: " << quality_threshold_def << "." << endl <<
    "-s <arg>    ASCII value of end-marker symbol. Default: " << int('#') << " (#)." << endl <<
    "-T <arg>    Number of threads. Default: " << num_threads_def << "." << endl <<
    "-D          Print debug info for each BWT position." << endl <<
    "-H <arg>    Original headers from file." << endl << endl <<
    
    "\nTo run bfq_ext, you must first build the extended Burrows-Wheeler Transform " <<
    "of the input DNA sequences and the corresponding permutation of base quality scores, and the associated LCP array." << endl;
    
    exit(0);
}


void buildFreq() {
    
    //Open BWT file
    FILE *InBWT = fopen(input_dna.c_str(), "rb");
    if (InBWT==NULL) {
        std::cerr << "Error opening " << input_dna << std::endl;
        exit (EXIT_FAILURE);
    }
    fseek(InBWT, 0, SEEK_END);
    BCRdec->n = ftell(InBWT)/sizeof(uchar);
    
    //Read BWT
    dataTypedimAlpha *bufferBWT = new uchar[BUFFERSIZE];
    dataTypeNChar numcharBWT;
    
    fseek(InBWT, 0, SEEK_SET);
    numcharBWT=fread(bufferBWT,sizeof(uchar), BUFFERSIZE,InBWT);
    
    while(numcharBWT>0){
        for (dataTypeNChar i=0; i < numcharBWT; i++){
            BCRdec->freq[(unsigned int)(bufferBWT[i])]++;
            if(bufferBWT[i] == TERM)
                BCRdec->numSeq++;
        }
        numcharBWT=fread(bufferBWT,sizeof(uchar), BUFFERSIZE,InBWT);
    }
    
    fVect.push_back(0);
    BCRdec->sizeAlpha = 0;
    
    //set alpha, alphaInverse freq and fVect
    for (uchar i = 0; i < UCHAR_MAX; ++i) {
        if (BCRdec->freq[i] > 0) {
            BCRdec->alpha[i] = BCRdec->sizeAlpha;
            BCRdec->alphaInverse.push_back(i);
            dataTypeNChar sum_pos = fVect[BCRdec->sizeAlpha]+BCRdec->freq[i];
            fVect.push_back(sum_pos);
            BCRdec->sizeAlpha++;
        }
    }
    if (BCRdec->freq[UCHAR_MAX] > 0) {
        BCRdec->alpha[UCHAR_MAX] = BCRdec->sizeAlpha;
        BCRdec->alphaInverse.push_back(UCHAR_MAX);
        dataTypeNChar sum_pos = fVect[BCRdec->sizeAlpha]+BCRdec->freq[UCHAR_MAX];
        fVect.push_back(sum_pos);
        BCRdec->sizeAlpha++;
    }
    #if DEBUG
    cout << "\nfVect: [";
    for(dataTypedimAlpha i = 0; i < fVect.size(); i++)
        cout << fVect[i] << ",";
    cout << "]\n";
    #endif
    fclose(InBWT);
    delete [] bufferBWT;
}


void splitIntoPartial() {
    
    //Open BWT, QS files
    vector < FILE * > InBWT, InQS;
    
    int t=0;
    
#if OMP
    InBWT.resize(num_threads);
    InQS.resize(num_threads);
    for(t=0;t<num_threads; t++)
#else
    InBWT.resize(1);
    InQS.resize(1);
#endif
    {
        InBWT[t] = fopen(input_dna.c_str(), "rb");
        InQS[t] = fopen(input_qual.c_str(), "rb");
        if (InBWT[t]==NULL) {
            std::cerr << "Error opening " << input_dna << "!" << std::endl;
            exit (EXIT_FAILURE);
        }
        if (InQS[t]==NULL) {
            std::cerr << "Error opening " << input_qual << "." << std::endl;
            exit (EXIT_FAILURE);}
    }
    
    //allocate tableOcc
    BCRdec->tableOcc = new dataTypeNChar*[BCRdec->sizeAlpha];
    //Counting for each pile, es. $-pile, A-pile, C-pile, G-pile, N-pile, T-pile
    for (dataTypedimAlpha j = 0 ; j < BCRdec->sizeAlpha; j++) {
            BCRdec->tableOcc[j] = new dataTypeNChar[BCRdec->sizeAlpha];
    }
    for (dataTypedimAlpha j = 0 ; j < BCRdec->sizeAlpha; j++){
            for (dataTypedimAlpha h = 0 ; h < BCRdec->sizeAlpha; h++)
                BCRdec->tableOcc[j][h]=0;
    }
    
    vector < FILE *> OutFileBWT;
    vector < FILE *> OutFileQS;
#if OMP
    OutFileBWT.resize(num_threads);
    OutFileQS.resize(num_threads);
#else
    OutFileBWT.resize(1);
    OutFileQS.resize(1);
#endif
    
    dataTypeNChar numEle=0;
    
    dataTypedimAlpha currentPile;
#if OMP
#pragma omp parallel for default(shared) private(currentPile) firstprivate(t) num_threads(num_threads) schedule(dynamic, 1) reduction(+:numEle)
#endif
    for (currentPile = 0 ; currentPile < BCRdec->sizeAlpha; ++currentPile) {
        
#if OMP
        t = omp_get_thread_num();//id_thread
        double start = omp_get_wtime();
#endif
        
        //Open BCR partial files
        dataTypeNChar numcharBWT, numcharQS, numWrite;
        dataTypedimAlpha *buffer = new dataTypedimAlpha[BUFFERSIZE];
        
        string fnOutBWT = "tmp/bwt_" + to_string((int)currentPile)+".aux";
        OutFileBWT[t] = fopen(fnOutBWT.c_str(), "wb");
        string fnOutQS = "tmp/bwt_qs_" + to_string((int)currentPile)+".aux";
        OutFileQS[t] = fopen(fnOutQS.c_str(), "wb");
        
        if (OutFileBWT[t]==NULL)
#if OMP
#pragma omp critical
#endif
        {
            std::cerr << "Error opening: " << fnOutBWT << std::endl;
            exit (EXIT_FAILURE);
        }
        fseek(InBWT[t], fVect[currentPile]*sizeof(dataTypedimAlpha), SEEK_SET);
        
        if (OutFileQS[t] ==NULL)
#if OMP
#pragma omp critical
#endif
        {
            std::cerr << "Error opening: " << fnOutQS << std::endl;
            exit (EXIT_FAILURE);
        }
        fseek(InQS[t], fVect[currentPile]*sizeof(uchar), SEEK_SET);
        
        //Build BCR partial files
        dataTypeNChar j = fVect[currentPile];
        
        dataTypeNChar toRead = BUFFERSIZE;  //read BUFFERSIZE symbols per time
        
        while (j < fVect[currentPile+1]) {
            
            if( (fVect[currentPile+1] - j) < toRead) //check remaining symbols to read
                toRead = fVect[currentPile+1] - j;
            
            //read toRead symbols
            numcharBWT = fread(buffer, sizeof(uchar),toRead, InBWT[t]);
            //write
            numWrite = fwrite(buffer, sizeof(uchar), toRead, OutFileBWT[t]);
            assert(numcharBWT == numWrite);
            
            //counting the number of occurrences in BWT of the currentPile
            for(dataTypeNChar i = 0; i < numcharBWT; i++ ){
                BCRdec->tableOcc[(unsigned int)currentPile][BCRdec->alpha[(unsigned int)buffer[i]]]++;
                numEle++;
            
            }
            
            //read toRead symbols
            numcharQS = fread(buffer, sizeof(uchar),toRead, InQS[t]);
            assert(numcharQS==numcharBWT);
            //write
            numWrite = fwrite (buffer, sizeof(uchar),toRead, OutFileQS[t]);
            assert(numcharQS == numWrite);
            
            j+=toRead;
            
        }  //end-while
        
        //Close partial files
        fclose(OutFileBWT[t]);
        fclose(OutFileQS[t]);
	    
	delete [] buffer;
        
#if OMP
#pragma omp critical
        {
            uchar c = (BCRdec->alphaInverse[currentPile]>32)?BCRdec->alphaInverse[currentPile]:'#';
            std::cerr << "splitIntoPartial: Pile = " << c << "\tTHREAD = " << t << " tooks " << omp_get_wtime()-start << " seconds " << "\n\n";
        }
#endif
        
    }  //end-for
    
    //Close BWT, QS files
#if OMP
    for(t=0;t<num_threads; t++)
#endif
    {
        fclose(InBWT[t]);
        fclose(InQS[t]);
    }
    
    assert(numEle==(BCRdec->n));
    
    
    //Build partial bitvectors and BWT_MOD
    BCRdec->rbv.resize(BCRdec->sizeAlpha,NULL);
    BCRdec->BWT_MOD.resize(BCRdec->sizeAlpha);
    
    for(dataTypedimAlpha i = 0; i < BCRdec->sizeAlpha; i++)
        BCRdec->rbv[i] = rankbv_create(BCRdec->freq[BCRdec->alphaInverse[i]],2);
}

void buildBitVectors(dataTypeNChar len){ //build LCP_threshold and LCP_minima
    
    // Open LCP array file
    FILE *InLCP = fopen(input_lcp.c_str(), "rb");
    if (InLCP==NULL){
        cerr << "Error opening " << input_lcp << ".";
        exit(EXIT_FAILURE);
    }
    
    //Inizialize LCP_threshold and LCP_minima
    LCP_threshold = vector<bool>(len,false);
    LCP_minima = vector<bool>(len,false);
    
    //Read LCP array
    dataTypelenSeq *bufferLCP = new dataTypelenSeq[BUFFERSIZE+1];
    bufferLCP[0]=0;
    dataTypeNChar numcharLCP, toRead = BUFFERSIZE, size = 0;

    fseek(InLCP,0,SEEK_SET);
    
    if(len < toRead) toRead = len;
    
    numcharLCP=fread(&bufferLCP[1],sizeof(dataTypelenSeq), toRead,InLCP);
    len -= numcharLCP;
    
    bool decreasing = false;
    
    while(numcharLCP>0){
        for (dataTypeNChar i=1; i <= numcharLCP; i++){
            //Set LCP_threshold: LCP_threshold[i] == 1 iff LCP[i] >= K
            if(bufferLCP[i] >= K)
                LCP_threshold[size] = 1;
            
            //Set LCP_minima: LCP_minima[i] == 1 iff LCP[i] is a local minimum
	    if( decreasing && (bufferLCP[i-1] <= bufferLCP[i]) ){
                decreasing = false;
                LCP_minima[size-1] = 1; //we have a local minimum
            }
            else if( bufferLCP[i-1] > bufferLCP[i] )
                decreasing=true;
                
            size++;
        }
        
        #if DEBUG
            for (dataTypeNChar i=1; i <= numcharLCP; i++)
                cout << "i: " << size-numcharLCP+i-1 << "\tLCP: " << (int)bufferLCP[i] << "\tLCP_threshold: " << LCP_threshold[size-numcharLCP+i-1] << "\tLCP_minima: " << LCP_minima[size-numcharLCP+i-1] << endl;
        #endif
        
        bufferLCP[0] = bufferLCP[numcharLCP];
        
        if( len < toRead) toRead = len;
        
        numcharLCP=fread(&bufferLCP[1],sizeof(dataTypelenSeq), toRead,InLCP);
        len -=numcharLCP;
        
    }
    assert (len==0);
	
    fclose(InLCP);
    delete [] bufferLCP;
    
}


void buildLF_Cluster(uchar *LF_cluster, dataTypedimAlpha currentPile, dataTypeNChar begin, dataTypeNChar len){
    
    dataTypeNChar *counters = new dataTypeNChar[BCRdec->sizeAlpha];
    for(dataTypedimAlpha i = 0; i<BCRdec->sizeAlpha; i++)
        counters[i]=0;
    
    dataTypeNChar posToRead = fVect[currentPile];
    
    //Open BWT file
    FILE *InBWT = fopen(input_dna.c_str(), "rb");
    if (InBWT==NULL) {
        std::cerr << "Error opening " << input_dna << std::endl;
        exit (EXIT_FAILURE);
    }
    
    //1 - Find the block where the cluster starts and count the number of symbols appearing in it up to begin-1
    dataTypeNChar numBlock = 0;
    dataTypeNChar rem = begin;
    
    //findBlockToRead updates counters using vectorOcc to know how many symbols are in currentPile-BWT
    if(begin>0)
        BCRdec->findBlockToRead(counters,currentPile,&rem,&numBlock);
    
    //updates counters using tableOcc to know how many symbols are in the previous currentPile-BWTs
    for(dataTypedimAlpha i = 0; i<BCRdec->sizeAlpha; i++){
        for(dataTypedimAlpha g = 0; g < currentPile; g++)
            counters[i]+=BCRdec->tableOcc[g][i];
    }
    //update pos to read BWT
    posToRead += numBlock*DIMBLOCK;
    
    #if DEBUG
    cout << "posToRead: " << posToRead << "\tlen: " << len << "\tcurrentPile: " << (int)currentPile << "\tnumBlock: " << numBlock << "\trem: " << rem << endl;
    #endif
    
    //Read BWT
    fseek(InBWT, posToRead*sizeof(uchar), SEEK_SET);
    dataTypedimAlpha *bufferBWT = new uchar[DIMBLOCK+len];
    dataTypeNChar numcharBWT = fread(bufferBWT,sizeof(uchar), DIMBLOCK+len, InBWT);
    
    dataTypeNChar i;
    
    for (i = 0; i < rem; i++) //update counters
        counters[BCRdec->alpha[bufferBWT[i]]]++;
    
    //2 - Compute the position in F --> find symbol in L
    for (dataTypeNChar j = 0; j < len; j++)
    {
        posToRead = fVect[BCRdec->alpha[bufferBWT[i]]] + counters[BCRdec->alpha[bufferBWT[i]]];
        fseek(InBWT, posToRead*sizeof(uchar), SEEK_SET);
        numcharBWT = fread(LF_cluster +j,sizeof(uchar), 1, InBWT);
        assert(numcharBWT == 1);
        
        counters[BCRdec->alpha[bufferBWT[i]]]++;
        i++;
        
    }
	
    fclose(InBWT);
    delete [] bufferBWT;
    delete [] counters;
}

//computing the average quality score in a cluster
dataTypedimAlpha avg_qs(dataTypedimAlpha *BWT, dataTypedimAlpha *QUAL, dataTypeNChar len){

    int sum=0;
    dataTypeNChar num=0;

    for(dataTypeNChar j=0; j< len; j++){

        if(BWT[j] != TERM){

            sum=sum+QUAL[j];
            num++;
        }


    }

    if(sum==0) return 0;
    return (dataTypedimAlpha)round((float)sum/num);

}


//computing the max quality score in a cluster
dataTypedimAlpha max_qs(dataTypedimAlpha *BWT, dataTypedimAlpha *QUAL, dataTypeNChar len){

    dataTypedimAlpha max = 0;
    for(dataTypeNChar j=0; j < len; j++){

        if(BWT[j] != TERM){

            if( QUAL[j] > max){
                max = QUAL[j];
            }
        }


    }
    return max;

}


//computing the mean error 
dataTypedimAlpha mean_error(dataTypedimAlpha *BWT, dataTypedimAlpha *QUAL, dataTypeNChar len){

    double avg_err = 0.0;
    dataTypeNChar num = 0;
    double sum_err = 0.0;

    for(dataTypeNChar j=0; j< len; j++){

        if(BWT[j] != TERM){
            num++;
            sum_err = sum_err + pow(10, -((double)QUAL[j]-33.0)/10.0);
        }

    }

    avg_err = sum_err/num;
    
    dataTypedimAlpha qs = (dataTypedimAlpha)round(-10*log10(avg_err));
    
    return qs+33;
    
}

//if mod == 0, modBasesSmoothQS just substitutes QS symbols
void modBasesSmoothQS(dataTypedimAlpha *BWT, dataTypedimAlpha *QUAL, dataTypedimAlpha i, dataTypeNChar begin, dataTypeNChar len, dataTypedimAlpha newSymb, dataTypedimAlpha newqs){
    for(dataTypeNChar j = 0; j < len; ++j)
    {
        if( BWT[j] != TERM ){
            if ( (BWT[j] != newSymb) and (lowQS[BCRdec->alpha[BWT[j]]]==0) )
            {
                #if DEBUG
                    cout << "j: " << j << "\tBWT: " << BWT[j] << "\tBWT_MOD: " << newSymb << endl;
                #endif
                BCRdec->BWT_MOD[i].push_back(newSymb);
                rankbv_setbit(BCRdec->rbv[i],begin+j);
                modified++;
            }
            else if(BWT[j] == newSymb){
		QUAL[j] = newqs;
		qs_smoothed++;
  	    }
	    else //(BWT[j] != newSymb) and (lowQS[ord(BWT[j])]==1)
	    {	
		if (newqs < QUAL[j])//only high quality scores are smoothed
		{
			QUAL[j] = newqs;    
			qs_smoothed++;
		} 
            }
        }
    }//end-for
}



/*
 * START PROCEDURE TO ANALYZE A CLUSTER [begin, len]
 */

dataTypeNChar process_cluster(FILE *InBWT, FILE *InQS, dataTypedimAlpha currentPile, dataTypeNChar start, dataTypeNChar len){

    //cluster is too short
    if((int)len < m) return 0;
    
    dataTypedimAlpha newqs;
  
    //Load bwt and qs
    uchar *bwt = new uchar[len];
    uchar *qs = new uchar[len];
    
    dataTypeNChar numcharBWT, numcharQS;
    
    fseek(InBWT,start*sizeof(uchar),SEEK_SET);
    numcharBWT = fread(bwt,sizeof(uchar),len,InBWT);
    fseek(InQS,start*sizeof(uchar),SEEK_SET);
    numcharQS = fread(qs,sizeof(uchar),len,InQS);
    assert(numcharQS==numcharBWT);
    
    #if DEBUG
        //if(verbose) cout << "----\n";
        //printing bases+QS in the cluster to look them up
        cout << "----\n";
    #endif
    
    dataTypeNChar base_num = 0;
	
    for(dataTypedimAlpha i=0; i<BCRdec->sizeAlpha; i++)
	    lowQS[i]=0;
	
    //freqs is a vector of length sizeAlpha inizialized to 0
    for(dataTypeNChar j = 0; j < len; ++j)
    {
        #if DEBUG
            cout << bwt[j] << "\t" << (int)qs[j]-33 << endl;
        #endif
        
        // Counts the frequency of each symbol (different from TERM) and stores it in a vector
        if(bwt[j] != TERM){
            freqs[BCRdec->alpha[bwt[j]]]++;
            base_num++;
	    if(qs[j] >= quality_threshold + 33)
		lowQS[BCRdec->alpha[bwt[j]]]=1;
        }

        //if(verbose) cout << bwt[j] << "\t" << (int)QUAL[j]-33 << endl;
    }
    #if DEBUG
        cout << "----\n";
    #endif
    
    //Only end-markers --> do nothing
    if(base_num == 0) return 0;
	bases_inside+=base_num;

   /* Set newqs 
	
	max_qs finds the highest qs in the cluster, 
	mean_error finds the mean of error rates, 
	avg_qs finds the average qs in the cluster
   */
	
    #if M==0
	newqs = max_qs(bwt,qs,len);
    #elif M==1
	newqs = mean_error(bwt,qs,len);
    #elif M==2
	newqs = default_value;
    #elif M==3
	newqs = avg_qs(bwt,qs,len);
    #else
	cout << "WARNING: unsupported choice. The process will use M=0." << endl;
	newqs = max_qs(bwt,qs,len);
    #endif
    
    //if(verbose) cout << "****\n";

    
    //a frequent symbol has frequency percentage greater than rare_threshold
    vector < dataTypedimAlpha > FreqSymb;
    #if DEBUG
        cout << "Symbol\t perc" << endl;
    #endif
    
    for(dataTypedimAlpha i=0; i<BCRdec->sizeAlpha; i++){
        if(freqs[i]>0){
            uchar perc = (100*freqs[i])/(base_num); //integer division
            #if DEBUG
                fi (verbose) cout << BCRdec->alphaInverse[i] << "\t" << (int)perc << endl;
            #endif
            if( perc >= freq_threshold )
                FreqSymb.push_back(BCRdec->alphaInverse[i]);
            
            //reset freqs
            freqs[i]=0;
        }
    }
    #if DEBUG
        if (verbose) cout << "FreqSymb.size: " << FreqSymb.size() << "\tbase_num: " << base_num << endl;
    #endif
    
    //We expect  to  have  no  more  than  two  frequent  symbols  in  any  cluster
    assert(FreqSymb.size() < 3 );
    
    /* Noise reduction and QS smoothing */
    if (FreqSymb.size() == 0)    //--> no information to modify clusters
        num_clust_discarded++;
    else if(FreqSymb.size()==1)
    {    //There is a unique most frequent symbol --> modify bases according to it
        if( (FreqSymb[0] == 'N') ) //it is 'N' --> no info
            num_clust_discarded++;
        else{
            modBasesSmoothQS(bwt,qs,currentPile,start,len,FreqSymb[0],newqs);
	    num_clust_mod++;
		}	
    }
    else if(base_num < (uint)m) //There are less than m bases in the cluster and FreqSymb.size() == 2 --> no information to modify bases
        num_clust_discarded++;
    else
    {   //FreqSymb.size() == 2 && base_num >= m
        
        //one of them is 'N' --> modify according to the other
        if (FreqSymb[0] == 'N') {
            //FreqSymb[1] cannot be TERM neither 'N'
            modBasesSmoothQS(bwt,qs,currentPile,start,len,FreqSymb[1],newqs);
	    num_clust_mod++;
        }
        else if (FreqSymb[1] == 'N'){
            //FreqSymb[0] cannot be TERM neither 'N'
            modBasesSmoothQS(bwt,qs,currentPile,start,len,FreqSymb[0],newqs);
            num_clust_mod++;
        }
        else //(FreqSymb[0] != TERM) and (FreqSymb[0] != 'N') and (FreqSymb[1] != TERM) and (FreqSymb[1] != 'N')
        {
            //perform modification according to the two most frequent bases
            //1 - cluster LF mapping
            uchar *LF_cluster = new uchar[len];
            buildLF_Cluster(LF_cluster,currentPile,start,len);
            
            #if DEBUG
                cout << "LF_Cluster: " << LF_cluster << endl;
            #endif
            
            //2 - Find the symbol preceding FreqSymb[0] (resp. FreqSymb[1])
            dataTypedimAlpha symbPrec_0 = '\0', symbPrec_1 = '\0';
            vector <dataTypeNChar> freq_0, freq_1;
            freq_0.resize(BCRdec->sizeAlpha,0);
            freq_1.resize(BCRdec->sizeAlpha,0);
            
            for(dataTypeNChar j = 0; j < len; ++j){
                if(bwt[j] == FreqSymb[0]){
		    if(LF_cluster[j]!=TERM && LF_cluster[j]!='N'){
			freq_0[ord(LF_cluster[j])]=1; //ord: A->0,C->1,G->2,T->3,N->4
			symbPrec_0=LF_cluster[j];
		    }
		}
                if(bwt[j] == FreqSymb[1]){
		    if(LF_cluster[j]!=TERM && LF_cluster[j]!='N'){
			freq_1[ord(LF_cluster[j])]=1; //ord: A->0,C->1,G->2,T->3,N->4
			symbPrec_1=LF_cluster[j];
		    }
		}
            }
            for(int i=0;i<4;i++){ //no Ns symbols are candidates
		freq_0[5]+=freq_0[i];
		freq_1[5]+=freq_1[i]; //freqs_[1][5]=sum_i freqs_[1][i]
	    }
		  
	     //3- Modify the cluster only if both predominant symbols are precedeed by a unique DIFFERENT dna symbol
	     if(freq_0[5]==1 && freq_1[5]==1 && symbPrec_0 != symbPrec_1){
            	num_clust_mod++;
                for(dataTypeNChar j = 0; j < len; ++j){
                    if(bwt[j] != TERM){
                        if ( bwt[j]!= FreqSymb[0] and bwt[j]!= FreqSymb[1] and (lowQS[BCRdec->alpha[bwt[j]]]==0) )
                        {
                            //Check if symbol preceding bwt[j] is equal either to symbPrec_0 or to symbPrec_1
                            if(LF_cluster[j]==symbPrec_0){
                            #if DEBUG
                                cout << "j: " << j << "\tBWT: " << bwt[j] << "\tBWT_MOD: " << FreqSymb[0] << endl;
                            #endif
                                rankbv_setbit(BCRdec->rbv[currentPile],start+j);
                                BCRdec->BWT_MOD[currentPile].push_back(FreqSymb[0]);
                                modified++;
                            }
                            else if(LF_cluster[j]==symbPrec_1){
                            #if DEBUG
                                cout << "j: " << j << "\tBWT: " << bwt[j] << "\tBWT_MOD: " << FreqSymb[1] << endl;
                            #endif
                                rankbv_setbit(BCRdec->rbv[currentPile],start+j);
                                BCRdec->BWT_MOD[currentPile].push_back(FreqSymb[1]);
                                modified++;
                            }
                            //else --> no information to modify bases
                            
                        }//end-if modify bases
                        else if(bwt[j]== FreqSymb[0] or bwt[j]== FreqSymb[1]){
				qs[j] = newqs;
				qs_smoothed++;
			}
			else{
				if(newqs < qs[j]){
					qs[j] = newqs;
					qs_smoothed++;
				}
			}
                        
                    }//end-if (bwt[j] != TERM)
                    
                }//end-for

            }
	    else
	    {  //--> no information to modify
		num_clust_amb_discarded++;
	    }
            delete [] LF_cluster;
	    freq_0.clear();
	    freq_0.shrink_to_fit();
	    freq_1.clear();
	    freq_1.shrink_to_fit();
                       
        }//end if-else
                       
    }//end if-else FreqSymb.size() == 2 && base_num >= 5
                   
    FreqSymb.clear();
    FreqSymb.shrink_to_fit();
	
    //Write QS smoothing on partial QS files
    fseek(InQS,start*sizeof(uchar),SEEK_SET);
    numcharQS = fwrite(qs,sizeof(uchar),len,InQS);
    
    #if DEBUG
	if (verbose){
    cout << "**QS_MOD**" << endl;
        for(dataTypeNChar j = 0; j < len; ++j)
            cout << bwt[j] << "\t" << (int)qs[j]-33 << "\t" << qs[j] << endl;
    cout << "*****" << endl;}
    #endif
    
    delete [] bwt;
    delete [] qs;
    return 1;
}

/*
 * END PROCEDURE TO ANALYZE A CLUSTER
 */


                   
                   
/*
 * PROCEDURE run DETECTS clusters and EXECUTES process_cluster
 */
                   
dataTypeNChar run(){
    
    dataTypeNChar begin = 0;//begin position
    
    dataTypeNChar clust_len=0;
    dataTypeNChar num_clust=0;
    bases_inside = 0;
    freqs.resize(BCRdec->sizeAlpha,0); //temporary vector used to count frequency of bases inside each cluster
    lowQS.resize(BCRdec->sizeAlpha,0); //bitvector for qs symbols associated with dna symbols in any cluster (if 0 all occ. have low qs)

    
    bool cluster_open=false;

    //Open BWT and QS files
    vector <FILE *> BWT;
    vector <FILE *> BWT_QS;
    
    BWT.resize(BCRdec->sizeAlpha);
    BWT_QS.resize(BCRdec->sizeAlpha);
    
    for(dataTypedimAlpha i = 0; i < BCRdec->sizeAlpha; i++){
        string partialFile = "tmp/bwt_" + to_string(i)+".aux";
        BWT[i] = fopen(partialFile.c_str(), "rb");
        if (BWT[i]==NULL){
            cerr << "Error opening " << partialFile << ".";
            exit(EXIT_FAILURE);
        }
        fseek(BWT[i],0,SEEK_SET);
        
        partialFile = "tmp/bwt_qs_" + to_string(i)+".aux";
        BWT_QS[i] = fopen(partialFile.c_str(), "r+b");
        if (BWT_QS[i]==NULL){
            cerr << "Error opening " << partialFile << ".";
            exit(EXIT_FAILURE);
        }
        fseek(BWT_QS[i],0,SEEK_SET);
    }
    
    
    //procedure that identifies clusters by looking at LCP_threshold and LCP_minima
    for(dataTypeNChar i=1; i < BCRdec->n; ++i){
        
        // 1 and 0
        if(LCP_threshold[i] and not LCP_minima[i] ) {
            
            if(cluster_open){//extend current cluster
                
                clust_len++;
                
            }else
            {//open new cluster
                
                cluster_open=true;
                begin=i-1;
                
                clust_len=2; //pos i and i-1
                
            }
            
        }
        else //0 and 1, 0 and 0, 1 and 1
        {
            
            if(cluster_open){//close current cluster
                
                clusters_size += clust_len;
                //clust_len is without current position i
                
                //find the pile the cluster belongs to
                dataTypedimAlpha j = (dataTypedimAlpha)(((std::upper_bound (fVect.begin(), fVect.end(), i-1))-fVect.begin())-1);
            #if DEBUG
                if (verbose) cout << "Cluster " << num_clust << " in Pile " << (int)j << " = [ " << begin-fVect[j] << "," << clust_len << " ]" << endl;
            #endif
                
                num_clust+=process_cluster(BWT[j], BWT_QS[j], j ,begin-fVect[j], clust_len);

            }
            
            cluster_open=false;
            clust_len = 0;
            
        }
        
    }
    //Last cluster
    if(cluster_open){
        
        clusters_size += clust_len;
        //clust_len is without position i

        dataTypedimAlpha j = BCRdec->sizeAlpha -1;
        #if DEBUG
        if (verbose) cout << "Cluster " << num_clust << " in Pile " << (int)j << " = [ " << begin-fVect[j] << "," << clust_len << " ]" << endl;
        #endif
        
        num_clust+=process_cluster(BWT[j], BWT_QS[j], j, begin - fVect[j], clust_len);
        
    }
    
    //Close BWT and QS files
    for(dataTypedimAlpha i = 0; i < BCRdec->sizeAlpha; i++){
        fclose(BWT[i]);
        fclose(BWT_QS[i]);
    }
    
    LCP_threshold.clear();
    LCP_threshold.shrink_to_fit();
    LCP_minima.clear();
    LCP_minima.shrink_to_fit();
	
    return num_clust;
}

/*
 * END run
 */


bool file_exists(string fileName)
{
    std::ifstream infile(fileName);
    return infile.good();
}

int main(int argc, char** argv){
    
    srand(time(NULL));

    
    if(argc < 3) help();
    
    int opt;
    while ((opt = getopt(argc, argv, "he:q:o:a:l:k:m:v:f:t:s:T:H:DV")) != -1){
        switch (opt){
            case 'h':
                help();
                break;
            case 'e':
                input_dna = string(optarg);
                break;
            case 'q':
                input_qual = string(optarg);
                break;
            case 'o':
                output = string(optarg);
                break;
            case 'a':
                input_lcp = string(optarg);
                break;
            case 'l':
                maxLengthRead = atoi(optarg);
                break;
            case 'k':
                K = atoi(optarg);
                break;
            case 'm':
                m = atoi(optarg);
                break;
	    case 'v':
                default_value = atoi(optarg);
                break;
	    case 't':
		quality_threshold = atoi(optarg);
		break;
	    case 'f':
		freq_threshold = atoi(optarg);
		break;
            case 's':
                TERM = atoi(optarg);
                break;
            case 'T':
                num_threads = atoi(optarg);
                break;
            case 'H':
                input_titles = string(optarg);
                isFileHeaders = true;
                break;
            case 'D':
                debug=true;
                break;
            case 'V':
            verbose=true;
            break;
            default:
                help();
                return -1;
        }
    }
    
    K = K == -1 ? K_def : K;
    m = m == 0 ? m_def : m;
    default_value = default_value == '\0' ? default_value_def : default_value;
    num_threads = num_threads == 0 ? num_threads_def : num_threads;
    quality_threshold = quality_threshold == -1 ? quality_threshold_def : quality_threshold;
    freq_threshold = freq_threshold == 0.0 ? freq_threshold_def : freq_threshold;
    
    if( (input_dna.compare("")==0) || (input_qual.compare("")==0) || (output.compare("")==0) ||(input_lcp.compare("")==0) || (maxLengthRead==0 ) || (isFileHeaders && input_titles.compare("")==0)
       ) help();
    
    if(not file_exists(input_dna)){
        cout << "Error: could not find file " << input_dna << "." << endl << endl;
        help();
    }
    
    if(not file_exists(input_qual)){
        cout << "Error: could not find file " << input_qual << endl << endl;
        help();
    }
    
    if(not file_exists(input_lcp)){
        cout << "Error: could not find file " << input_lcp << endl << endl;
        help();
    }
    
    if(isFileHeaders and not file_exists(input_titles)){
        cout << "Error: could not find file " << input_titles << "." << endl << endl;
        help();
    }
        
    cout << "Running bfq_ext..." << endl;
	cout << "\tMode: " << M;
	#if M==2
		cout << "\treplacing QS with symbol: " <<  default_value ;
	#endif
	cout << "\n\tK: " << K << endl;
	cout << "\tm: " << m << endl;
	cout << "\tFrequency threshold: " << freq_threshold << "\%" << endl << endl;
	cout << "Output fastq file: " << output << endl;
	cout << endl;
    
    BCRdec = new BCRdecode(num_threads,isFileHeaders);
    BCRdec->TERMINATE_CHAR=TERM;

    cout << "\n\nPhase 1/5: Split Files..." << endl;
    buildFreq();
    splitIntoPartial();
    cout << "done." << endl;
    /**/
    fprintf(stdout, "##\nmalloc_count ### current peak: %'zu\n##\n", malloc_count_peak_curr());
    /**/
    
    cout << "\n\nPhase 2/5: Compute VectorUnbuildBCR..." << endl;
    BCRdec->computeVectorUnbuildBCR();
    cout << "done." << endl;
    /**/
    fprintf(stdout, "##\nmalloc_count ### current peak: %'zu\n##\n", malloc_count_peak_curr());
    /**/
    
    cout << "\n\nPhase 3/5: Reading LCP array ..." << endl;
    buildBitVectors(BCRdec->n);
    cout << "done." << endl;
    /**/
    fprintf(stdout, "##\nmalloc_count ### current peak: %'zu\n##\n", malloc_count_peak_curr());
    /**/
    
    //start procedure run
    cout << "\n\nPhase 4/5: Process clusters ..." << endl;
    dataTypeNChar tot_clust = run();
    
    for(dataTypedimAlpha i=0; i < BCRdec->sizeAlpha; i++)
        rankbv_build(BCRdec->rbv[i]);
    
    cout << "done." << endl;
    /**/
    fprintf(stdout, "##\nmalloc_count ### current peak: %'zu\n##\n", malloc_count_peak_curr());
    /**/
    
    cout << "\n\nPhase 5/5: Inverting eBWT ..." << endl;

    BCRdec->BCRInverse(input_titles, output, maxLengthRead);
    
    std::cerr << "\nThe rebuilt file is ready!\n\n";
    /**/
    fprintf(stdout, "##\nmalloc_count ### current peak: %'zu\n##\n", malloc_count_peak_curr());
    /**/
    
    cout << "**** Cluster statistics ****" << endl;
    cout << tot_clust << " is the total number of clusters" << endl;
    cout << "Discarded: " << num_clust_discarded << "(" << (double(100*num_clust_discarded)/(tot_clust)) <<  "%)" << endl;
    cout << "Ambiguous discarded: " << num_clust_amb_discarded << "(" << (double(100*num_clust_amb_discarded)/(tot_clust)) <<  "%)" << endl;
    cout << "Processed: " << num_clust_mod << "(" << (double(100*num_clust_mod)/(tot_clust)) <<  "%)" << endl << endl;
    cout << bases_inside << "/" << (BCRdec->n - BCRdec->nText) << " (" << (double(100*bases_inside)/(BCRdec->n - BCRdec->nText)) <<  "%) bases fall inside clusters" << endl;   
	
    cout << "**** Quality statistics ****" << endl;
    cout << qs_smoothed << "/" << (BCRdec->n - BCRdec->nText) << " qualities have been modified (" << 100*double(qs_smoothed)/(BCRdec->n - BCRdec->nText) << "% of all qs and " <<
		100*double(qs_smoothed)/bases_inside << "% of qs inside clusters)." << endl;
    cout << "***********************" << endl << endl;
  
    cout << "**** Bases statistics ****" << endl;
    cout << modified << "/" << BCRdec->n - BCRdec->nText << " bases have been modified (" << 100*double(modified)/(BCRdec->n - BCRdec->nText) << "% of all bases and " <<
    100*double(modified)/bases_inside << "% of bases inside clusters)." << endl;
    cout << "***********************" << endl << endl;

    for(dataTypedimAlpha i=0; i < BCRdec->sizeAlpha; i++)
        rankbv_free(BCRdec->rbv[i]);
    
    delete BCRdec;
    
    return 0;
}
