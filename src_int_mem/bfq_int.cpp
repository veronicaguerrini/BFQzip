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
#include <cstring>
#include "../external/bwt2lcp/include.hpp"
#include "../external/bwt2lcp/dna_string_n.hpp"
#include "../external/bwt2lcp/dna_bwt_n.hpp"
#include "../external/rankbv/rankbv.hpp"
#include "../external/malloc_count/malloc_count.h"

#ifndef DEBUG
  #define DEBUG 1
#endif


#define LONGEST 10000 //longest read

using namespace std;

string input_dna;
string input_qual;
string input_titles;
string output;

bool ignore_headers = true; //ignore headers

/*
 * Debug mode variables: print the BWT and read names/qualities for each base
 */

bool debug = false; //print debug info
bool verbose = false; //verbose output
int max_id_len=20; //in read_info, store at most this number of chars for the IDs
vector<string> read_info;//if debug, store read coordinate for every BWT position
//------------

vector<string> read_ids;//ID of each read

uint64_t modified = 0;//count how many bases have been modified
uint64_t qs_smoothed; //count how many qs have been modified

uint64_t bases_inside; //total number of bases inside analyzed clusters

uint64_t num_clust; //total number of clusters
uint64_t num_clust_discarded; //total number of cluster discareded for modifications
uint64_t num_clust_amb_discarded; //total number of ambiguous cluster discareded for modifications
uint64_t num_clust_mod; //total number of cluster modified
uint64_t num_clust_alleq; //total number of cluster with only one symbol

vector<uint64_t> statistics_qual_before(256,0);//count absolute frequencies of qualities in reads, before modifying
vector<uint64_t> statistics_qual_after(256,0);//count absolute frequencies of qualities in reads, after modifying

int border = 1;//exclude/include this number of bases at the borders of the analyzed cluster

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

//terminator character at the end of the reads
char TERM = '#';

char *QUAL;//string of length |BWT| that contains the base qualities, for each BWT position

string BWT_MOD;
rankbv_t* rbv = NULL;

vector<bool> LCP_minima;//bitvector that stores the LCP minima
vector<bool> LCP_threshold;//bitvector that stores LCP values that exceed the threshold: >= K

dna_bwt_n_t bwt;//the BWT data structure

uint64_t freqs[5]{0}; //#occ. of dna symbols in any cluster
uint64_t lowQS[5]{0}; //bitvector for qs symbols associated with dna symbols in any cluster (if 0 all occ. have low qs)


void help(){
    
    cout << "bfq_int [options]" << endl <<
    "Options:" << endl <<
    "-h          Print this help." << endl <<
    "-e <arg>    Input eBWT file (A,C,G,T,#) of DNA (REQUIRED)." << endl <<
    "-q <arg>    Qualities permuted according to the DNA's ebwt (REQUIRED)." << endl <<
    "-o <arg>    Output fastq (REQUIRED)." << endl <<
	"-k <arg>    Minimum LCP required in clusters. Default: " << K_def << "." << endl <<
    "-m <arg>    Minimum length of cluster to be processed. Default: " << m_def << "." << endl <<
    "-v <arg>    Quality score for constant replacement (if M=2). Default: " << (int)default_value_def-33 << "." << endl <<
	"-f <arg>    Percentage threshold for frequent bases in clusters. Default: " << freq_threshold_def << "." << endl <<
	"-t <arg>    Quality score threshold for trusted bases. Default: " << quality_threshold_def << "." << endl <<
    "-s <arg>    ASCII value of terminator character. Default: " << int('#') << " (#)." << endl <<
    "-H <arg>    List of original headers." << endl << 
	"-D          Print debug info for each BWT position." << endl << endl <<
    
    "\nTo run bfq_int, you must first build the extended Burrows-Wheeler Transform " <<
    "of the input DNA sequences and the corresponding permutation of base quality scores." << endl;
    
    exit(0);
}

/*
 *  START PROCEDURES TO NAVIGATE SUFFIX TREE AND COMPUTE LCP MINIMA
 */

void update_LCP_leaf(sa_leaf L, uint64_t & lcp_values){
    
  for(uint64_t i = L.rn.first+1; i<L.rn.second; ++i){ 
    LCP_threshold[i] = (L.depth >= K);
    lcp_values++;
  }
}

void update_lcp_minima(sa_node_n x, uint64_t & n_min){
    
  /*
   * we have a minimum after the end of each child (that is different than #) of size at least 2 of the input node x, except
   * if the candidate minimum position is the last or exceeds the interval of x
   */
  
  if( x.first_C - x.first_A >= 2 and     // there are at least 2 'A'
     x.first_C < x.last-1             // candidate min in x.first_C is not >= last position
     ){
      LCP_minima[x.first_C] = true;
      n_min++;
  }
  
  if( x.first_G - x.first_C >= 2 and     // there are at least 2 'C'
     x.first_G < x.last-1             // candidate min in x.first_G is not >= last position
     ){
      LCP_minima[x.first_G] = true;
      n_min++;
  }

  if( x.first_N - x.first_G >= 2 and     // there are at least 2 'G'
     x.first_N < x.last-1             // candidate min in x.first_N is not >= last position
     ){
      LCP_minima[x.first_N] = true;
      n_min++;
  }

  if( x.first_T - x.first_N >= 2 and     // there are at least 2 'N'
     x.first_T < x.last-1             // candidate min in x.first_T is not >= last position
     ){
      LCP_minima[x.first_T] = true;
      n_min++;
  }
}

void detect_minima(){
    
  uint64_t n = bwt.size();
  
  cout << "\nPhase 2/5: navigating suffix tree leaves." << endl;
  
  /*
   * LCP_threshold[i] == 1 iff LCP[i] >= K
   */
  LCP_threshold = vector<bool>(n,false);
  
  uint64_t leaves = 0;//number of visited leaves
  uint64_t max_stack = 0;
  uint64_t lcp_values = 1;//number of computed LCP values
  
  {
    auto TMP_LEAVES = vector<sa_leaf>(5);
    
    stack<sa_leaf> S;
    S.push(bwt.first_leaf());
    
    int last_perc_lcp = -1;
    int perc_lcp = 0;
    
    while(not S.empty()){
        
      sa_leaf L = S.top();
      S.pop();
      leaves++;
      
      assert(leaf_size(L)>0);
      max_stack = S.size() > max_stack ? S.size() : max_stack;
      
      update_LCP_leaf(L,lcp_values);
      
      int t = 0;//number of children leaves
      bwt.next_leaves(L, TMP_LEAVES, t, 2);
      
      for(int i=t-1;i>=0;--i) S.push(TMP_LEAVES[i]);
      perc_lcp = (100*lcp_values)/n;
      
      if(perc_lcp > last_perc_lcp){
        #if DEBUG
          if(verbose) cout << "LCP: " << perc_lcp << "%." <<  endl;
        #endif  
        last_perc_lcp = perc_lcp;
      }
    }
  }
  
  cout << "Computed " << lcp_values << "/" << n << " LCP threshold values." << endl;
  
  cout << "Max stack depth = " << max_stack << endl;
  cout << "Processed " << leaves << " suffix-tree leaves." << endl << endl;
  
  /**/
  fprintf(stdout, "##\nmalloc_count ### current peak: %'zu\n##\n", malloc_count_peak_curr());
  /**/
  
  cout << "\n\nPhase 3/5: computing LCP minima." << endl;
  
  LCP_minima = vector<bool>(n,false);
  //LCP_minima = new bool[n]{false};
  
  auto TMP_NODES = vector<sa_node_n>(5);
  
  uint64_t nodes = 0;//visited ST nodes
  max_stack = 0;
  
  stack<sa_node_n> S;
  S.push(bwt.root());
  
  int last_perc_lcp = -1;
  int perc_lcp = 0;
  uint64_t n_min = 0;//number of LCP minima
  
  while(not S.empty()){
      
      max_stack = S.size() > max_stack ? S.size() : max_stack;
      
      sa_node_n N = S.top();
      S.pop();
      nodes++;
      
      //compute LCP values at the borders of N's children
      update_lcp_threshold(N, LCP_threshold, lcp_values, K);
      
      update_lcp_minima(N, n_min);
      
      //follow Weiner links
      int t = 0;
      bwt.next_nodes(N, TMP_NODES, t);
      
      for(int i=t-1;i>=0;--i) S.push(TMP_NODES[i]);
      
      perc_lcp = (100*lcp_values)/n;
      
      if(perc_lcp > last_perc_lcp){
          
          #if DEBUG
            if(verbose) cout << "LCP: " << perc_lcp << "%." << endl;
          #endif
          
          last_perc_lcp = perc_lcp;
          
      }
      
  }
  
  cout << "Computed " << lcp_values << "/" << n << " LCP values." << endl;
  cout << "Found " << n_min << " LCP minima." << endl;
  cout << "Max stack depth = " << max_stack << endl;
  cout << "Processed " << nodes << " suffix-tree nodes." << endl;
  
  /**/
  fprintf(stdout, "##\nmalloc_count ### current peak: %'zu\n##\n", malloc_count_peak_curr());
  /**/
}

/*
 *  END PROCEDURES TO NAVIGATE SUFFIX TREE AND COMPUTE LCP MINIMA
 */

//return a quality score according to Illumina 8-Level Binning
int illumina_8_level_binning(int newqs){

  if (newqs >= 40) newqs = 40;
  else if (newqs >= 35) newqs = 37;
  else if (newqs >= 30) newqs = 33;
  else if (newqs >= 25) newqs = 27;
  else if (newqs >= 20) newqs = 22;
  else if (newqs >= 10) newqs = 15;
  else if (newqs >= 2) newqs = 6;

  return newqs+33;

}


//computing the average quality score in a cluster
int avg_qs(uint64_t start, uint64_t end){

  int sum=0;
  uint64_t num=0;
  
  for(uint64_t j=start; j<=end; j++){
  
  	if(bwt[j] != TERM){
  		sum=sum+(int)QUAL[j];
  		num++;
  	}
  }

  if(sum==0) return 0;
  return (sum/num);
}


//computing the max quality score in a cluster
char max_qs(uint64_t start, uint64_t end){

  char max=0;
  for(uint64_t j=start; j<=end; j++){
    if(bwt[j] != TERM){
      if(QUAL[j] > max){
        max = QUAL[j];
      }
    }
  }
  return max;
}


//computing the mean error
int mean_error(uint64_t start, uint64_t end){

double avg_err=0;
uint64_t num = 0;
double sum_err = 0;

  for(uint64_t j=start; j<=end; j++){
    if(bwt[j] != TERM){
      num++;
      sum_err = sum_err + pow(10, -((double)QUAL[j]-33)/10);
    }
  }
  avg_err = sum_err/num;
  int qs = round(-10*log10(avg_err));

return qs+33;
}

//analysing clusters
void modBasesSmoothQS(uint64_t begin, uint64_t end, char newSymb, char newqs){
    for(uint64_t j = begin; j <= end; ++j)
    {
        if( bwt[j] != TERM ){
            if ((bwt[j] != newSymb) and (lowQS[ord(bwt[j])]==0))
            {
                #if DEBUG
                    if (verbose) cout << "j: " << j << "\tBWT: " << bwt[j] << "\tBWT_MOD: " << newSymb << endl;
                #endif
                
                BWT_MOD.push_back(newSymb);
                rankbv_setbit(rbv,j);
                modified++;
            }
            else if(bwt[j] == newSymb){
				QUAL[j] = newqs;
				qs_smoothed++;
			}
			else //(bwt[j] != newSymb) and (lowQS[ord(bwt[j])]==1)
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
 *
 * START PROCEDURE TO ANALYZE A CLUSTER [begin, i]
 *
 * We may include/exclude some symbols at cluster beginning by changing the variable border
 */

uint64_t process_cluster(uint64_t begin, uint64_t i){

  uint64_t start=(begin>=border?begin-border:0);
  uint64_t end=(i>border?i-border:0);
  
  uint64_t size = (end-start+1);

  //cluster is too short
  if(size < m) return size;

  char newqs;
 
  //printing bases+QS in the cluster to look them up
  #if DEBUG
    if(verbose) cout << "----\n";
  #endif

  uint64_t base_num = 0;//number of bases inside cluster
  for(int i=0; i<5; i++){
	freqs[i]=0;
	lowQS[i]=0;
  }
  
  for(uint64_t j = start; j <= end; ++j){
    /*Counts the frequency of each base and its QS */
    if(bwt[j] != TERM){
      freqs[ord(bwt[j])]++;
      base_num++;

	  if(QUAL[j] >= quality_threshold + 33)
		lowQS[ord(bwt[j])]=1;
    }
    #if DEBUG
      if(verbose) cout << bwt[j] << "\t" << (int)QUAL[j]-33 << endl;
    #endif
  }
  
	num_clust++;
	
  if(base_num == 0) return size;
	
	bases_inside+=base_num;
	
  /*
  max_qs finds the highest qs in the cluster, 
  mean_error finds the mean of error rates, 
  avg_qs finds the average qs in the cluster
  */
  #if M==0
    newqs = max_qs(start,end);
  #elif M==1
    newqs = (char)mean_error(start,end);
  #elif M==2
    newqs = default_value;
  #elif M==3
    newqs = (char)avg_qs(start,end);
  #else
    cout << "WARNING: unsupported choice. The process will use M=0." << endl;
    newqs = max_qs(start,end);
  #endif

  #if DEBUG
    if(verbose) cout << "****\n";
  #endif

  //a frequent symbol has frequency percentage greater than freq_threshold
  vector < char > FreqSymb;
  #if DEBUG
     if(verbose) cout << "Symbol\t perc" << endl;
  #endif

	unsigned char nnn = 0;
	
  for(int i=0; i<5; i++){
      if(freqs[i]>0){
		  nnn++;
         unsigned char perc = (100*freqs[i])/(base_num); //integer division
         #if DEBUG
             if(verbose) cout << dna(i) << "\t" << (int)perc << endl;
		#endif
         if( perc >= freq_threshold )
            FreqSymb.push_back(dna(i));
     }
  }
  if(nnn==1)
		num_clust_alleq++;
  #if DEBUG
      if(verbose) cout << "FreqSymb.size: " << FreqSymb.size() << "\tbase_num: " << base_num << endl;	
  #endif

  //We expect  to  have  no  more  than  two  frequent  symbols  in  any  cluster
  assert(FreqSymb.size() < 3 );
	
  /* Noise reduction and QS smoothing */
  if (FreqSymb.size() == 0)
  { //--> no information to modify bases 
	num_clust_discarded++;
  }
  else if(FreqSymb.size()==1)
  {    //There is a unique most frequent symbol --> modify bases according to it 
     if( (FreqSymb[0] == 'N') ) //it is 'N' 
		num_clust_discarded++;
     else{
		  modBasesSmoothQS(start,end,FreqSymb[0],newqs);
	 }
  }
  else if(base_num < m)
  {    //There are less than 5 bases in the cluster and FreqSymb.size() == 2 
		num_clust_discarded++;
  }
  else
  {  //FreqSymb.size() == 2 && base_num >= m
	
	//one of them is 'N' --> modify according to the other
        if (FreqSymb[0] == 'N') {
            //FreqSymb[1] cannot be TERM neither 'N'
            modBasesSmoothQS(start,end,FreqSymb[1],newqs);
			num_clust_mod++;
        }
        else if (FreqSymb[1] == 'N'){
            //FreqSymb[0] cannot be TERM neither 'N'
            modBasesSmoothQS(start,end,FreqSymb[0],newqs);
			num_clust_mod++;
        }
        else //(FreqSymb[0] != TERM) and (FreqSymb[0] != 'N') and (FreqSymb[1] != TERM) and (FreqSymb[1] != 'N')
        {  
	     //perform modification according to the two most frequent bases
	     //find the symbols preceding FreqSymb[0] (resp. FreqSymb[1])
	     char c, symbPrec_[2];
	     char freqs_[2][6]{0};

	     for(uint64_t j = start; j <= end; ++j){
			 if(bwt[j] == FreqSymb[0]){
				 c=bwt[bwt.LF(j)]; //bwt[bwt.LF(j)] can be TERM
				 if(c!=TERM){
					freqs_[0][ord(bwt[bwt.LF(j)])]=1; //ord: A->0,C->1,G->2,T->3,N->4
					symbPrec_[0]=c;
				 }
			 }
			 else if(bwt[j] == FreqSymb[1]){
				 c=bwt[bwt.LF(j)]; //bwt[bwt.LF(j)] can be TERM
				 if(c!=TERM){
					freqs_[1][ord(bwt[bwt.LF(j)])]=1; //ord: A->0,C->1,G->2,T->3,N->4
					symbPrec_[1]=c;
				 }
			 }
	     }   
	
		 for(int i=0;i<4;i++){ //no Ns symbols are candidates
			freqs_[0][5]+=freqs_[0][i];
			freqs_[1][5]+=freqs_[1][i]; //freqs_[1][5]=sum_i freqs_[1][i]
		 }
		  
		 //Modify the cluster only if both predominant symbols are precedeed by a unique DIFFERENT dna symbol
		 if(freqs_[0][5]==1 && freqs_[1][5]==1 && symbPrec_[0] != symbPrec_[1]){
			
			num_clust_mod++;
			
			for(uint64_t j = start; j <= end; ++j){
				if(bwt[j] != TERM){
					if ( bwt[j]!= FreqSymb[0] and bwt[j]!= FreqSymb[1] and lowQS[ord(bwt[j])]==0)
					{
						//Check if symbol preceding bwt[j] is equal either to symbPrec_0 or to symbPrec_1
						c=bwt[bwt.LF(j)];
						if(c==symbPrec_[0]){
						   #if DEBUG
							if(verbose) cout << "j: " << j << "\tBWT: " << bwt[j] << "\tBWT_MOD: " << FreqSymb[0] << endl;
						   #endif
						   rankbv_setbit(rbv,j);
						   BWT_MOD.push_back(FreqSymb[0]);
						   modified++;
						 }
						 else if(c==symbPrec_[1]){
						   #if DEBUG
						    if(verbose) cout << "j: " << j << "\tBWT: " << bwt[j] << "\tBWT_MOD: " << FreqSymb[1] << endl;
						   #endif
						   rankbv_setbit(rbv,j);
						   BWT_MOD.push_back(FreqSymb[1]);
						   modified++;
						 }
			     
					}//end-if modify bases
					else if(bwt[j]== FreqSymb[0] or bwt[j]== FreqSymb[1]){
						QUAL[j] = newqs;
						qs_smoothed++;
					}
					else{
						if(newqs < QUAL[j]){
							QUAL[j] = newqs;
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

		 
		}//end else
	
    }//end if-else FreqSymb.size() == 2 && base_num >= 5
	
	FreqSymb.clear();
	FreqSymb.shrink_to_fit();
	
	return size;
}

/*
 * END PROCEDURE TO ANALYZE A CLUSTER
 */


/*
 * PROCEDURE run NAVIGATES suffix tree, and computes LCP minima, EXECUTES process_cluster for each detected cluster.
 */
void run(){
  
	uint64_t n = bwt.size();
  
	QUAL = (char*) malloc (sizeof(char)*n);
	
	//read base qualities (permuted according to the BWT) into QUAL
	FILE *f_qual;
	
	f_qual = fopen(input_qual.c_str(), "r");
  	if(!f_qual) perror("run");
 
	uint64_t m = fread(QUAL,sizeof(char),n,f_qual);
	assert(m == n);
	
	fclose(f_qual);
  
    #if DEBUG
        for(uint64_t i=0; i<n; i++) statistics_qual_before[QUAL[i]]++;
    #endif

	/**/
	//fprintf(stdout, "##\nmalloc_count ### INSIDE RUN current peak: %'zu\n##\n", malloc_count_peak_curr());
	/**/
	
	//Set bitvector to modify BWT symbols
	{
		rbv = rankbv_create(n,2);
	}	
  
  uint64_t begin = 0;//begin position
  
  uint64_t clust_len=0;
  bool cluster_open=false;

  num_clust = 0;
  num_clust_discarded = 0;
  num_clust_amb_discarded = 0;
  num_clust_mod = 0;
  num_clust_alleq = 0;
  
  #if DEBUG
    //used only to compute and visualize cluster statistics
    uint64_t MAX_CLUST_LEN = 200;
    //auto CLUST_SIZES = vector<uint64_t>(MAX_CLUST_LEN+1,0);
    uint64_t CLUST_SIZES[MAX_CLUST_LEN+1]{0};
  #endif
  
  //procedure that identifies clusters by looking at LCP_threshold and LCP_minima
  for(uint64_t i=0;i<n;++i){
      
    if(LCP_threshold[i] and not LCP_minima[i]){
          
      if(not cluster_open){//open new cluster
        cluster_open=true;
        begin=i;
      }
          
    }else{
          
        if(cluster_open){//close current cluster
		  clust_len = process_cluster(begin, i);//position i included
          #if DEBUG
            if(clust_len<=MAX_CLUST_LEN) CLUST_SIZES[clust_len]++;
          #endif
        }
        cluster_open=false;
      }
  }

  if(cluster_open){//close last cluster
	clust_len = process_cluster(begin, n);//position i included
     #if DEBUG
      if(clust_len<=MAX_CLUST_LEN) CLUST_SIZES[clust_len]++;
     #endif
  }
	

  rankbv_build(rbv);
  

  //Remove bitvectors
  LCP_minima.clear();
  LCP_minima.shrink_to_fit();
  LCP_threshold.clear();
  LCP_threshold.shrink_to_fit();
	
	/**/
	//fprintf(stdout, "##\nmalloc_count ### INSIDE RUN current peak: %'zu\n##\n", malloc_count_peak_curr());
	/**/
	
  //print clusters statistics (i.e. number of bases that fall inside each cluster of a fixed size)
  #if DEBUG
    uint64_t scale = *max_element(&CLUST_SIZES[0], &CLUST_SIZES[MAX_CLUST_LEN]+1);
    for(int i=0;i<=MAX_CLUST_LEN;++i){
      cout << i << ( i < 10 ? "   " : (i<100 ? "  " : " "));
      for(uint64_t j = 0; j < (100*CLUST_SIZES[i])/scale; ++j) cout << "-";
      cout << " " << CLUST_SIZES[i] << endl;
    }
  #endif
}

/*
 * END run
 */

/*
 * PROCEDURE invert INVERT BWT AND WRITE A NEW FASTQ FILE.
 *
 * If the input file consisted of reads and their reverse, we must define the behaviour.
 */
void invert(){
	
	FILE *f_in;
	if(not ignore_headers){
		f_in= fopen(input_titles.c_str(), "r");
		if(!f_in) perror("invert");
	}
	FILE *f_out = fopen(output.c_str(), "w");
	if(!f_out) perror("invert");
	
	char header[]="@\n";
	char plus[]="+\n";
	size_t len = 0;

	char BASES[LONGEST]{0};
	char QS[LONGEST]{0};

	BASES[LONGEST-1] = '\n';
	QS[LONGEST-1] = '\n';


	//number of reads in the file
	uint64_t N = bwt.get_number_of_strings();
  
	string line;
  

	for(uint64_t i = 0;i < N;++i){//for each read 
      
      int nbases = LONGEST-1;
      uint64_t j = i;//bwt[j] = current read character

      while(bwt[j] != TERM){
         
        BASES[--nbases] = (rankbv_access(rbv, j)==1)?(BWT_MOD[rankbv_rank1(rbv,j)-1]):(bwt[j]);

        #if B==1
          QUAL[j] = illumina_8_level_binning((int)QUAL[j]-33);
        #endif

        QS[nbases] = QUAL[j];

        j = bwt.LF(j); //backward search
      }

      #if DEBUG
        for(int i=nbases; i<LONGEST-1; i++) statistics_qual_after[QS[i]]++;
      #endif

      //write output FASTQ file
      char *buf = NULL;
    
      if(not ignore_headers){
        ssize_t size = getline(&buf, &len, f_in); // @'s line
        fwrite(buf, sizeof(char), size, f_out);
      }
      else{
        fwrite(header, sizeof(char), 2, f_out);
      }

      fwrite(&BASES[nbases], sizeof(char), LONGEST-nbases, f_out);
      fwrite(plus, sizeof(char), 2, f_out);
      fwrite(&QS[nbases], sizeof(char), LONGEST-nbases, f_out);

      free(buf);
    }

	if(not ignore_headers)
		fclose(f_in);
	fclose(f_out);

}
/*
 * END invert
 */

/*
 * PROCEDURES TO DEBUG
 */

// for each bwt position, the read coordinate, bwt base, modified base, and modified quality score
void print_info(){
    
  if(debug){
      
    //number of reads
    uint64_t N = bwt.get_number_of_strings();
    
    read_info = vector<string>(bwt.size());
    

    for(uint64_t i = 0;i < N;++i){//for each read (if RC=true, only first half of reads)
        int x=0;
        {
			uint64_t j = i + x*(N/2);
			uint64_t off=0;//offset from the end of the read
			
			while(bwt[j] != TERM){
			  read_info[j] = read_ids[i].substr(0,max_id_len);
			  read_info[j].append(string("\t"));
			  read_info[j].append(to_string(off));
			  j = bwt.LF(j);
			  off++;
			}
		}
    }
    
    cout << "ID\tposition\toriginal\tmodified\tmodified.quality\tLCP>=K\tminimum?" << endl;
    for(uint64_t i=0;i<bwt.size();++i){
		cout << read_info[i] << "\t" << bwt[i] << "\t" << QUAL[i] << "\t" << (LCP_threshold[i]?"+\t":"\t") << (LCP_minima[i]?"*":"") << endl;
    }
      
  }
    
}

bool file_exists(string fileName){

  std::ifstream infile(fileName);

return infile.good();
}

/*
 * END PROCEDURES TO DEBUG
 */

int main(int argc, char** argv){
    
  srand(time(NULL));

  
  if(argc < 3) help();
  
  int opt;
  while ((opt = getopt(argc, argv, "he:q:o:k:m:v:f:t:s:DVH:")) != -1){
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
      case 'k':
        K = atoi(optarg);
        break;
      case 'm':
        m = atoi(optarg);
        break;
      case 'v':
        default_value = atoi(optarg);
        break;
	  case 'f':
        freq_threshold = atoi(optarg);
        break;
	  case 't':
        quality_threshold = atoi(optarg);
        break;
      case 's':
        TERM = atoi(optarg);
        break;
      case 'D':
        debug=true;
        break;
	  case 'V':
        verbose=true;
        break;
      case 'H':
	    input_titles = string(optarg);
        ignore_headers=false;
        break;
      default:
        help();
        return -1;
    }
  }

  K = K == -1 ? K_def : K;
  m = m == 0 ? m_def : m;
  default_value = default_value == '\0' ? default_value_def : default_value;
  quality_threshold = quality_threshold == -1 ? quality_threshold_def : quality_threshold;
  freq_threshold = freq_threshold == 0.0 ? freq_threshold_def : freq_threshold;

  if( input_dna.compare("")==0 or
      input_qual.compare("")==0 or
      output.compare("")==0
    ) help();

  if(not file_exists(input_dna)){
    cout << "Error: could not find file " << input_dna << "." << endl << endl;
    help();
  }

  if(not file_exists(input_qual)){
    cout << "Error: could not find file " << input_qual << endl << endl;
    help();
  }

  if(not ignore_headers and not file_exists(input_titles)){
    cout << "Error: could not find file " << input_titles << "." << endl << endl;
    help();
  }

  cout << "Running bfq_int..." << endl;
  cout << "\tMode: " << M;
  #if M==2
    cout << "\treplacing QS with symbol: " <<  default_value ;
  #endif
  cout << "\n\tIllumina 8-level binning: " << B << endl;
  cout << "\tK: " << K << endl;
  cout << "\tm: " << m << endl;
  cout << "\tFrequency threshold: " << freq_threshold << "\%" << endl << endl;
  cout << "Output fastq file: " << output << endl;
  cout << endl;

  cout << "Phase 1/5: loading and indexing eBWT ... " << flush;
  bwt = dna_bwt_n_t(input_dna,TERM);

  uint64_t N = bwt.get_number_of_strings();
  cout << "Number of reads: " << N << endl;
  cout << "done." << endl;
  /**/
    fprintf(stdout, "##\nmalloc_count ### current peak: %'zu\n##\n", malloc_count_peak_curr());
  /**/
  
  //detects clusters through local LCP minima
  //Phase 2/4: navigating suffix tree leaves
  //Phase 3/4: computing LCP minima
  detect_minima();
  cout << "done." << endl;
    /**/
    fprintf(stdout, "##\nmalloc_count ### current peak: %'zu\n##\n", malloc_count_peak_curr());
    /**/
	
  //start procedure run		
  cout << "\nPhase 4/5: process clusters ... " << endl;
  run();
  cout << "done." << endl;
  /**/
  fprintf(stdout, "##\nmalloc_count ### current peak: %'zu\n##\n", malloc_count_peak_curr());
  /**/
  
  //invert BWT
  cout << "\nPhase 5/5: inverting eBWT ... " << endl;
  invert();
  cout << "done." << endl;
   /**/
  fprintf(stdout, "##\nmalloc_count ### current peak: %'zu\n##\n", malloc_count_peak_curr());
  /**/
  
  cout << "**** Cluster statistics ****" << endl;
  cout << "Tot: " << num_clust << endl;
  cout << bases_inside << " (" << (double(100*bases_inside)/(bwt.size()-N)) <<  "%) bases fall inside clusters" << endl;
  cout << "Discarded: " << num_clust_discarded << "(" << (double(100*num_clust_discarded)/(num_clust)) <<  "%)" << endl;
  cout << "Ambiguous discarded: " << num_clust_amb_discarded << "(" << (double(100*num_clust_amb_discarded)/(num_clust)) <<  "%)" << endl;
  cout << "Processed: " << num_clust_mod << "(" << (double(100*num_clust_mod)/(num_clust)) <<  "%)" << endl;
  cout << "Clusters with only one symbol: " << num_clust_alleq << "(" << (double(100*num_clust_alleq)/(num_clust)) <<  "%)" << endl << endl;
  
  cout << "**** Quality statistics ****" << endl;
  cout << qs_smoothed << "/" << bwt.size() - N << " qualities have been modified (" << 100*double(qs_smoothed)/(bwt.size()-N) << "% of all qs and " <<
		100*double(qs_smoothed)/bases_inside << "% of qs inside clusters)." << endl;
  cout << "***********************" << endl << endl;
  
  cout << "**** Bases statistics ****" << endl;
  cout << modified << "/" << bwt.size() - N << " bases have been modified (" << 100*double(modified)/(bwt.size()-N) << "% of all bases and " <<
    100*double(modified)/bases_inside << "% of bases inside clusters)." << endl;
  cout << "***********************" << endl << endl;

  #if DEBUG
     cout << "Distribution of base qualities before: " << endl;
     //uint64_t sum_tot = 0;
     //for(auto x : statistics_qual_before) sum_tot+=x;

	uint64_t sum = 0;

    for(int i=32; i<127;++i){
		if(statistics_qual_before[i]>0){
			sum += statistics_qual_before[i];
			cout << (char)i << "\t" << statistics_qual_before[i] << "\t" << 100*(double(statistics_qual_before[i]))/(bwt.size() - N) << endl;
		}
    }
    cout << endl;

	cout << "sum: " << sum << endl;
	sum = 0;
	 
	cout << "Distribution of base qualities after: " << endl;

	for(int i=32;i<127;++i){
		if(statistics_qual_after[i]>0){
			sum += statistics_qual_after[i];
			cout << (char)i << "\t" << statistics_qual_after[i] << "\t" << 100*(double(statistics_qual_after[i]))/(bwt.size() - N) << endl;
		}
	}

    cout << endl;
		
	cout << "sum: " << sum << endl;
	
  #endif


  if(debug)
    print_info();

  rankbv_free(rbv);

return 0;
}
