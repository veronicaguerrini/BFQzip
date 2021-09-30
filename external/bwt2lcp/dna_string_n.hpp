// Copyright (c) 2018, Nicola Prezza.  All rights reserved.
// Use of this source code is governed
// by a MIT license that can be found in the LICENSE file.

/*
 * dna_string_n.hpp
 *
 *  Created on: Dec 6, 2018
 *      Author: nico
 *
 *  Optimized string with rank on DNA alphabet: {A,C,G,N,T,TERM}, where the terminator TERM is defined in include.hpp
 *
 *  One access or a parallel rank for all 5 letters A,C,G,T,N causes only 1 cache miss in the worst case
 *
 *  Max string length: 2^64
 *
 *  Supports very efficient (1 cache miss) parallel rank for (A,C,G,N,T), and (1 cache miss) single rank for TERM
 *
 *  Data is stored and cache-aligned in blocks of 512 bits (64 bytes)
 *
 *  Size of the string: 512/117 < 4.38n bits, where n = string length
 *
 *  512-bits Block layout:
 *
 *  | 32-bit rank A | 32-bit rank C | 32-bit rank G | 32-bit rank T|
 *  | 117-bit 1st bits | 11-bit partial rank N | 117-bit 2nd bits | 11-bit partial rank N | 117-bit 3rd bits | 11-bit partial rank N|
 *
 *
 */

#ifndef INTERNAL_DNA_STRING_N_HPP_
#define INTERNAL_DNA_STRING_N_HPP_

#define SUPERBLOCK_SIZE_N_N 3925868544 		//number of characters in a superblock = 117*2^25 characters
#define BLOCKS_PER_SUPERBLOCK_N 33554432	//blocks in a superblock = 2^25
#define BYTES_PER_SUPERBLOCK_N 2147483648	//bytes in a superblock = 64*BLOCKS_PER_SUPERBLOCK_N
#define BLOCK_SIZE_N 117 					//number of characters inside a block
#define BYTES_PER_BLOCK_N 64				//bytes in a block of 512 bits
#define ALN_N 64							//alignment

#include "include.hpp"

class dna_string_n{

public:

	dna_string_n(){}

	/*
	 * constructor from ASCII file
	 */
	dna_string_n(string path, char TERM = '#'){

		this->TERM = TERM;

		n = uint64_t(filesize(path));

		n_superblocks = (n+1)/SUPERBLOCK_SIZE_N_N + ((n+1)%SUPERBLOCK_SIZE_N_N != 0);
		n_blocks = (n+1)/BLOCK_SIZE_N + ((n+1)%BLOCK_SIZE_N != 0);
		nbytes = (n_blocks * BYTES_PER_BLOCK_N);//number of bytes effectively filled with data

		superblock_ranks = vector<p_rank_n>(n_superblocks);

		/*
		 * this block of code ensures that data is aligned by 64 bytes = 512 bits
		 */
		memory = vector<uint8_t>(nbytes+ALN_N,0);
		data = memory.data();
		while(uint64_t(data) % ALN_N != 0) data++;

		//cout << "alignment of data: " << (void*)data << endl;

		{

			ifstream ifs(path);

			string BUF(BLOCK_SIZE_N,'A');

			for(uint64_t i = 0; i<n; ++i){

				char c;

				ifs.read((char*)&c, sizeof(char));

				BUF[i%BLOCK_SIZE_N] = c;

				if(c!='A' and c!='C' and c!='G' and c!='N' and c!='T' and c!=TERM){
					cout << "Error while reading file: read forbidden character '" <<  c << "' (ASCII code " << int(c) << ")." <<
					"Only A,C,G,N,T, and " << TERM << " are admitted in the input BWT!" << endl <<
					"Possible solution: if the unknown character is the terminator, you can solve the problem by adding option \"-t " << int(c) << "\"." << endl;

					exit(1);
				}

				//buffer is full
				if((i%BLOCK_SIZE_N) == BLOCK_SIZE_N-1) set(i/BLOCK_SIZE_N, BUF);

			}

			if(n % BLOCK_SIZE_N != 0) set(n/BLOCK_SIZE_N, BUF);

		}

		build_rank_support();

		assert(check_content(path));
		assert(check_rank());

	}

	//return i-th character
	char operator[](uint64_t i){

		/*
		 * internal encoding (does not reflect lexicographic ordering, which is the standard alphabetical one)
		 *
		 * A     000
		 * C     001
		 * G     010
		 * T     011
		 * TERM  100
		 * N     101
		 *
		 */

		assert(i<n);

		uint64_t superblock_number = i / SUPERBLOCK_SIZE_N_N;
		uint64_t superblock_off = i % SUPERBLOCK_SIZE_N_N;
		uint64_t block_number = superblock_off / BLOCK_SIZE_N;
		uint64_t block_off = superblock_off % BLOCK_SIZE_N;

		//chars[2,1,0] contains 1st, 2nd, 3rd most significant bits of the 117 characters (in the 117-bits prefix of each block of 128 bits)
		__uint128_t* chars = (__uint128_t*)(data + superblock_number*BYTES_PER_SUPERBLOCK_N + block_number*BYTES_PER_BLOCK_N);

		uint64_t b =	((chars[0]>>(128-(block_off+1)))&0x1) +
						(((chars[1]>>(128-(block_off+1)))&0x1)<<1) +
						(((chars[2]>>(128-(block_off+1)))&0x1)<<2);

		return 	(b == 0)*'A' +
				(b == 1)*'C' +
				(b == 2)*'G' +
				(b == 3)*'T' +
				(b == 4)*TERM +
				(b == 5)*'N';

	}

	/*
	 * Parallel rank of (A,C,G,N,T) at position i.
	 */
	p_rank_n parallel_rank(uint64_t i){

		uint64_t superblock_number = i / SUPERBLOCK_SIZE_N_N;
		uint64_t superblock_off = i % SUPERBLOCK_SIZE_N_N;
		uint64_t block_number = superblock_off / BLOCK_SIZE_N;
		uint64_t block_off = superblock_off % BLOCK_SIZE_N;

		p_rank_n superblock_r = superblock_ranks[superblock_number];
		p_rank_n block_r = get_counters(superblock_number,block_number);

		return superblock_r + block_r + block_rank(superblock_number, block_number, block_off);

	}

	/*
	 * standard rank. c can be A,C,G,T, or TERM
	 */
	uint64_t rank(uint64_t i, uint8_t c){

		p_rank_n pr = parallel_rank(i);

		if(c==TERM) return rank_non_dna(i);

		switch(c){
			case 'A' : return pr.A; break;
			case 'C' : return pr.C; break;
			case 'G' : return pr.G; break;
			case 'N' : return pr.N; break;
			case 'T' : return pr.T; break;
		}

		return 0;

	}

	/*
	 * return number of non-dna symbols in the prefix of length i of the text. At most 1 cache miss!
	 */
	uint64_t rank_non_dna(uint64_t i){

		assert(i<=n);
		auto r = parallel_rank(i);

		assert(r.A + r.C + r.G + r.T + r.N <= i);

		return i - (r.A + r.C + r.G + r.T + r.N);

	}

	uint64_t serialize(std::ostream& out){

		uint64_t w_bytes = 0;

		out.write((char*)&n,sizeof(n));
		out.write((char*)&nbytes,sizeof(nbytes));
		out.write((char*)&n_superblocks,sizeof(n_superblocks));
		out.write((char*)&n_blocks,sizeof(n_blocks));

		w_bytes += sizeof(n) + sizeof(nbytes) + sizeof(n_superblocks) + sizeof(n_blocks);

		out.write((char*)superblock_ranks.data(),n_superblocks*sizeof(p_rank_n));
		w_bytes += n_superblocks*sizeof(p_rank_n);

		out.write((char*)data,nbytes*sizeof(uint8_t));
		w_bytes += nbytes*sizeof(uint8_t);

		return w_bytes;

	}

	void load(std::istream& in) {

		in.read((char*)&n,sizeof(n));
		in.read((char*)&nbytes,sizeof(nbytes));
		in.read((char*)&n_superblocks,sizeof(n_superblocks));
		in.read((char*)&n_blocks,sizeof(n_blocks));

		superblock_ranks = vector<p_rank_n>(n_superblocks);
		in.read((char*)superblock_ranks.data(),n_superblocks*sizeof(p_rank_n));

		memory = vector<uint8_t>(nbytes+ALN_N,0);
		data = memory.data();
		while(uint64_t(data) % ALN_N != 0) data++;
		in.read((char*)data,nbytes*sizeof(uint8_t));

		assert(check_rank());

	}

	uint64_t size(){
		return n;
	}

private:

	void build_rank_support(){

		p_rank_n superblock_r = {};
		p_rank_n block_r = {};

		for(uint64_t bl = 0; bl < n_blocks-1; ++bl){

			uint64_t superblock_number = bl/BLOCKS_PER_SUPERBLOCK_N;
			uint64_t block_number = bl%BLOCKS_PER_SUPERBLOCK_N;

			if(block_number == 0){

				superblock_ranks[superblock_number]=superblock_r;
				block_r = {};

			}

			set_counters(superblock_number, block_number,block_r);

			p_rank_n local_rank = block_rank(superblock_number, block_number);

			block_r = block_r + local_rank;
			superblock_r = superblock_r + local_rank;

		}

		uint64_t superblock_number = (n_blocks-1)/BLOCKS_PER_SUPERBLOCK_N;
		uint64_t block_number = (n_blocks-1)%BLOCKS_PER_SUPERBLOCK_N;

		if(block_number == 0){

			superblock_ranks[superblock_number]=superblock_r;
			block_r = {};

		}

		set_counters(superblock_number, block_number,block_r);

	}


	/*
	 * set i-th block to s. Assumption: s.length() == BLOCK_SIZE_N
	 */
	void set(uint64_t i, string & s){

		assert(s.length()==BLOCK_SIZE_N);
		assert(i<n_blocks);

		uint64_t superblock_number = i / BLOCKS_PER_SUPERBLOCK_N;
		uint64_t block_number = i % BLOCKS_PER_SUPERBLOCK_N;

		//a block contains 117 characters
		//chars[2,1,0] contains 1st, 2nd, 3rd most significant bits of the 117 characters
		__uint128_t* chars = (__uint128_t*)(data + superblock_number*BYTES_PER_SUPERBLOCK_N + block_number*BYTES_PER_BLOCK_N);

		chars[2] = 0;
		chars[1] = 0;
		chars[0] = 0;

		/*
		 * internal encoding (does not reflect lexicographic ordering, which is the standard alphabetical one)
		 *
		 * A     000
		 * C     001
		 * G     010
		 * T     011
		 * TERM  100
		 * N     101
		 *
		 */

		for(auto c : s){

			if(c==TERM){

				chars[2] = (chars[2] << 1) | __uint128_t(0x1);
				chars[1] = chars[1] << 1;
				chars[0] = chars[0] << 1;

			}else{

				switch(c){

				case 'A' : 	chars[2] = chars[2] << 1;
							chars[1] = chars[1] << 1;
							chars[0] = chars[0] << 1;break;

				case 'C' : 	chars[2] = chars[2] << 1;
							chars[1] = chars[1] << 1;
							chars[0] = (chars[0] << 1) | __uint128_t(0x1);break;

				case 'G' : 	chars[2] = chars[2] << 1;
							chars[1] = (chars[1] << 1) | __uint128_t(0x1);
							chars[0] = chars[0] << 1;break;

				case 'T' : 	chars[2] = chars[2] << 1;
							chars[1] = (chars[1] << 1) | __uint128_t(0x1);
							chars[0] = (chars[0] << 1) | __uint128_t(0x1);break;

				case 'N' : 	chars[2] = (chars[2] << 1) | __uint128_t(0x1);
							chars[1] = chars[1] << 1;
							chars[0] = (chars[0] << 1) | __uint128_t(0x1);break;
				}

			}

		}

		//leave 11 (=128-117) bits free at the end of each 128-bit block. We will store N's partial rank here.

		chars[2] = chars[2] << 11;
		chars[1] = chars[1] << 11;
		chars[0] = chars[0] << 11;

	}

	/*
	 * rank in block given as coordinates: superblock, block, offset in block
	 */
	inline p_rank_n block_rank(uint64_t superblock_number, uint64_t block_number, uint64_t block_off=BLOCK_SIZE_N){

		/*
		 * internal encoding (does not reflect lexicographic ordering, which is the standard alphabetical one)
		 *
		 * A     000
		 * C     001
		 * G     010
		 * T     011
		 * TERM  100
		 * N     101
		 *
		 */

		assert(block_off<=BLOCK_SIZE_N);

		//starting address of the block
		uint8_t* start = data + superblock_number*BYTES_PER_SUPERBLOCK_N + block_number*BYTES_PER_BLOCK_N;

		//chars[0..3] contains 1st, 2nd, 3rd bits of the 128 characters
		__uint128_t* chars = (__uint128_t*)(start);

		__uint128_t PAD = ((~__uint128_t(0))>>block_off);

		//no character's code begins with 11, so we pad most 2 significant bits
		__uint128_t b2 = chars[2] | PAD;
		__uint128_t b1 = chars[1] | PAD;
		__uint128_t b0 = chars[0];

		return {

			popcount128((~b2) & (~b1) & (~b0)), // A = 000
			popcount128((~b2) & (~b1) & (b0)), // C = 001
			popcount128((~b2) & (b1) & (~b0)), // G = 010
			popcount128((b2) & (~b1) & (b0)), // N = 101
			popcount128((~b2) & (b1) & (b0)), // T = 011

		};

	}

	bool check_rank(){

		p_rank_n p = {};

		bool res = true;

		for(uint64_t i=0;i<size();++i){

			auto r = parallel_rank(i);

			if(p != r){

				res = false;
				/*cout << "Error in local rank at position " << n << " "
				<< p.A << "/" << r.A << " "
				<< p.C << "/" << r.C << " "
				<< p.G << "/" << r.G << " "
				<< p.T << "/" << r.T << endl;*/

			}

			p.A += (operator[](i)=='A');
			p.C += (operator[](i)=='C');
			p.G += (operator[](i)=='G');
			p.T += (operator[](i)=='T');
			p.N += (operator[](i)=='N');

		}

		auto r = parallel_rank(n);

		if(p != r){

			res = false;
			/*cout << "Error in local rank at position " << n << " "
			<< p.A << "/" << r.A << " "
			<< p.C << "/" << r.C << " "
			<< p.G << "/" << r.G << " "
			<< p.T << "/" << r.T << endl;*/

		}

		if(res){

			cout << "rank is correct" << endl;

		}else{

			cout << "rank is not correct" << endl;

		}

		return res;

	}

	/*
	 * check that the string contains exactly the same characters as the file in path
	 */
	bool check_content(string path){

		ifstream ifs(path);

		bool res = true;

		for(uint64_t i=0;i<n;++i){

			char c;
			ifs.read((char*)&c, sizeof(char));

			if(operator[](i) != c) res = false;

		}

		if(res){

			cout << "string content is valid" << endl;

		}else{

			cout << "string content is not valid" << endl;

		}

		return res;

	}

	/*
	 * set counters in the i-th block to r
	 */
	void set_counters(uint64_t superblock_number, uint64_t block_number, p_rank_n r){

		//block start
		uint8_t* start = data + superblock_number*BYTES_PER_SUPERBLOCK_N + block_number*BYTES_PER_BLOCK_N;
		uint32_t * block_ranks = (uint32_t*)(start+48);

		block_ranks[0] = r.A;
		block_ranks[1] = r.C;
		block_ranks[2] = r.G;
		block_ranks[3] = r.T;

		uint64_t * chars = (uint64_t*)(start);

		chars[0] += (r.N & MASK);
		chars[2] += ((r.N>>11) & MASK);
		chars[4] += ((r.N>>22) & MASK);

		//the 32 bits of N counter are stored in the length-11 suffixes of chars[0,2,4]
		uint64_t rank_N = (chars[0]&MASK) + ((chars[2]&MASK)<<11) + ((chars[4]&MASK)<<22);

		assert(get_counters(superblock_number,block_number) == r);

	}

	/*
	 * get counters of the i-th block
	 */
	inline p_rank_n get_counters(uint64_t superblock_number, uint64_t superblock_off){

		//block start
		uint8_t* start = data + superblock_number*BYTES_PER_SUPERBLOCK_N + superblock_off*BYTES_PER_BLOCK_N;
		uint32_t * block_ranks = (uint32_t*)(start+48);

		uint64_t * chars = (uint64_t*)(start);

		//the 32 bits of N counter are stored in the length-11 suffixes of chars[0,2,4]
		uint64_t rank_N = (chars[0]&MASK) + ((chars[2]&MASK)<<11) + ((chars[4]&MASK)<<22);

		return {
			block_ranks[0],
			block_ranks[1],
			block_ranks[2],
			rank_N,
			block_ranks[3]
		};

	}

	char TERM = '#';

	static const uint64_t MASK = (uint64_t(1)<<11)-1;

	uint64_t n_superblocks = 0;
	uint64_t n_blocks = 0;

	vector<uint8_t> memory; //allocated memory

	//data aligned with blocks of 64 bytes = 512 bits
	uint8_t * data = NULL;

	vector<p_rank_n> superblock_ranks;

	uint64_t nbytes = 0; //bytes used in data
	uint64_t n = 0;

};


#endif /* INTERNAL_DNA_STRING_N_HPP_ */
