/* fqc: what this program does
*
* Copyright (C) 2019 Guilherme De Sena Brandine
*
* Authors: Guilherme De Sena Brandine
*
* This program is free software: you can redistribute it and/or
* modify it under the terms of the GNU General Public License as
* published by the Free Software Foundation, either version 3 of the
* License, or (at your option) any later version.
*
* This program is distributed in the hope that it will be useful, but
* WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
* General Public License for more details.
*/
#include <stdio.h>
#include <unistd.h>
#include <fcntl.h>

#include <sys/stat.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <vector>
#include <ctime>

#include <string>
#include <algorithm>
#include <iostream>

using std::string;
using std::runtime_error;
using std::cerr;
using std::endl;
using std::vector;
using std::reverse;
static string
int_to_seq(size_t v, const size_t &kmer_size) {
  string ans;
  for(size_t i = 0; i < kmer_size; ++i) {
    switch (v & 7) {
      case 0: ans.push_back('A'); break;
      case 1: ans.push_back('C'); break;
      case 2: ans.push_back('T'); break;
      case 3: ans.push_back('G'); break;
      case 7: ans.push_back('N'); break;
    }
    v >>= 3;
  }
  if(v > 0)
    throw std::runtime_error("bad kmer");

  reverse(ans.begin(),ans.end());
  return ans;
}

int main(int argc, char **argv){
  clock_t begin = clock();

  // Hardcoded until I add option parser
  const string filename = string(argv[1]);

  ////////////////////////////////////////////////////////////////////////////
  // MEMORY MAP
  ////////////////////////////////////////////////////////////////////////////

  // open the file
  int fd = open(filename.c_str(), O_RDONLY, 0);
  if (fd == -1)
    throw runtime_error("failed to open fastq file: " + filename);

  // get the file size
  struct stat st;
  fstat(fd, &st);

  // execute mmap
  void* mmap_data = mmap(NULL, st.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
  close(fd);
  if (mmap_data == MAP_FAILED)
    throw runtime_error("failed to mmap fastq file: " + filename);

  cerr << "Fastq size: " << st.st_size /(1024*1024)<< "Mb" << endl;
  char *first = static_cast<char*>(mmap_data);
  char *last = first + st.st_size - 1;

  ////////////////////////////////////////////////////////////////////////////
  // STATS
  ////////////////////////////////////////////////////////////////////////////

  // This can be a very large number, the maximum illumina read size. Does not
  // affect memory very much
  const size_t num_bases = 1000;
  const double ascii_to_quality = 33.0;
  const size_t num_chars = 8; // A = 000, C = 001, T = 010, G = 011, N = 111

  /***********ALL****************/
  vector<size_t> total_char (num_chars, 0);
  vector<size_t> total_quality (num_chars, 0);
  size_t total_bases = 0;

  /***********BASE****************/

  // counts the number of bases in every read position
  vector<vector<size_t>> base_count(num_chars, vector<size_t>(num_bases,0));

  // Sum of base qualities in every read position
  vector<vector<size_t>> base_quality(num_chars, vector<size_t>(num_bases,0));

  /***********READ***************/

  // Distribution of read lengths
  vector<size_t> read_length_freq(num_bases,0);

  /**********KMER****************/

  const size_t kmer_size = 8;
  const size_t kmer_mask = (1 << (3*kmer_size)) - 1;

  // A 3^(K+1) vector to count all possible kmers
  vector<size_t> kmer_count(kmer_mask + 1, 0);

  ////////////////////////////////////////////////////////////////////////////
  // PASS THROUGH FILE
  ////////////////////////////////////////////////////////////////////////////
  size_t base_ind = 0;  // 0,1,2,3 or 7
  size_t read_pos = 0; // which base we are at in the read
  size_t nreads = 0; // nlines/4
  char c; // to avoid accessing *curr multiple times
  char buff[num_bases]; // I'll store the sequence to know the quality afterwards
  size_t cur_kmer = 0; //kmer hash

  cerr << "Started analysis of " << filename << "\n";

  // read character by character
  for(char *curr = first; curr < last;) {
    ++nreads;

    // fast forward first line
    for(; *curr != '\n'; ++curr);
    ++curr; //skips \n from line 1

    //*************READ LINE*******************/
    read_pos = 0;
    for(; *curr != '\n'; ++curr){
      c = *curr;

      // Transforms base into 3-bit index
      // Bits 2,3 and 4 of charcters A,C,G,T and N are distinct so we can just
      // use them to get the index instead of doing if-elses.
      base_ind = (c >> 1) & 7;

      // Increments count of base
      base_count[base_ind][read_pos]++;

      // Need this to know what base the quality is associated to, therefore
      // we will store each read in memory
      buff[read_pos] = c;
      cur_kmer = ((cur_kmer << 3) | base_ind);

      // If we already read >= k bases, increment the k-mer count
      if(read_pos >= kmer_size - 1)
        kmer_count[cur_kmer & kmer_mask]++;
      
      read_pos++;

    }
    ++curr; // skips \n from line 2

    // count read length frequency
    read_length_freq[read_pos]++;

    // fast forward third line 
    for(; *curr != '\n'; ++curr);
    ++curr; //skips \n from line 3

    //*************QUALITY LINE*******************/
    read_pos = 0;
    for(; (*curr != '\n') && (curr < last); ++curr){
      // Gets the base from buffer, finds index, adds quality value
      base_quality[(buff[read_pos] >> 1) & 7][read_pos] += *curr;
      read_pos++;
    }
    ++curr; // skip \n or increments from last which is not a problem

  }

  // Deallocates memory
  munmap(mmap_data, st.st_size);

  cerr << "Finished reading " << nreads << " reads\n";
  cerr << "Summarizing totals...\n";

  for(size_t j = 0; j < num_chars; ++j){
    for(size_t i = 0; i < 60; ++i){
      total_char[j] += base_count[j][i];
      if(i == 1)
        total_quality[j] += base_quality[j][i];
      total_bases += base_count[j][i];
    }
  }

  cerr << "Average A quality: ";
  for(size_t i = 0; i < 60; ++i)
    cerr << base_quality[0][i] / static_cast<double>(base_count[0][i]) - ascii_to_quality << " ";
  cerr << "\n\n";

  cerr << "Average A frequency: ";
  for(size_t i = 0; i < 60; ++i) {
    const double denom = (base_count[0][i] + base_count[1][i] +
                          base_count[2][i] + base_count[3][i] +
                          base_count[7][i]);
    cerr << base_count[0][i]/denom << " ";
  }
  cerr << "\n";

  cerr << "Number of reads: " << nreads << "\n";
  cerr << "gc % " << 100.0 * (total_char[1] + total_char[3]) / static_cast<double>(total_bases) << "\n";
  cerr << "n % " << 100.0 * total_char[7] / static_cast<double>(total_bases) << "\n";
  cerr << "Average length " << total_bases / static_cast<double>(nreads) << "\n";

  for(size_t i = 0; i < kmer_mask; ++i)
    if(kmer_count[i] > 100000000)
      cerr << int_to_seq(i, kmer_size) << "\t" << kmer_count[i] << "\n";

  cerr << "Elapsed time: " << (clock() - begin) / CLOCKS_PER_SEC << " seconds\n";
}
