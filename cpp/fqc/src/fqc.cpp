#include <stdio.h>
#include <unistd.h>
#include <fcntl.h>

#include <sys/stat.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <vector>
#include <ctime>

#include <string>
#include <iostream>
using std::string;
using std::runtime_error;
using std::cerr;
using std::vector;

string int_to_seq(size_t v){
  string ans;
  while(v > 0){
    switch (v & 7){
      case 0: ans.push_back('A'); break;
      case 1: ans.push_back('C'); break;
      case 2: ans.push_back('T'); break;
      case 3: ans.push_back('G'); break;
      case 7: ans.push_back('N'); break;
    }
    v >>= 3;
  }
  return ans;
}

int main(){
  clock_t begin = clock();
  const string filename = "SRR1853178_pass_2.fastq";

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

  cerr << "Fastq size: " << st.st_size /(1024*1024)<< "Mb \n";
  char *curr = static_cast<char*>(mmap_data);
  char *first = curr;
  char *last = curr + st.st_size;

  ////////////////////////////////////////////////////////////////////////////
  // STATS
  ////////////////////////////////////////////////////////////////////////////
  const size_t num_bases = 1000;
  const double ascii_to_quality = 33.0;
  short num_chars = 8; // A = 000, C = 001, T = 010, G = 011, N = 111

  /***********ALL****************/
  vector<size_t> total_char (num_chars, 0);
  vector<size_t> total_quality (num_chars, 0);
  size_t total_bases = 0;

  /***********BASE****************/
  vector<vector<size_t>> base_count(num_chars, vector<size_t>(num_bases,0)); // avg per base gc
  vector<vector<size_t>> base_quality(num_chars, vector<size_t>(num_bases,0)); // avg per base gc

  /***********READ***************/
  vector<size_t> read_length_freq(num_bases,0); //read length frequency

  /**********KMER****************/
  const size_t kmer_size = 8;
  const size_t kmer_lookup_len = 1<< (3*(kmer_size + 1));
  const size_t kmer_mask = 1 << (3*kmer_size) - 1;
  vector<size_t> kmer_count(kmer_lookup_len, 0);

  ////////////////////////////////////////////////////////////////////////////
  // PASS THROUGH FILE 
  ////////////////////////////////////////////////////////////////////////////
  short line = 0,kmer_base = 0, base = 0, base_ind = 0;
  size_t nreads = 1;
  char c; 
  char buff[num_bases]; // I'll store the sequence to know the quality afterwards

  // Hash of kmer 
  size_t cur_kmer = 0;

  cerr << "Started analysis of " << filename << "\n";

  // read character by character
  for(curr = first; curr < last - 1; ++curr) {
    // Process line
    if (*curr !=  '\n'){
     // base read
      if(line == 1){
        c = *curr;

        // Transforms base into 3-bit index
        base_ind = (c >> 1) & 7;

        // Increments count of base
        base_count[base_ind][base]++;

        // Need this to know what base the quality is associated to
        buff[base] = c;
        cur_kmer = ((cur_kmer << 3) | base_ind);

        base++;
        if(kmer_base == kmer_size){
          kmer_count[cur_kmer]++;
          cur_kmer &= kmer_mask;
        } else {
          kmer_base++;
        }
      }

      // quality read
      else if(line == 3){
        base_ind = (buff[base] >> 1) & 7;
        base_quality[base_ind][base++] += *curr;
      }
    
    // Start new line 
    } else {
      // count read length frequency
      if(line == 1)
        read_length_freq[base]++;

      line = (line+1) & 3;
      base = 0;
      cur_kmer = 0;
      kmer_base = 0;
      nreads++;
    }
  }

  // Deallocates memory
  munmap(mmap_data, st.st_size);

  nreads /= 4;
  cerr << "Finished reading " << nreads << " reads\n";
  cerr << "Summarizing totals...\n";

  for(size_t j = 0; j < num_chars; ++j){
    for(size_t i = 0; i < num_bases; ++i){
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
  for(size_t i = 0; i < 60; ++i)  
    cerr << base_count[0][i] / static_cast<double>(base_count[0][i] + base_count[1][i] + base_count[2][i] + base_count[3][i] + base_count[7][i]) << " ";
  cerr << "\n";


  cerr << "Average gc % " << 100.0 * (total_char[1] + total_char[3]) / static_cast<double>(total_bases) << "\n";
  cerr << "Average n % " << 100.0 * total_char[7] / static_cast<double>(total_bases) << "\n";
  cerr << "Average length " << total_bases / static_cast<double>(nreads) << "\n";
  cerr << "Elapsed time: " << (clock() - begin) / CLOCKS_PER_SEC << " seconds\n";
}


