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

int main(){
  clock_t begin = clock();
  const string filename = "test.fq";

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

  ////////////////////////////////////////////////////////////////////////////
  // PASS THROUGH FILE 
  ////////////////////////////////////////////////////////////////////////////

  short line = 0,base = 0, base_ind = 0;
  size_t nreads = 1;
  char c, buff[num_bases];

  // to log progress
  size_t max_block = 20, cur_block = 1;

  cerr << "Started analysis of " << filename << "\n";

  // read character by character
  for(curr = first; curr < last - 1; ++curr) {
    // Log progress
    if(curr - first >= (last - first)*cur_block/max_block){
      cerr << "Approximately " << 100*cur_block/max_block <<
              "% complete for " << filename << "\n";
      cur_block++;
    }

    // Process line
    if (*curr !=  '\n'){
     // base read
      if(line == 1){
        base_ind = (*(curr) >> 1) & 7;
        base_count[base_ind][base]++;
        buff[base++] = *curr;
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
      nreads++;
    }
  }

  munmap(mmap_data, st.st_size);

  cerr << "Finished reading " << nreads << " reads\n";
  cerr << "Summarizing totals...\n";

  nreads /= 4;
  for(size_t j = 0; j < num_chars; ++j){
    for(size_t i = 0; i < num_bases; ++i){
      total_char[j] += base_count[j][i];
      if(i == 1)
        total_quality[j] += base_quality[j][i];
      total_bases += base_count[j][i];
    }
  }

  cerr << "Average A quality on read 2: ";
  for(size_t i = 0; i < 60; ++i)  
    cerr << base_quality[0][i] / static_cast<double>(base_count[0][i]) - ascii_to_quality << " ";
  cerr << "\n";

  cerr << "Average gc % " << 100.0 * (total_char[1] + total_char[3]) / static_cast<double>(total_bases) << "\n";
  cerr << "Average n % " << 100.0 * total_char[7] / static_cast<double>(total_bases) << "\n";
  cerr << "Average length " << total_bases / static_cast<double>(nreads) << "\n";

  cerr << "Elapsed time: " << (clock() - begin) / CLOCKS_PER_SEC << " seconds\n";
}


