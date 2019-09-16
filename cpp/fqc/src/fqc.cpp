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

#include <iostream>
#include <iomanip>
#include <fstream>

#include <vector>
#include <array>
#include <ctime>
#include <string>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "OptionParser.hpp"
#include "zlib_wrapper.hpp"

using std::string;
using std::runtime_error;
using std::cerr;
using std::endl;
using std::vector;
using std::array;
using std::reverse;
using std::ostream;
using std::ofstream;
using std::cout;
using std::unordered_map;
using std::unordered_set;
using std::pair;
using std::make_pair;
using std::ifstream;
using std::sort;
using std::max;

/*************************************************************
 ******************** AUX FUNCTIONS **************************
 *************************************************************/

// converts 64 bit integer to a sequence string by reading 2 bits at a time and
// converting back to ACTG
static inline string
int_to_seq(size_t v, const size_t &kmer_size) {
  string ans;
  for (size_t i = 0; i < kmer_size; ++i) {
    switch (v & 3) {
      case 0: ans.push_back('A'); break;
      case 1: ans.push_back('C'); break;
      case 2: ans.push_back('T'); break;
      case 3: ans.push_back('G'); break;
    }
    v >>= 2;
  }
  if (v > 0)
    throw std::runtime_error("bad kmer");

  reverse(ans.begin(), ans.end());
  return ans;
}

// Converts A,T,G,C to 2-bit values
static inline size_t
atgc_to_2bit(const char &c){
  return ((c >> 1) & 3);
}

// a 64-bit hash of the last 32 characters of a read
static inline size_t
get_suffix_hash(const string &s) {
  size_t ans = 0;
  for (size_t i = 0; i < s.size(); ++i)
    if (i >= s.size() - 32)
      ans = ((ans << 2) | (atgc_to_2bit(s[i])));

  return ans;
}

/*************************************************************
 ******************** FASTQ STATS ****************************
 *************************************************************/

struct FastqStats {
 public:
  /************************************************************
   *************** UNIVERSAL CONSTANTS ************************
   ************************************************************/

  // This can be a very large number, the maximum illumina read size.
  // Affects base memory footprint. Reads whose size is larger than this value
  // will have their last nucleotides ignored and not summarized. 
  static const size_t kNumBases = 1e4;

  // Value to subtract quality characters to get the actual quality value
  static const size_t kBaseQuality = 33;  // The ascii for the lowest quality

  // Smallest power of two that comprises all possible Illumina quality values.
  // Illumina gives qualities from 0 to 40, therefore we set it as 64. Power of
  // is to avoid double pointer jumps and to get indices with bit shifts. 
  static const size_t kNumQualityValues = 64; 
  static const size_t kBitShiftQuality = 6;  // log 2 of value above

  // How many possible nucleotides (must be power of 2!)
  static const size_t kNumNucleotides = 4;  // A = 00,C = 01,T = 10,G = 11
  static const size_t kBitShiftNucleotide = 2;  // log 2 of value above

  // threshold for a sequence to be considered  poor quality
  static const size_t kPoorQualityThreshold = 20;

  // Number of reads in which to do slower operations such as hashing strings
  static const size_t kNumSlowReads = 1e6;

  // fraction of the number of slow reads a sequence needs to be seen to be
  // considered a candiate for overrepresentation
  constexpr static double kMinFracFrequency = 0.001;

  // Position in which to store prefix for duplication estimate
  static const size_t kPrefixPosition = 32;

  // Kmer size given as input
  size_t kmer_size;

  // mask to get only the first 2*k bits of the sliding window
  size_t kmer_mask;

  /*********************************************************
   *********** METRICS COLLECTED DURING IO *****************
   *********************************************************/
  /*********** PER BASE METRICS ****************/

  // counts the number of bases in every read position
  array<size_t, kNumNucleotides * kNumBases> base_count;

  // counts the number of N bases in every read position
  array<size_t, kNumNucleotides * kNumBases> n_base_count;

  // Sum of base qualities in every read position
  array<size_t, kNumQualityValues * kNumBases> base_quality;

  // Sum of N nucleotide base qualities in every position
  array<size_t, kNumNucleotides * kNumBases> n_base_quality;

  /*********** PER QUALITY VALUE METRICS ****************/
  // Counts of quality in each base position
  array<size_t, kNumQualityValues * kNumBases> position_quality_count;

  // Counts of average quality (truncated) per sequence
  array<size_t, kNumQualityValues> quality_count;

  /*********** PER GC VALUE METRICS ****************/
  array<size_t, 101> gc_count;

  /********** KMER FREQUENCY ****************/
  // A 2^K + 1 vector to count all possible kmers
  vector<size_t> kmer_count;

  /*********** PER READ METRICS ***************/
  // Distribution of read lengths
  array<size_t, kNumBases> read_length_freq;

  // I need this to know what to divide each base by
  // when averaging content, bases, etc. It stores, for every element i, how
  // many reads are of length >= i, ie, how many reads actually have a
  // nucleotide at position i
  array<size_t, kNumBases> cumulative_read_length_freq;

  /*********** SLOW STUFF *******************/
  unordered_map <string, size_t> seq_count;
  unordered_map <size_t, vector <size_t> > suffix_count;
  unordered_map <size_t, vector <string> > suffix_seq;
  unordered_set <size_t> suffix_lookup;
  unordered_map <size_t, size_t> prefix_count;

  /*********************************************************
   *********** METRICS SUMMARIZED AFTER IO *****************
   *********************************************************/

  /*********** PASS = 0 WARN = 1 FAIL = 2 **************/
  string pass_basic_statistics,
         pass_per_base_sequence_quality,
         pass_per_tile_sequence_quality,
         pass_per_sequence_quality_scores,
         pass_per_base_sequence_content,
         pass_per_sequence_gc_content,
         pass_per_base_n_content,
         pass_sequence_length_distribution,
         pass_sequence_duplication_levels,
         pass_overrepresented_sequences,
         pass_duplicate_sequences,
         pass_kmer_content,
         pass_adapter_content;

  /*********** SINGLE NUMBERS FROM THE ENTIRE FASTQ ****************/
  size_t total_bases;  // sum of all bases in all reads
  size_t avg_read_length;  // average of all read lengths
  size_t avg_gc;  // sum of g bases + c bases
  size_t avg_n;  // sum of n bases
  size_t num_reads;  // total number of lines read
  size_t num_poor;  // reads whose average quality was <= poor

  // Quantiles for the position_quality_count
  array<size_t, kNumBases> ldecile, lquartile, median, uquartile, udecile;
  array<double, kNumBases> mean;

  // Percentages for per base sequence content
  array<double, kNumBases> a_pct,
                           c_pct,
                           t_pct,
                           g_pct,
                           n_pct;

  // For sequence duplication levels
  //1 to 9, >10, >50, >100, >500, >1k, >5k, >10k+
  array<double, 16> percentage_deduplicated;
  array<double, 16> percentage_total;

  // Overrepresented sequences and the number of times they were seen
  vector<pair<string, size_t>> overrep_sequences;
  /*********************************************************
   ********************* NON-WRITE FUNCTIONS ***************
   *********************************************************/

  // Default constructor that zeros everything
  explicit FastqStats(const size_t _kmer_size);

  // Summarize all statistics we need before writing
  void calc_summary();

  /******* DUPLICATION AND OVERREPRESENTATION *******/
  // Makes a hash map with keys as 32-bit suffixes and values as all the
  // candidate frequent sequences with that given suffix
  void make_suffix_lookup();
  
  /*********************************************************
   ********************* WRITING FUNCTIONS *****************
   *********************************************************/

  // HEADER
  void write_header(ostream &os);

  // BASIC STATISTICS: columns = measure, value
  void write_basic_statistics(ostream &os, string filename);

  // PER BASE SEQUENCE QUALITY columns = base pos, quality
  void write_per_base_quality(ostream &os);

  // PER TILE SEQUENCE QUALITY: tile x position matrix
  void write_per_tile_sequence_quality(ostream &os);

  // PER SEQUENCE QUALITY SCORES: columns = quality value, count
  void write_per_sequence_quality(ostream &os);

  // PER BASE SEQUENCE CONTENT: columns = Base,G,A,T,C
  void write_per_base_sequence_content(ostream &os);

  // PER SEQUENCE GC CONTENT: columns = GC content, Count
  void write_per_sequence_gc_content(ostream &os);

  // PER BASE N CONTENT: columns = Base, N-count
  void write_per_base_n_content(ostream &os);

  // SEQUENCE LENGTH DISTRIBUTION: Length, count
  void write_sequence_length_distribution(ostream &os);

  // SEQUENCE DUPLICATION LEVELS: Dup_level, perc_dedup, perc_total
  void write_sequence_duplication_levels(ostream &os);

  // OVERREPRESENTED SEQUENCES: Sequence, count, percentage, source
  void write_overrepresented_sequences(ostream &os);

  // ADAPTERT CONTENT: Position, list of adapters
  void write_adapter_content(ostream &os);
};

// Default constructor
FastqStats::FastqStats(const size_t _kmer_size) {
  total_bases = 0;
  avg_read_length = 0;
  avg_gc = 0;
  avg_n = 0;
  num_reads = 0;
  num_poor = 0;

  // Initialize IO arrays
  base_count.fill(0);
  n_base_count.fill(0);
  base_quality.fill(0);
  n_base_quality.fill(0);
  read_length_freq.fill(0);

  // Defines k-mer mask, length and allocates vector
  kmer_size = _kmer_size;
  kmer_mask = (1ll << (2*kmer_size)) - 1;
  kmer_count = vector<size_t>(kmer_mask + 1, 0);
}

// function that makes a hash map, where keys are 32 character suffixes of reads
// and values are all frequent strings that were seen in the first kNumSlowReads
// reads
void
FastqStats::make_suffix_lookup() {
  for (auto v : seq_count) {
    // Candidates have at least half the fraction of frequency given as input
    if (v.second > kMinFracFrequency * kNumSlowReads / 2.0) {
      size_t suffix_hash = get_suffix_hash(v.first);

      // Registers the sequence associated to suffix
      suffix_seq[suffix_hash].push_back(v.first);

      // Registers how many times we have currently seen the string
      suffix_count[suffix_hash].push_back(v.second);

      // Adds the hash to look up
      suffix_lookup.insert(suffix_hash);
    }
  }
}

// Calculates all summary statistics and pass warn fails
void
FastqStats::calc_summary() {

  /******************* BASIC STATISTICS **********************/
  cerr << "Basic statistics\n";
  pass_basic_statistics = "pass";  // in fastqc, basic statistics is always pass
  // Average read length
  avg_read_length = 0;
  for (size_t i = 0; i < kNumBases; ++i)
    avg_read_length += i * read_length_freq[i];
  avg_read_length /= num_reads;

  // GC %
  avg_gc = 0;
  for (size_t i = 0; i < kNumBases; ++i)
    avg_gc += gc_count[i];
  avg_gc /= num_reads;

  // Poor quality reads
  num_poor = 0;
  for (size_t i = 0; i < kPoorQualityThreshold; ++i)
    num_poor += quality_count[i];

  // Cumulative read length frequency
  size_t cumulative_sum = 0;
  for (long long int i = kNumBases - 1; i >= 0; --i) {
    cumulative_sum += read_length_freq[i];
    cumulative_read_length_freq[i] = cumulative_sum;
  }


  /******************* PER BASE SEQUENCE QUALITY **********************/
  cerr << "Per base sequence quality\n";
  pass_per_base_sequence_quality = "pass";

  // Quality quantiles for base positions
  size_t cur;  // for readability, get the quality x position count
  for (size_t i = 0; i < kNumBases; ++i) {
    if (cumulative_read_length_freq[i] > 0) {
      mean[i] = 0;
      size_t counts = 0;
      size_t nreads = cumulative_read_length_freq[i];

      // Number of counts I need to see to know in which bin each *ile is
      double ldecile_thresh = 0.1*nreads;
      double lquartile_thresh = 0.25*nreads;
      double median_thresh = 0.5*nreads;
      double uquartile_thresh = 0.75*nreads;
      double udecile_thresh = 0.9*nreads;

      // Iterate through quality values to find quantiles in increasing order
      for (size_t j = 0; j < kNumQualityValues; ++j) {
        cur = position_quality_count[(i << kBitShiftQuality) + j];

        // Finds in which bin of the histogram reads are
        if (counts < ldecile_thresh && counts + cur >= ldecile_thresh)
          ldecile[i] = j;
        if (counts < lquartile_thresh && counts + cur >= lquartile_thresh)
          lquartile[i] = j;
        if (counts < median_thresh && counts + cur >= median_thresh)
          median[i] = j;
        if (counts < uquartile_thresh && counts + cur >= uquartile_thresh)
          uquartile[i] = j;
        if (counts < udecile_thresh && counts + cur >= udecile_thresh)
          udecile[i] = j;

        mean[i] += cur*j;
        counts += cur;
      }

      // Normalize mean
      mean[i] = mean[i] / nreads;

      // Pass warn fail criteria
      if(pass_per_base_sequence_quality != "fail"){
        if (lquartile[i] < 5)
          pass_per_base_sequence_quality = "fail";
        else if (lquartile[i] < 10)
          pass_per_base_sequence_quality = "warn";

        if (median[i] < 20)
          pass_per_base_sequence_quality = "fail";
        else if (median[i] < 25)
          pass_per_base_sequence_quality = "warn";
 
      }
    }
  }

  /******************* PER SEQUENCE QUALITY SCORE **********************/
  pass_per_sequence_quality_scores = "pass";
  size_t mode_val = 0;
  size_t mode_ind = 0;
  for (size_t i = 0; i < kNumQualityValues; ++i) {
    if (quality_count[i] > mode_val) {
      mode_val = quality_count[i];
      mode_ind = i;
    }
  }

  if (mode_ind < 27)
    pass_per_sequence_quality_scores = "warn";

  else if (mode_ind < 20)
    pass_per_sequence_quality_scores = "fail";

  /******************* PER BASE SEQUENCE CONTENT **********************/
  pass_per_base_sequence_content = "pass";
  size_t a, t, g, c, n;
  double total;
  double max_diff = 0.0;
  for (size_t i = 0; i < kNumBases; ++i) {
    if (cumulative_read_length_freq[i] > 0) {
      a = base_count[(i << kBitShiftNucleotide)];
      c = base_count[(i << kBitShiftNucleotide) + 1];
      t = base_count[(i << kBitShiftNucleotide) + 2];
      g = base_count[(i << kBitShiftNucleotide) + 3];
      n = n_base_count[i];

      // turns above values to percent
      total = static_cast<double>(a+c+t+g+n);
      g_pct[i] = 100.0*(g + 0.25*n) / total;
      a_pct[i] = 100.0*(a + 0.25*n) / total;
      t_pct[i] = 100.0*(t + 0.25*n) / total;
      c_pct[i] = 100.0*(c + 0.25*n) / total;
      n_pct[i] = 100.0*n / total;

      max_diff = max(max_diff, fabs(a-c));
      max_diff = max(max_diff, fabs(a-t));
      max_diff = max(max_diff, fabs(a-g));
      max_diff = max(max_diff, fabs(c-t));
      max_diff = max(max_diff, fabs(c-g));
      max_diff = max(max_diff, fabs(t-g));

      if (pass_per_base_sequence_content != "fail") {
        if (max_diff > 0.2)
          pass_per_base_sequence_content = "fail";
        else if (max_diff > 0.1)
          pass_per_base_sequence_content = "warn";
      }
    }
  }

  /******************* PER BASE N CONTENT **********************/
  pass_per_base_n_content = "pass";
  for (size_t i = 0; i < kNumBases; ++i) {
    if (cumulative_read_length_freq[i] > 0) {
      if (pass_per_base_n_content != "fail") {
        if (n_pct[i] > 20.0) {
          cerr << "n failed at i = " << i << "\t" << n_pct[i] << "\n";
          pass_per_base_n_content = "fail";
        }
        else if (n_pct[i] > 10.0)
          pass_per_base_n_content = "warn";
      }
    }
  }

  /************** SEQUENCE LENGTH DISTRIBUTION *****************/
  pass_sequence_length_distribution = "pass";
  if (read_length_freq[avg_read_length] != num_reads)
    pass_sequence_length_distribution = "warn";
  if (read_length_freq[0] > 0)
    pass_sequence_length_distribution = "fail";


  /************** DUPLICATE SEQUENCES **************************/
  pass_duplicate_sequences = "pass";
  double seq_total = 0, seq_dedup = 0;

  // Seq counts has the exact count for the first kNumSlowReas sequences
  for (auto v : prefix_count) {
    if(v.second < 10){
      percentage_deduplicated[v.second - 1]++;
      percentage_total[v.second - 1] += v.second; 

    // Between 10 and 49
    } else if (v.second >= 10 && v.second < 50) {
      percentage_deduplicated[9]++;
      percentage_total[9] += v.second; 

    // Between 50 and 99
    } else if (v.second >= 10 && v.second < 100) {
      percentage_deduplicated[10]++;
      percentage_total[10] += v.second; 

    // Between 100 and 499
    } else if (v.second >= 100 && v.second < 500) {
      percentage_deduplicated[11]++;
      percentage_total[11] += v.second; 

    // Between 500 and 999
    } else if (v.second >= 500 && v.second < 1000) {
      percentage_deduplicated[12]++;
      percentage_total[12] += v.second; 

    // Between 1000 and 4999
    } else if (v.second >= 1000 && v.second < 5000) {
      percentage_deduplicated[13]++;
      percentage_total[13] += v.second; 

    // Between 5000 and 9999
    } else if (v.second >= 5000 && v.second < 10000) {
      percentage_deduplicated[14]++;
      percentage_total[14]++;

    // Over 10000
    } else {
      percentage_deduplicated[15]++;
      percentage_total[15]++;
    }
    seq_total += v.second;
    seq_dedup += 1.0;
  }

  // Convert to percentage
  for (auto &v : percentage_deduplicated)
    v = 100.0 * v / seq_dedup;  // Percentage of unique sequences in bin

   // Convert to percentage
  for (auto &v : percentage_total)
    v = 100.0 * v / seq_total;  // Percentage of sequences in bin
 
  /************** OVERREPRESENTED SEQUENCES ********************/
  pass_overrepresented_sequences = "pass";
  for (auto v : suffix_lookup) {
    for (size_t i = 0; i < suffix_count[v].size(); ++i)
      if (suffix_count[v][i] > kMinFracFrequency * num_reads)
        overrep_sequences.push_back(make_pair(suffix_seq[v][i]
                                  , suffix_count[v][i]));
  }

  // Sort strings by frequency
  sort(overrep_sequences.begin(), overrep_sequences.end(), 
      [](auto &a, auto &b){
        return a.second > b.second;
      });

  for (auto seq : overrep_sequences) {
    if (pass_overrepresented_sequences != "fail") {
      if (seq.second > 0.01 * num_reads)
        pass_overrepresented_sequences = "fail";
      else if (seq.second > 0.001 * num_reads)
        pass_overrepresented_sequences = "warn";
    }
  }

  /************** ADAPTER CONTENT ******************************/
  pass_adapter_content = "pass";

  /************** KMER CONTENT *********************************/
  pass_kmer_content = "pass";

  /************** PER TILE SEQUENCE QUALITY *********************************/
}

/*************************************************************
 ******************** FASTQ READER ***************************
 *************************************************************/

struct FastqReader{
 private:
  // Memory map variables
  char *curr;  // current position in file
  char *last;  // last position in file
  struct stat st;

  // Temp variables to be updated as you pass through the file
  size_t base_ind;  // 0,1,2 or 3
  size_t read_pos;  // which base we are at in the read
  size_t quality_value;  // to convert from ascii to number
  size_t cur_gc_count;
  size_t cur_quality;
  size_t kmer_pos;
  void *mmap_data;

  // Temporarily store line 2 out of 4 to know the base to which
  // quality characters are associated
  char *buff;

  size_t cur_kmer;  // 32-mer hash as you pass through the sequence line

  /************ FUNCTIONS TO READ LINES IN DIFFERENT WAYS ***********/
  inline void read_fast_forward_line();  // run this to ignore a line
  inline void read_tile_line(FastqStats &stats);  // get tile from read name
  inline void read_sequence_line(FastqStats &stats);  // parse sequence
  inline void read_quality_line(FastqStats &stats);  // parse quality
  inline void process_slow_read(FastqStats &stats);  // hash string
  inline void process_fast_read(FastqStats &stats); 

 public:
  explicit FastqReader(const size_t buffer_size);
  void memorymap(const string &filename,
                 const bool VERBOSE);
  void memoryunmap();
  inline bool operator >> (FastqStats &stats);
};

FastqReader::FastqReader(const size_t buffer_size) {
  base_ind = 0;  // 0,1,2 or 3
  read_pos = 0;  // which base we are at in the read
  cur_gc_count = 0;
  cur_quality = 0;
  kmer_pos = 0;
  cur_kmer = 0;

  // Allocates buffer to temporarily store reads
  buff = new char[buffer_size];
}

// Open file
void
FastqReader::memorymap(const string &filename,
                       const bool VERBOSE) {
  int fd = open(filename.c_str(), O_RDONLY, 0);
  if (fd == -1)
    throw runtime_error("failed to open fastq file: " + filename);

  // get the file size
  fstat(fd, &st);

  // execute mmap
  mmap_data = mmap(NULL, st.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
  close(fd);
  if (mmap_data == MAP_FAILED)
    throw runtime_error("failed to mmap fastq file: " + filename);

  if (VERBOSE) {
    // Display fastq size
    cerr << "Fastq file size: ";
    if (st.st_size > (1ll << 40))
      cerr << st.st_size / ((1ll << 40)) << " Tb\n";
    else if (st.st_size > (1ll << 30))
      cerr << st.st_size / ((1ll << 30)) << " Gb\n";
    else if (st.st_size > (1ll << 20))
      cerr << st.st_size / ((1ll << 20)) << " Mb\n";
    else if (st.st_size > (1ll << 10))
      cerr << st.st_size / ((1ll << 10)) << " Kb\n";
    else
      cerr << st.st_size << " b\n";
  }

  // Initialize position pointer
  curr = static_cast<char*>(mmap_data);
  last = curr + st.st_size - 1;
}

// Gets the tile from the sequence name (if applicable)
inline void
FastqReader::read_tile_line(FastqStats &stats){
  for (; *curr != '\n'; ++curr) {}
  ++curr;  // skips \n from line 1
}

// Skips lines that are not relevant
inline void
FastqReader::read_fast_forward_line(){
  for (; *curr != '\n'; ++curr) {}
  ++curr;  // skips \n from line 1
}

// Reads the line that has the biological sequence
inline void
FastqReader::read_sequence_line(FastqStats &stats){
  read_pos = 0;
  cur_gc_count = 0;

  /*********************************************************
   ********** THIS LOOP MUST BE ALWAYS OPTIMIZED ***********
   *********************************************************/
  for (; *curr != '\n';) {
    // only collect statistics for the first kNumBases of each read
    if (read_pos < stats.kNumBases) {
      // Need to store the character to know which base the quality read is
      // associated to
      buff[read_pos] = *curr;

      // Gets data for ATGC
      if (*curr != 'N') {
        // Transforms base into 3-bit index
        // Bits 2,3 and 4 of charcters A,C,G,T and N are distinct so we can just
        // use them to get the index instead of doing if-elses.
        base_ind = atgc_to_2bit(*curr);

        // Increments count of base
        stats.base_count[(read_pos << stats.kBitShiftNucleotide) | base_ind]++;

        // Hash the current kmer
        cur_kmer = ((cur_kmer << stats.kBitShiftNucleotide) | base_ind);

        // number of gcs in the read, G and C have odd last bit
        cur_gc_count += (base_ind & 1);

        // If we already read >= k bases since the last N, 
        // increment the k-mer count
        if (kmer_pos >= stats.kmer_size - 1)
          stats.kmer_count[cur_kmer & stats.kmer_mask]++;

        kmer_pos++;
      } else {
        // Gets data for N
        stats.n_base_count[read_pos]++;
        kmer_pos = 0;  // start over the current kmer
      }

      // Increment char and read position
      curr++;
      read_pos++;
      
      if(stats.num_reads < stats.kNumSlowReads)
        if(read_pos == stats.kPrefixPosition)
          stats.prefix_count[cur_kmer]++;
      
    } else {
      // Otherwise fast forward subsequent nucleotides
      while (*curr != '\n')
        ++curr;

    }
  }

  ++curr;  // skips \n from line 2
  stats.read_length_freq[read_pos - 1]++;  // read length distribution
  stats.gc_count[100 *cur_gc_count / read_pos]++;  // gc histogram

}

// Reads the quality line of each base. 
inline void
FastqReader::read_quality_line(FastqStats &stats){
  cur_quality = 0;
  read_pos = 0;
  quality_value = 0;
  for (; (*curr != '\n') && (curr < last);) {
    // Converts quality ascii to zero-based
    quality_value = *curr - stats.kBaseQuality;

    if(read_pos < stats.kNumBases){

      // Gets the base from buffer, finds index, adds quality value for
      // that respective base
      if (buff[read_pos] == 'N') {
        stats.n_base_quality[read_pos] += quality_value;
      } else {
        base_ind = atgc_to_2bit(buff[read_pos]);
        stats.base_quality[(read_pos << stats.kBitShiftNucleotide) | base_ind] 
              += quality_value;
      }

      // Sums quality value so we can bin the average at the end
      cur_quality += quality_value;

      // We store for every position x quality value pair the count
      stats.position_quality_count[
            (read_pos << stats.kBitShiftQuality) | quality_value]++;

      ++curr;
      ++read_pos;
    } else {
      // Fast forward if read is too large
      while (*curr != '\n')
        ++curr;
    }
  }
  ++curr;  // skip \n

  // Average quality approximated to the nearest integer. Used to make a
  // histogram in the end of the summary. 
  stats.quality_count[cur_quality / read_pos]++;  // avg quality histogram
}

// Every slow process done in the first kNumSlowReads reads
inline void
FastqReader::process_slow_read(FastqStats &stats){
  buff[read_pos] = '\0';
  stats.seq_count[string(buff)]++;
}

// Every fast lookup done in all remaining reads
inline void
FastqReader::process_fast_read(FastqStats &stats){
  // cur_kmer has the last 32 bases stored. Look it up
  auto got = stats.suffix_lookup.find(cur_kmer);
  if (got != stats.suffix_lookup.end()) {
    buff[read_pos] = '\0';

    for (size_t i = 0; i < stats.suffix_seq[cur_kmer].size(); ++i)

      // Possible optimization is to make this string comparison faster
      if (string(buff) == stats.suffix_seq[cur_kmer][i]) {
        stats.suffix_count[cur_kmer][i]++;
        break;  // get rid of this later
      }
  }
}

inline bool
FastqReader::operator >> (FastqStats &stats) {
  // Get tile number on 1st line
  read_tile_line(stats);
  
  // Process sequence
  read_sequence_line(stats);

  // fast forward third line
  read_fast_forward_line();

  // Process quality values
  read_quality_line(stats);

  // Hash the entire read string to check for overrepresentation
  if (stats.num_reads < stats.kNumSlowReads) {
    process_slow_read(stats);
  }

  // Keep the frequent suffixes to check the future sequences
  else if (stats.num_reads == stats.kNumSlowReads) {
    stats.make_suffix_lookup();
  }

  // Otherwise we just lookup the overrepresented sequence
  else {
    process_fast_read(stats);
  }

  // Successful read, increment number in stats
  stats.num_reads++;

  // Returns if file should keep being checked
  if (curr < last - 1) {
    return true;  // file should keep being read
  }
  return false;  // EOF
}

void
FastqReader::memoryunmap() {
  munmap(mmap_data, st.st_size);
}

/******************************************************
 ********************* WRITE STATS ********************
 ******************************************************/

void
FastqStats::write_header(ostream &os) {
  os << "##FastQC\t0.11.8\n";
}

void
FastqStats::write_basic_statistics(ostream &os, string filename) {
  // TODO(gui): Check what defines pass or fail
  os << ">>Basic Statistics\t" << pass_basic_statistics << "\n";
  os << "#Measure\tValue\n";
  os << "Filename\t" << filename << "\n";
  os << "File type\tConventional base calls\n";
  os << "Total Sequences\t" << num_reads << "\n";
  os << "Sequences flagged as poor quality \t" << num_poor << "\n";
  os << "%GC \t" << avg_gc << "\n";
  os << ">>END_MODULE\n";
}

void
FastqStats::write_per_base_quality(ostream &os) {
  os << ">>Per base sequence quality\t" << 
         pass_per_base_sequence_quality << "\n";

  os << "#Base\tMean\tMedian\tLower Quartile\tUpper Quartile" <<
        "\t10th Percentile 90th Percentile\n";
  for (size_t i = 0; i < kNumBases; ++i) {
    if (cumulative_read_length_freq[i] > 0) {
      // Write distribution to new line
      os << i + 1 << "\t" 
         << mean[i] << "\t" 
         << median[i] << "\t" 
         << lquartile[i] << "\t"
         << uquartile[i] << "\t" 
         << ldecile[i] << "\t" 
         << udecile[i] << "\n";
    }
  }
  os << ">>END_MODULE\n";
}

void
FastqStats::write_per_sequence_quality(ostream &os) {
  os << ">>Per sequence quality scores\t" << 
        pass_per_sequence_quality_scores << "\n";

  os << "#Quality\tCount\n";

  for (size_t i = 0; i < kNumQualityValues; ++i) {
    if (quality_count[i] > 0)
      os << i << "\t" << quality_count[i] << "\n";
  }
  os << ">>END_MODULE\n";
}

void
FastqStats::write_per_base_sequence_content(ostream &os) {
  os << ">>Per base sequence content\t" << 
        pass_per_base_sequence_content << "\n";

  os << "#Base\tG\tA\tT\tC\n";

  for (size_t i = 0; i < kNumBases; ++i) {
    if (cumulative_read_length_freq[i] > 0) {
      os << i+1 << "\t" <<
            g_pct[i] << "\t" <<
            a_pct[i] << "\t" <<
            t_pct[i] << "\t" <<
            c_pct[i] << "\n";
    }
  }
  os << ">>END_MODULE\n";
}

void
FastqStats::write_per_sequence_gc_content(ostream &os) {
  os << ">>Per sequence gc content\t" << pass_per_sequence_gc_content << "\n";
  os << "#GC Content\tCount\n";
  for (size_t i = 0; i <= 100; ++i) {
    if (gc_count[i] > 0) {
      os << i << "\t" << gc_count[i] << "\n";
    }
  }
  os << ">>END_MODULE\n";
}

void
FastqStats::write_per_base_n_content(ostream &os) {
  os << ">>Per base N concent\t" << pass_per_base_n_content << "\n";
  os << "#Base\tN-Count\n";

  for (size_t i = 0; i < kNumBases; ++i) {
    if (cumulative_read_length_freq[i] > 0) 
      os << i+1 << "\t" << n_pct[i] << "\n";
  }

  os << ">>END_MODULE\n";
}

void
FastqStats::write_sequence_length_distribution(ostream &os) {
  os << "Sequence Length Distribution\t" << 
        pass_sequence_length_distribution << "\n";

  os << "Length\tCount\n";
  for (size_t i = 0; i < kNumBases; ++i)
    if (read_length_freq[i] > 0)
      os << i+1 << "\t" << read_length_freq[i] << "\n";
  os << ">>END_MODULE\n";
}

void
FastqStats::write_sequence_duplication_levels(ostream &os) {
  os << ">>Sequence Duplication Levels\t" << 
         pass_sequence_duplication_levels << "\n";

  os << "#Duplication Level  Percentage of deduplicated  Percentage of total\n";
  for(size_t i = 0; i < 9; ++i)
    os << i+1 << "\t" << percentage_deduplicated[i] << "\t" 
       << percentage_total[i] << "\n";

  os << ">10\t" << percentage_deduplicated[9] << "\t" << percentage_total[9] << "\n";
  os << ">50\t" << percentage_deduplicated[10] << "\t" << percentage_total[10] << "\n";
  os << ">100\t" << percentage_deduplicated[11] << "\t" << percentage_total[11] << "\n";
  os << ">500\t" << percentage_deduplicated[12] << "\t" << percentage_total[12] << "\n";
  os << ">1k\t" << percentage_deduplicated[13] << "\t" << percentage_total[13] << "\n";
  os << ">5k\t" << percentage_deduplicated[14] << "\t" << percentage_total[14] << "\n";
  os << ">10k+\t" << percentage_deduplicated[15] << "\t" << percentage_total[15] << "\n";
  os << ">>END_MOUDLE\n";
}

void
FastqStats::write_overrepresented_sequences(ostream &os) {
  os << ">>Overrepresented sequences\t" << 
        pass_overrepresented_sequences << "\n";
  os << "#Sequence\tCount\tPercentage\tPossible Source\n";

  for (auto seq : overrep_sequences)
      os << seq.first << "\t" << seq.second <<  "\t" <<
        100.0 * seq.second / num_reads << "\t???\n";
  os << ">>END_MODULE\n";
}

void
FastqStats::write_adapter_content(ostream &os) {
  os << ">>Adapter Content\t" << pass_adapter_content << "\n";
  os << ">>END_MODULE\n";
}

/******************************************************
 ********************* MAIN ***************************
 ******************************************************/

int main(int argc, const char **argv) {
  clock_t begin = clock();  // register ellapsed time
  /****************** COMMAND LINE OPTIONS ********************/
  string filename;  // fastq file
  string outfile;  // optional filename to save
  bool VERBOSE = false;  // print more run info

  // file containing contaminants to be checked for
  string contaminants_file = "misc/contaminants.tsv";

  // Length of k-mers to count. Grows memory exponentially so it is bound to 10
  size_t kmer_size = 7;
  const size_t MAX_KMER_SIZE = 10;

  OptionParser opt_parse(strip_path(argv[0]),
                         "Quality control metrics for fastq files",
                         "<fastq-file>");

  opt_parse.add_opt("kmer", 'k',
                    "k-mer size (default = 7, max = 10)", false, kmer_size);

  opt_parse.add_opt("outfile", 'o',
                    "filename to save results (default = stdout)",
                    false, outfile);

  opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
  vector<string> leftover_args;

  opt_parse.parse(argc, argv, leftover_args);
  if (argc == 1 || opt_parse.help_requested()) {
    cerr << opt_parse.help_message() << endl
    << opt_parse.about_message() << endl;
    return EXIT_SUCCESS;
  }
  if (opt_parse.about_requested()) {
    cerr << opt_parse.about_message() << endl;
    return EXIT_SUCCESS;
  }
  if (opt_parse.option_missing()) {
    cerr << opt_parse.option_missing_message() << endl;
    return EXIT_SUCCESS;
  }

  if (leftover_args.size() != 1) {
    cerr << opt_parse.help_message() << endl;
    return EXIT_SUCCESS;
  }

  if (kmer_size > MAX_KMER_SIZE) {
    cerr << "K-mer size should not exceed << " << MAX_KMER_SIZE << "\n";
    cerr << opt_parse.help_message() << endl;
    return EXIT_SUCCESS;
  }

  if (kmer_size < 2) {
    cerr << "K-mer size should be smaller than 2\n";
    cerr << opt_parse.help_message() << endl;
    return EXIT_SUCCESS;
  }
  filename = leftover_args.front();

  /****************** END COMMAND LINE OPTIONS *****************/
  // Allocates vectors to summarize data
  FastqStats stats(kmer_size);

  if (VERBOSE)
    cerr << "Started reading file " << filename << ".\n";

  // Initializes a reader given the maximum number of bases to summarie
  FastqReader in (stats.kNumBases);
  in.memorymap(filename, VERBOSE);

  if (VERBOSE)
    cerr << "Started processing file " << filename << ".\n";

  while (in >> stats) {}

  in.memoryunmap();
  if (VERBOSE)
    cerr << "Finished reading file.\n";

  if (VERBOSE)
    cerr << "Summarizing data.\n";

  // This function has to be called before writing to output
  stats.calc_summary();

  /************************ WRITE TO OUTPUT *****************************/
  if (VERBOSE)
    cerr << "Writing data.\n";

  // define output
  ofstream of;
  if (!outfile.empty())
    of.open(outfile.c_str(), ofstream::binary);

  ostream os(outfile.empty() ? cout.rdbuf() : of.rdbuf());

  // Write
  stats.write_header(os);
  stats.write_basic_statistics(os, filename);
  stats.write_per_base_quality(os);
  stats.write_per_sequence_quality(os);
  stats.write_per_base_sequence_content(os);
  stats.write_per_sequence_gc_content(os);
  stats.write_per_base_n_content(os);
  stats.write_sequence_length_distribution(os);
  stats.write_sequence_duplication_levels(os);
  stats.write_overrepresented_sequences(os);
  stats.write_adapter_content(os);

  /************** TIME SUMMARY *********************************/
  // TODO(gui) : find adapters with significant kmer enrichment
  if (VERBOSE)
    cerr << "Elapsed time: "
         << (clock() - begin) / CLOCKS_PER_SEC
         << " seconds\n";
}
