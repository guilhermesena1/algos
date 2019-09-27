/* fqc: quality control for fastq fiCes
*
* Copyright (C) 2019 Guilherme De Sena Brandine and
*                    Andrew D. Smith
* Authors: Guilherme De Sena Brandine, Andrew Smith
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
#include <sstream>

#include <vector>
#include <array>
#include <string>
#include <cstring>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <zlib.h>

#include <htslib/sam.h>

#include <chrono>
#include <ctime>

#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "OptionParser.hpp"

using std::string;
using std::runtime_error;
using std::cerr;
using std::endl;
using std::vector;
using std::array;
using std::find;
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
using std::min;
using std::stringstream;
using std::istringstream;
using std::ostringstream;
using std::put_time;
using std::to_string;

using std::chrono::system_clock;
using std::ctime;

/*************************************************************
 ******************** AUX FUNCTIONS **************************
 *************************************************************/

static inline void
log_process (const string &s) {
  auto tmp = system_clock::to_time_t(system_clock::now());
  string time_fmt = string(ctime(&tmp));
  time_fmt.pop_back();
  cerr << "[" << time_fmt << "] " << s << "\n";
}

// converts 64 bit integer to a sequence string by reading 2 bits at a time and
// converting back to ACTG
static inline string
size_t_to_seq(size_t v, const size_t seq_length) {
  string ans;
  for (size_t i = 0; i < seq_length; ++i) {
    switch (v & 3) {
      case 0: ans.push_back('A'); break;
      case 1: ans.push_back('C'); break;
      case 2: ans.push_back('T'); break;
      case 3: ans.push_back('G'); break;
    }
    v >>= 2;
  }

  reverse(ans.begin(), ans.end());
  return ans;
}

// Converts A,T,G,C to 2-bit values
static inline size_t
actg_to_2bit(const char &c) {
  return ((c >> 1) & 3);
}

// log of a power of two, to use in bit shifting for fast index acces
constexpr size_t
log2exact(size_t powerOfTwo) {
  if (powerOfTwo & (powerOfTwo - 1))
    throw std::runtime_error("not a power of two!");

  size_t ans = 0;
  while (powerOfTwo > 0) {
    ans++;
    powerOfTwo >>= 1;
  }

  return ans - 1;
}

// Check if a string ends with another, to be use to figure out the file format
inline bool 
endswith(std::string const & value, std::string const & ending) {
      if (ending.size() > value.size()) return false;
          return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
}

/******************* IMPLEMENTATION OF FASTQC FUNCTIONS **********************/
// FastQC extrapolation of counts to the full file size
double get_corrected_count (size_t count_at_limit,
                            size_t num_reads,
                            size_t dup_level,
                            size_t num_obs) {

  // See if we can bail out early
  if (count_at_limit == num_reads)
    return num_obs;

  // If there aren't enough sequences left to hide another sequence with this 
  // count the we can also skip the calculation
  if (num_reads - num_obs < count_at_limit) 
    return num_obs;

  // If not then we need to see what the likelihood is that we had 
  // another sequence with this number of observations which we would 
  // have missed. We'll start by working out the probability of NOT seeing a 
  // sequence with this duplication level within the first count_at_limit 
  // sequences of num_obs.  This is easier than calculating
  // the probability of seeing it.
  double p_not_seeing = 1.0;

  // To save doing long calculations which are never going to produce anything 
  // meaningful we'll set a limit to our p-value calculation.  This is the 
  // probability below which we won't increase our count by 0.01 of an 
  // observation.  Once we're below this we stop caring about the corrected 
  // value since it's going to be so close to the observed value thatwe can 
  // just return that instead.
  double limit_of_caring = 1.0 - (num_obs/(num_obs + 0.01));
  for (size_t i = 0; i < count_at_limit; ++i) {
    p_not_seeing *= static_cast<double>((num_reads-i)-dup_level) /
                         static_cast<double>(num_reads-i);

    if (p_not_seeing < limit_of_caring) {
      p_not_seeing = 0;
      break;
    }
  }

  // Now we can assume that the number we observed can be 
  // scaled up by this proportion
  return num_obs/(1 - p_not_seeing);
}

// Function to calculate the deviation of a histogram with 100 bins from a
// theoretical normal distribution with same mean and standard deviation
double 
sum_deviation_from_normal (const array <size_t, 101> &gc_content,
                           array <double, 101> &theoretical) {
  double mode = 0.0,  num_reads = 0.0;

  // Mode and total number of reads (note that we "smoothed" the gc content by
  // filling up the 0 reads, so the number of reads is not the same as the
  // number of reads in the whole fastq
  for (size_t i = 0; i < 101; ++i) {
    mode += i*gc_content[i];
    num_reads += gc_content[i];
  }
  mode /= num_reads;

  // stadard deviation from the modde
  double stdev = 0.0;
  for (size_t i = 0; i < 101; ++i)
    stdev += (mode - i) * (mode - i) * gc_content[i];
  stdev = sqrt(stdev / (num_reads - 1));

  // theoretical sampling from a normal distribution with mean = mode and stdev
  // = stdev to the mode from the sampled gc content from the data
  double ans = 0.0, theoretical_sum = 0.0, z;
  theoretical.fill(0);
  for (size_t i = 0; i < 101; ++i) {
    z = i - mode;
    theoretical[i] = exp(- (z*z)/ (2 * stdev *stdev));
    theoretical_sum += theoretical[i];
  }

  // Normalize theoretical so it sums to the total of reads 
  for (size_t i = 0; i < 101; ++i) {
    theoretical[i] = theoretical[i] * num_reads / theoretical_sum;
  }
  for (size_t i = 0; i < 101; ++i)
    ans += fabs(gc_content[i] - theoretical[i]);
  
  // Fractional deviation
  return 100.0 * ans / num_reads;
}

/*************************************************************
 ******************** CUSTOM CONFIGURATION *******************
 *************************************************************/

// config from options, constants, magic numbers, etc
struct Config {
 public:
  /************************************************************
   *************** MY UNIVERSAL CONSTANTS *********************
   ************************************************************/
  // threshold for a sequence to be considered  poor quality
  size_t kPoorQualityThreshold;

  /************ OVERREPRESENTATION ESTIMTES **********/
  // fraction of the number of slow reads a sequence needs to be seen to be
  // considered a candiate for overrepresentation

  double kOverrepMinFrac;
  /************************************************************
   *************** FASTQC OPTION PARSER************************
   ************************************************************/
  bool casava;  // files from raw casava output
  bool nanopore;  // fast5 format
  bool nofilter;  // if running with --casava don't remove read flagged by casava
  bool extract;  // if set the zipped file will be uncompressed
  bool nogroup;  // disable grouping of bases for reads >50bp
  bool compressed;  // whether or not to inflate file
  bool quiet;
  size_t min_length;  // lower limit in sequence length to be shown in report
  size_t threads;  // number of threads to read multiple files in parallel
  size_t kmer_size;  // kmer size
  string format;  // force file format
  string contaminants_file;  // custom contaminants file
  string adapters_file;  // adapters file
  string limits_file;  // file with limits and options and custom analyses
  string tmpdir;  // dir for temp files when generating report images

  // config on how to handle reads
  bool do_duplication,
       do_kmer,
       do_n_content,
       do_overrepresented,
       do_quality_base,
       do_sequence,
       do_gc_sequence,
       do_quality_sequence,
       do_tile,
       do_adapter,
       do_sequence_length;

  /************************************************************
   *************** FASTQC LIMITS *******************************
   ************************************************************/
  // These will be checked in summary, so hashing is not a probem
  unordered_map <string, unordered_map <string, double> > limits;
  static const vector <string> values_to_check;

  /*************** CONTAMINANTS ******************************/
  vector <pair <string, string> > contaminants;  // first = name, scond = seq
  vector <pair <string, size_t> > adapters;  // kmer of the adapter prefix

  // IO
  string filename;

  /*********** FUNCTIONS TO READ FILES *************/
  Config();  // set magic defaults
  void define_file_format();
  void read_limits();  // populate limits hash map
  void read_adapters();
  void read_contaminants();

  void setup();
  string get_matching_contaminant(string seq) const;
};
const vector<string> Config::values_to_check ({
    "duplication",
    "kmer",
    "n_content",
    "overrepresented",
    "quality_base",
    "sequence",
    "gc_sequence",
    "quality_sequence",
    "tile",
    "sequence_length",
    "adapter",
    "duplication",
    "kmer",
    "n_content",
    "overrepresented",
    "quality_base_lower",
    "quality_base_median",
    "sequence",
    "gc_sequence",
    "quality_sequence",
    "tile",
    "sequence_length",
    "adapter"
  });

// Sets magic numbers
Config::Config (){
   kPoorQualityThreshold = 20;
   kOverrepMinFrac = 0.001;
   casava = false;
   nanopore = false;
   nofilter = false;
   extract = false;
   nogroup = false;
   min_length = 0;
   format = "";
   threads = 1;
   contaminants_file = "Configuration/contaminant_list.txt";
   adapters_file = "Configuration/adapter_list.txt";
   limits_file = "Configuration/limits.txt";
   kmer_size = 7;
   quiet = false;
   tmpdir = ".";
}

void
Config::setup() {
  define_file_format();
  read_limits();
  if (limits["adapter"]["ignore"] == 0.0)
    read_adapters();
  if (limits["adapter"]["ignore"] == 0.0)
    read_contaminants();
}

void
Config::define_file_format() {
  if (format == "") {
    if (endswith(filename, "sam")) {
      format = "sam";
      compressed = false;
    } 
    else if (endswith(filename, "bam")) {
      format = "bam";
      compressed = true;
    }
    else {
      format = "fastq";
      if (endswith(filename, "gz"))
        compressed = true;
      else
        compressed = false;
    }
  }
}

void
Config::read_limits (){
  ifstream in (limits_file);
  if (!in)
    throw runtime_error ("limits file does not exist: " + limits_file);

  // Variables to parse lines
  string line, limit, instruction;
  double value;
  while(getline (in, line)) {

    // Lines with # are comments and should be skipped
    if(line[0] != '#') {
      istringstream iss(line);

      // Every line is a limit, warn/error/ignore and the value
      iss >> limit >> instruction >> value;

      if (find(values_to_check.begin(), values_to_check.end(), limit) 
          == values_to_check.end())
        throw runtime_error ("unknown limit option: " + limit);

      if (instruction != "warn" && 
          instruction != "error" && 
          instruction != "ignore")
        throw runtime_error ("unknown instruction for limit " + limit +
                             ": " + instruction);

      limits[limit][instruction] = value;
    }
  }

  for (auto v : values_to_check) {
    if(limits.count(v) == 0)
      throw runtime_error ("instruction for limit " + v +
                           " not found in file " + limits_file);
  }
  in.close();

  // Get useful data from config that tells us which analyses to skip
  do_duplication = (limits["duplication"]["ignore"] == 0.0);
  do_kmer = (limits["kmer"]["ignore"] == 0.0);
  do_n_content = (limits["n_content"]["ignore"] == 0.0);
  do_overrepresented = (limits["overrepresented"]["ignore"] == 0.0);
  do_quality_base = (limits["quality_base"]["ignore"] == 0.0);
  do_sequence = (limits["sequence"]["ignore"] == 0.0);
  do_gc_sequence = (limits["gc_sequence"]["ignore"] == 0.0);
  do_quality_sequence= (limits["quality_sequence"]["ignore"] == 0.0);
  do_tile = (limits["tile"]["ignore"] == 0.0);
  do_adapter = (limits["adapter"]["ignore"] == 0.0);
  do_sequence_length = (limits["sequence_length"]["ignore"] == 0.0);

}

void 
Config::read_adapters (){
  ifstream in (adapters_file);
  if (!in)
    throw runtime_error("adapter file not found: " + adapters_file);

  string line, _tmp;
  vector<string> line_by_space;
  string adapter_name, adapter_seq;
  size_t adapter_hash;

  // The contaminants file has a space separated name, and the last instance is
  // the biological sequence
  while (getline(in,line)) {
    if(line[0] != '#'){
      adapter_name = "";
      adapter_seq = "";
      istringstream iss(line);
      while(iss >> _tmp) {
        line_by_space.push_back(_tmp);
      }

      if (line_by_space.size() > 1){
        for(size_t i = 0; i < line_by_space.size() - 1; ++i)
          adapter_name += line_by_space[i] + " ";
        adapter_seq = line_by_space.back();

        if(adapter_seq.size() > kmer_size)
          adapter_seq = adapter_seq.substr(0, kmer_size);

        adapter_hash = 0;
        char c;
        for (size_t i = 0; i < adapter_seq.size(); ++i) {
          c = adapter_seq[i];
          if (c != 'A' && c != 'C' && c != 'T' && c != 'G')
            throw runtime_error("Bad adapter (non-ATGC characters): " 
                                + adapter_seq);

          adapter_hash = (adapter_hash << 2) | actg_to_2bit(c);
        }
        adapters.push_back(make_pair(adapter_name, adapter_hash));
      }

      line_by_space.clear();
    }
  }
  in.close();

}

void
Config::read_contaminants(){
  ifstream in (contaminants_file);
  if (!in)
    throw runtime_error("contaminants file not found: " + contaminants_file);

  string line, _tmp;
  vector<string> line_by_space;
  string contaminant_name, contaminant_seq;

  // The contaminants file has a space separated name, and the last instance is
  // the biological sequence
  while (getline(in,line)) {
    if(line[0] != '#'){
      contaminant_name = "";
      contaminant_seq = "";
      istringstream iss(line);
      while(iss >> _tmp) {
        line_by_space.push_back(_tmp);
      }

      if (line_by_space.size() > 1){
        for(size_t i = 0; i < line_by_space.size() - 1; ++i)
          contaminant_name += line_by_space[i] + " ";
        contaminant_seq = line_by_space.back();
        contaminants.push_back(make_pair(contaminant_name,contaminant_seq));
      }

      line_by_space.clear();
    }
  }
  in.close();
}

// Find contaminant with highest overlap with sequence or return "No Hit" if
// there is none
string
Config::get_matching_contaminant (string seq) const {
  size_t best = 0;
  string ret;
  for (auto v : contaminants) {
    if (seq.size() > v.second.size()) {

      // contaminant contained in sequence
      if(seq.find(v.second) != string::npos) {
        if (v.second.size() > best) {
          best = v.second.size();
          ret = v.first;
        }
      }
    } else {
      // sequence contained in contaminant
      if(v.second.find(seq) != string::npos) {
        // In this case this is the best possible match so return it
        return v.first;
      }
    }
  }

  // If any sequence is a match, return the best one
  if (best > 0)
    return ret;
  return "No Hit";
}
/*************************************************************
 ******************** FASTQ STATS ****************************
 *************************************************************/

struct FastqStats {
  // number of bases for static allocation. 
  static const size_t kNumBases = 500;

  // Value to subtract quality characters to get the actual quality value
  static const size_t kBaseQuality = 33;  // The ascii for the lowest quality

  // Smallest power of two that comprises all possible Illumina quality values.
  // Illumina gives qualities from 0 to 40, therefore we set it as 64. Power of
  // is to avoid double pointer jumps and to get indices with bit shifts.
  static const size_t kNumQualityValues = 64;

  // How many possible nucleotides (must be power of 2!)
  static const size_t kNumNucleotides = 4;  // A = 00,C = 01,T = 10,G = 11

  // maximum tile value
  static const size_t kNumMaxTiles = 2048;
  static const size_t kMinTileNumber = 10000;

  // Maximum number of bases for which to do kmer statistics
  static const size_t kKmerMaxBases = 500;

  /************* DUPLICATION ESTIMATES *************/
  // Number of unique sequences to see before stopping counting sequences
  static const size_t kDupUniqueCutoff = 1e5;

  // Maximum read length to store the entire read in memory
  static const size_t kDupReadMaxSize = 75;

  // Prefix size to cut if read length exceeds the value above
  static const size_t kDupReadTruncateSize = 50;
 
  // Kmer size given as input
  static const size_t kmer_size = 7;

  // Bit shifts as instructions for the arrays
  static const size_t kBitShiftNucleotide = log2exact(kNumNucleotides);
  static const size_t kBitShiftTile = log2exact(kNumMaxTiles);
  static const size_t kBitShiftQuality = log2exact(kNumQualityValues);
  static const size_t kBitShiftKmer = 2 * kmer_size; // two bits per base

 public:
  /*********** SINGLE NUMBERS FROM THE ENTIRE FASTQ ****************/
    // Number of unique sequences seen thus far
  size_t num_unique_seen;

  // How many reads were processed before num_unique_seen = kDupUniqueCutoff
  size_t count_at_limit;

  size_t total_bases;  // sum of all bases in all reads
  size_t avg_read_length;  // average of all read lengths
  size_t num_reads;  // total number of lines read
  size_t num_reads_kmer; // number of reads in which kmer counts was performed
  size_t num_reads_tile; // number of reads in which tile count was performed
  size_t min_read_length; // minimum read length seen
  size_t max_read_length;  // total number of lines read
  size_t num_poor;  // reads whose average quality was <= poor
  size_t num_extra_bases;  // number of bases outside of buffer

  // mask to get only the first 2*k bits of the sliding window
  size_t kmer_mask;

  double avg_gc;  // (sum of g bases + c bases) / (num_reads)
  double total_deduplicated_pct;  // number of reads left if deduplicated

  /*********************************************************
   *********** METRICS COLLECTED DURING IO *****************
   *********************************************************/
  /*********** PER BASE METRICS ****************/

  // counts the number of bases in every read position
  array<size_t, kNumNucleotides * kNumBases> base_count;  // ATGC
  array<size_t, kNumNucleotides * kNumBases> n_base_count; // N

  /*********** PER QUALITY VALUE METRICS ****************/
  // Counts of quality in each base position
  array<size_t, kNumQualityValues * kNumBases> position_quality_count;

  // Counts of average quality (truncated) per sequence
  array<size_t, kNumQualityValues> quality_count;

  /*********** PER GC VALUE METRICS ****************/
  // histogram of GC fraction in each read from 0 to 100%
  array<size_t, 101> gc_count;
  array<double, 101> theoretical_gc_count;

    /*********** PER READ METRICS ***************/
  // Distribution of read lengths
  array<size_t, kNumBases> read_length_freq;
  array<size_t, kNumBases> cumulative_read_length_freq;

  /*********** PER TILE SEQUENCE QUALITY OVERSERQUENCES ********/
  array <double, kNumMaxTiles * kNumBases> tile_position_quality;
  array <size_t, kNumMaxTiles> tile_count;

  /*********** SUMMARY ARRAYS **************/
  // Quantiles for the position_quality_count
  array<size_t, kNumBases> ldecile, lquartile, median, uquartile, udecile;
  array<double, kNumBases> mean;

  // For sequence duplication levels
  // 1 to 9, >10, >50, >100, >500, >1k, >5k, >10k+
  array<double, 16> percentage_deduplicated;
  array<double, 16> percentage_total;

  // Percentages for per base sequence content
  array<double, kNumBases> a_pct,
                           c_pct,
                           t_pct,
                           g_pct,
                           n_pct;

  /*********** SLOW STUFF *******************/
  // Leftover memory using dynamic allocation
  vector<size_t> long_base_count;
  vector<size_t> long_n_base_count;
  vector<size_t> long_position_quality_count;
  vector<size_t> long_read_length_freq;
  vector<size_t> long_cumulative_read_length_freq;
  vector<double> long_tile_position_quality;
  vector<size_t> long_ldecile, long_lquartile, long_median,
                 long_uquartile,long_udecile;
  vector<double> long_mean;
  vector<double> long_a_pct,
                 long_c_pct,
                 long_t_pct,
                 long_g_pct,
                 long_n_pct;


  /********** KMER FREQUENCY ****************/
  // A 2^K + 1 vector to count all possible kmers
  vector<size_t> kmer_count;

  /********** ADAPTER COUNT *****************/
  // Percentage of times we saw each adapter in each position
  unordered_map <size_t, vector <double>> kmer_by_base;
  /*********** DUPLICATION ******************/
  unordered_map <string, size_t> sequence_count;

  /*********** OVERREPRESENTED SERQUENCES ********/
  vector <pair<string,size_t>> overrep_sequences;

  /*********************************************************
   *********** METRICS SUMMARIZED AFTER IO *****************
   *********************************************************/

  // I need this to know what to divide each base by
  // when averaging content, bases, etc. It stores, for every element i, how
  // many reads are of length >= i, ie, how many reads actually have a
  // nucleotide at position i

  /*********** PASS WARN FAIL MESSAGE FOR EACH METRIC **************/
  string pass_basic_statistics,
         pass_per_base_sequence_quality,
         pass_per_tile_sequence_quality,
         pass_per_sequence_quality_scores,
         pass_per_base_sequence_content,
         pass_per_sequence_gc_content,
         pass_per_base_n_content,
         pass_sequence_length_distribution,
         pass_overrepresented_sequences,
         pass_duplicate_sequences,
         pass_adapter_content;

    /**************** FUNCTIONS ****************************/

  // Default constructor that zeros everything
  explicit FastqStats();

  // Allocation of more read positions
  inline void allocate_new_base(const bool ignore_tile);

  /******* DUPLICATION AND OVERREPRESENTATION *******/
  // Makes a hash map with keys as 32-bit suffixes and values as all the
  // candidate frequent sequences with that given suffix

  // Summarize all statistics we need before writing
  void summarize(Config &config);

  // Writes to outpute fastqc-style
  void write(ostream &os, const Config &config);
};

// Default constructor
FastqStats::FastqStats() {
  total_bases = 0;
  num_extra_bases = 0;
  avg_read_length = 0;
  avg_gc = 0;
  num_reads = 0;
  num_reads_kmer = 0;
  num_reads_tile = 0;
  min_read_length = 0;
  max_read_length = 0;
  num_poor = 0;

  num_unique_seen = 0;
  count_at_limit = 0;

  // Initialize IO arrays
  base_count.fill(0);
  n_base_count.fill(0);
  read_length_freq.fill(0);
  quality_count.fill(0);
  gc_count.fill(0);
  position_quality_count.fill(0);
  tile_position_quality.fill(0);
  tile_count.fill(0);

  // Defines k-mer mask, length and allocates vector
  kmer_mask = (1ll << (2*kmer_size)) - 1;
  kmer_count = vector<size_t>(min(kNumBases, kKmerMaxBases)
                                  * (kmer_mask + 1), 0);
}

// When we read new bases, dynamically allocate new space for their statistics
inline void
FastqStats::allocate_new_base(const bool ignore_tile) {
  
  for (size_t i = 0; i < kNumNucleotides; ++i)
    long_base_count.push_back(0);
  long_n_base_count.push_back (0);

  // space for quality boxplot
  for (size_t i = 0; i < kNumQualityValues; ++i)
    long_position_quality_count.push_back(0);

  long_read_length_freq.push_back(0);

  // space for tile quality in each position. This takes A LOT of memory so 
  // if tiles are not being considered we will just ignore this
  if(!ignore_tile) {
    for (size_t i = 0; i < kNumMaxTiles; ++i)
      long_tile_position_quality.push_back(0);
  }

  // Successfully allocated space for a new base
  ++num_extra_bases;
}

// Calculates all summary statistics and pass warn fails
void
FastqStats::summarize(Config &config) {
  /******************* BASIC STATISTICS **********************/
  pass_basic_statistics = "pass";  // in fastqc, basic statistics is always pass
  // Average read length
  avg_read_length = 0;
  for (size_t i = 0; i < max_read_length; ++i) {
    if (i < kNumBases)
      total_bases += i * read_length_freq[i];
    else
      total_bases += i * long_read_length_freq[i - kNumBases];
  }

  avg_read_length = total_bases / num_reads;

  // counts bases G and C in each base position
  avg_gc = 0;
  for (size_t i = 0; i < max_read_length; ++i) {
    if (i < kNumBases) {
      avg_gc += base_count[(i << kBitShiftNucleotide) | 1];  // C
      avg_gc += base_count[(i << kBitShiftNucleotide) | 3];  // G
    } else {
        avg_gc += long_base_count[(i - kNumBases) << kBitShiftNucleotide | 1];  // C
      avg_gc += long_base_count[(i - kNumBases) << kBitShiftNucleotide | 3];  // G
    }
  }
  
  // GC %
  avg_gc = 100 * avg_gc / total_bases;

  // Poor quality reads
  num_poor = 0;
  for (size_t i = 0; i < config.kPoorQualityThreshold; ++i)
    num_poor += quality_count[i];

  // Cumulative read length frequency
  size_t cumulative_sum = 0;
  for (size_t i = 0; i < max_read_length; ++i) {
    if (i < kNumBases) {
      cumulative_sum += read_length_freq[i];
      if (read_length_freq[i] > 0)
        if (min_read_length == 0)
          min_read_length = i + 1;
    }
    else
      cumulative_sum += long_read_length_freq[i - kNumBases];
  }

  for (size_t i = 0; i < max_read_length; ++i) {
    if (i < kNumBases) {
      cumulative_read_length_freq[i] = cumulative_sum;
      cumulative_sum -= read_length_freq[i];
    }
    else {
      long_cumulative_read_length_freq.push_back(cumulative_sum);
      cumulative_sum -= long_read_length_freq[i - kNumBases];
    }
  }

  /******************* PER BASE SEQUENCE QUALITY **********************/
  if (config.do_quality_base) {
    pass_per_base_sequence_quality = "pass";

    // Quality quantiles for base positions
    size_t cur;  // for readability, get the quality x position count
    double ldecile_thresh, 
           lquartile_thresh,
           median_thresh,
           uquartile_thresh,
           udecile_thresh;

    size_t cur_ldecile = 0,
           cur_lquartile = 0,
           cur_median = 0,
           cur_uquartile = 0,
           cur_udecile = 0;

    double cur_mean;
    for (size_t i = 0; i < max_read_length; ++i) {
      cur_mean = 0;
      size_t counts = 0;

      // Number of counts I need to see to know in which bin each *ile is
      if (i < kNumBases) {
        ldecile_thresh = 0.1 * cumulative_read_length_freq[i];
        lquartile_thresh = 0.25 * cumulative_read_length_freq[i];
        median_thresh = 0.5 * cumulative_read_length_freq[i];
        uquartile_thresh = 0.75 * cumulative_read_length_freq[i];
        udecile_thresh = 0.9 * cumulative_read_length_freq[i];
      } else {
        ldecile_thresh = 0.1 * long_cumulative_read_length_freq[i - kNumBases];
        lquartile_thresh = 0.25 * long_cumulative_read_length_freq[i - kNumBases];
        median_thresh = 0.5 * long_cumulative_read_length_freq[i - kNumBases];
        uquartile_thresh = 0.75 * long_cumulative_read_length_freq[i - kNumBases];
        udecile_thresh = 0.9 * long_cumulative_read_length_freq[i - kNumBases];
      }
      // Iterate through quality values to find quantiles in increasing order
      for (size_t j = 0; j < kNumQualityValues; ++j) {
        if (i < kNumBases)
          cur = position_quality_count[(i << kBitShiftQuality) | j];
        else
          cur = long_position_quality_count[((i - kNumBases) << kBitShiftQuality) | j];

        // Finds in which bin of the histogram reads are
        if (counts < ldecile_thresh && counts + cur >= ldecile_thresh)
          cur_ldecile = j;

        if (counts < lquartile_thresh && counts + cur >= lquartile_thresh)
          cur_lquartile = j;

        if (counts < median_thresh && counts + cur >= median_thresh)
          cur_median = j;

        if (counts < uquartile_thresh && counts + cur >= uquartile_thresh)
          cur_uquartile = j;

        if (counts < udecile_thresh && counts + cur >= udecile_thresh)
          cur_udecile = j;

        cur_mean += cur*j;
        counts += cur;
      }

      // Normalize mean
      if (i < kNumBases)
        cur_mean = cur_mean / cumulative_read_length_freq[i];
      else
        cur_mean = cur_mean / long_cumulative_read_length_freq[i - kNumBases];

      if (i < kNumBases) {
        mean[i] = cur_mean;
        ldecile[i] = cur_ldecile;
        lquartile[i] = cur_lquartile;
        median[i] = cur_median;
        uquartile[i] = cur_uquartile;
        udecile[i] = cur_udecile;
      } else {
        long_mean.push_back(cur_mean);
        long_ldecile.push_back(cur_ldecile);
        long_lquartile.push_back(cur_lquartile);
        long_median.push_back(cur_median);
        long_uquartile.push_back(cur_uquartile);
        long_udecile.push_back(cur_udecile);
      }

      // Pass warn fail criteria
      if (pass_per_base_sequence_quality != "fail") {
        if (cur_lquartile < config.limits["quality_base_lower"]["error"])
          pass_per_base_sequence_quality = "fail";
        else if (cur_lquartile < config.limits["quality_base_lower"]["warn"])
          pass_per_base_sequence_quality = "warn";

        if (cur_median < config.limits["quality_base_median"]["error"])
          pass_per_base_sequence_quality = "fail";
        else if (cur_median < config.limits["quality_base_median"]["warn"])
          pass_per_base_sequence_quality = "warn";
      }
    }
  }
  /******************* PER SEQUENCE QUALITY SCORE **********************/
  if (config.do_quality_sequence) {
    pass_per_sequence_quality_scores = "pass";
    size_t mode_val = 0;
    size_t mode_ind = 0;
    for (size_t i = 0; i < kNumQualityValues; ++i) {
      if (quality_count[i] > mode_val) {
        mode_val = quality_count[i];
        mode_ind = i;
      }
    }

    if (mode_ind < config.limits["quality_sequence"]["warn"])
      pass_per_sequence_quality_scores = "warn";

    else if (mode_ind < config.limits["quality_sequence"]["error"])
      pass_per_sequence_quality_scores = "fail";
  }
  /******************* PER BASE SEQUENCE CONTENT **********************/
  if (config.do_sequence) {
    pass_per_base_sequence_content = "pass";
    double a, t, g, c, n;
    double total;
    double max_diff = 0.0;
    for (size_t i = 0; i < max_read_length; ++i) {
      if (i < kNumBases) {
        a = base_count[(i << kBitShiftNucleotide)];
        c = base_count[(i << kBitShiftNucleotide) | 1];
        t = base_count[(i << kBitShiftNucleotide) | 2];
        g = base_count[(i << kBitShiftNucleotide) | 3];
        n = n_base_count[i];
      } else {
        a = long_base_count[((i - kNumBases) << kBitShiftNucleotide)];
        c = long_base_count[((i - kNumBases) << kBitShiftNucleotide) | 1];
        t = long_base_count[((i - kNumBases) << kBitShiftNucleotide) | 2];
        g = long_base_count[((i - kNumBases) << kBitShiftNucleotide) | 3];
        n = long_n_base_count[i - kNumBases];
      }

      // turns above values to percent
      total = static_cast<double>(a + c + t + g + n);
      a = 100.0*a / total;
      c = 100.0*c / total;
      t = 100.0*t / total;
      g = 100.0*g / total;
      n = 100.0*n / total;
      if (i < kNumBases) {
        g_pct[i] = g;
        a_pct[i] = a;
        t_pct[i] = t;
        c_pct[i] = c;
        n_pct[i] = n;
      } else {
        long_g_pct.push_back(g);
        long_a_pct.push_back(a);
        long_t_pct.push_back(t);
        long_c_pct.push_back(c);
        long_n_pct.push_back(n);
      }

      max_diff = max(max_diff, fabs(a-c));
      max_diff = max(max_diff, fabs(a-t));
      max_diff = max(max_diff, fabs(a-g));
      max_diff = max(max_diff, fabs(c-t));
      max_diff = max(max_diff, fabs(c-g));
      max_diff = max(max_diff, fabs(t-g));

      if (pass_per_base_sequence_content != "fail") {
        if (max_diff > config.limits["sequence"]["error"])
          pass_per_base_sequence_content = "fail";
        else if (max_diff > config.limits["sequence"]["warn"])
          pass_per_base_sequence_content = "warn";
      }
    }
  }

  /******************* PER SEQUENCE GC CONTENT *****************/
  if (config.do_gc_sequence) {
    pass_per_sequence_gc_content = "pass";
    // fix the zero gcs by averaging around the adjacent values
    for (size_t i = 1; i < 99; ++i) 
      if (gc_count[i] == 0)
        gc_count[i] = (gc_count[i+1] + gc_count[i-1])/2;

    // Calculate pass warn fail statistics
    double gc_deviation = sum_deviation_from_normal(gc_count, 
                                                    theoretical_gc_count);
    if (gc_deviation >= config.limits["gc_sequence"]["error"])
      pass_per_sequence_gc_content = "fail";
    else if (gc_deviation >= config.limits["gc_sequence"]["warn"])
      pass_per_sequence_gc_content = "warn";
  }
  /******************* PER BASE N CONTENT **********************/
  if (config.do_n_content) {
    pass_per_base_n_content = "pass";
    double cur_n_pct;
    for (size_t i = 0; i < max_read_length; ++i) {
      if (pass_per_base_n_content != "fail") {
        if (i < kNumBases)
          cur_n_pct = n_pct[i];
        else
          cur_n_pct = long_n_pct[i - kNumBases];

        if (cur_n_pct > config.limits["n_content"]["error"]) {
          pass_per_base_n_content = "fail";
        } else if (cur_n_pct > config.limits["n_content"]["warn"])
          pass_per_base_n_content = "warn";
      }
    }
  }
  /************** SEQUENCE LENGTH DISTRIBUTION *****************/
  if (config.do_sequence_length) {
    pass_sequence_length_distribution = "pass";
    size_t freq_of_avg;
    
    if (avg_read_length < kNumBases)
      freq_of_avg = read_length_freq[avg_read_length];
    else
      freq_of_avg = long_read_length_freq[avg_read_length - kNumBases];

    if (config.limits["sequence_length"]["error"] == 1) {
      if (freq_of_avg != num_reads)
        pass_sequence_length_distribution = "warn";
      if (read_length_freq[0] > 0)
        pass_sequence_length_distribution = "fail";
    }
  }
  /************** DUPLICATE SEQUENCES **************************/
  if (config.do_duplication) {
    pass_duplicate_sequences = "pass";

    double seq_total = 0.0;
    double seq_dedup = 0.0;

    // Key is frequenccy (r), value is number of times we saw a sequence
    // with that frequency
    unordered_map <size_t, size_t> counts_by_freq;
    for (auto v : sequence_count) {
      if(counts_by_freq.count(v.second) == 0)
        counts_by_freq[v.second] = 0;
      counts_by_freq[v.second]++;
    }

    // Now we run the fastqc corrected extrapolation
    for (auto v: counts_by_freq) {
      counts_by_freq[v.first] = get_corrected_count (count_at_limit,
                                                     num_reads,
                                                     v.first,
                                                     v.second);
    }

    // Group in blocks similarly to fastqc
    for (auto v : counts_by_freq) {
      size_t dup_slot = v.first - 1;
      if (v.first >= 10000) dup_slot = 15;
      else if (v.first >= 5000) dup_slot = 14;
      else if (v.first >= 1000) dup_slot = 13;
      else if (v.first >= 500) dup_slot = 12;
      else if (v.first >= 100) dup_slot = 11;
      else if (v.first >= 50) dup_slot = 10;
      else if (v.first >= 10) dup_slot = 9;

      percentage_deduplicated[dup_slot] += v.second;
      percentage_total[dup_slot] += v.second * v.first;

      seq_total += v.second * v.first;
      seq_dedup += v.second;
    }

    // "Sequence duplication estimate" in the summary
    total_deduplicated_pct = 100.0 * seq_dedup / seq_total;

    // Convert to percentage
    for (auto &v : percentage_deduplicated)
      v = 100.0 * v / seq_dedup;  // Percentage of unique sequences in bin

     // Convert to percentage
    for (auto &v : percentage_total)
      v = 100.0 * v / seq_total;  // Percentage of sequences in bin

    // pass warn fail criteria : unique reads must be >80%
    // (otherwise warn) or >50% (otherwisefail)
    if (percentage_total[0] <= config.limits["duplication"]["error"])
      pass_duplicate_sequences = "fail";

    else if (percentage_total[0] <= config.limits["duplication"]["warn"])
      pass_duplicate_sequences = "warn";
  }
  /************** OVERREPRESENTED SEQUENCES ********************/
  if (config.do_overrepresented) {
    pass_overrepresented_sequences = "pass";

    // Keep only sequences that pass the input cutoff
    for (auto it = sequence_count.begin(); it != sequence_count.end(); ++it) {
      if (it->second > num_reads * config.kOverrepMinFrac) {
        overrep_sequences.push_back (*it);
      }

      // implment pass warn fail for overrep sequences
      if (pass_overrepresented_sequences != "fail") {
        // get percentage that overrep reads represent
        double pct = 100.0 * it->second / num_reads;
        if (pct > config.limits["overrepresented"]["error"])
          pass_overrepresented_sequences = "fail";
        else if (pct > config.limits["overrepresented"]["warn"])
          pass_overrepresented_sequences = "warn";
      }
    }

    // Sort strings by frequency
    sort(overrep_sequences.begin(), overrep_sequences.end(),
         [](auto &a, auto &b){
            return a.second > b.second;
         });
  }
  /************** ADAPTER CONTENT ******************************/
  if (config.do_adapter) {
    pass_adapter_content = "pass";

    // Cumulative count of adapter kmers by position
    size_t jj;
    for (size_t i = 0; i < min(kNumBases, kKmerMaxBases); ++i) {
      if (cumulative_read_length_freq[i] > 0) {

        // Makes the count of kmers by position cumulative by copying
        // the previous position count
        if (i == 0)
          kmer_by_base[i] = vector<double> (config.adapters.size(), 0.0);
        else
          kmer_by_base[i] = vector<double> (kmer_by_base[i-1].begin(),
                                            kmer_by_base[i-1].end());
        jj = 0;

        // Get count for adapter's k-prefix 
        for (auto v : config.adapters) {
          size_t kmer_ind = (i << kBitShiftKmer) | v.second;
          kmer_by_base[i][jj] += kmer_count[kmer_ind];
          ++jj;
        }
      }
    }

    for (size_t i = 0; i < min(kNumBases, kKmerMaxBases); ++i) {
      if (cumulative_read_length_freq[i] > 0) {
        jj = 0;
        for (auto v : config.adapters) {

          kmer_by_base[i][jj] = kmer_by_base[i][jj] * 100.0 / num_reads_kmer;

          // Update pass warn fail
          if (pass_adapter_content != "fail") {
            if (kmer_by_base[i][jj] > config.limits["adapter"]["error"])
              pass_adapter_content = "fail";
            else if (kmer_by_base[i][jj] > config.limits["adapter"]["warn"])
              pass_adapter_content = "warn";
          }
          ++jj;
        }
      }
    }
  }

  /************** PER TILE SEQUENCE QUALITY ********************/
  if (config.do_tile) {
    pass_per_tile_sequence_quality = "pass";

    // Normalize by average: We can reuse the mean we calculated in per base
    // quality distributions above.
    size_t tile_ind;
    for(size_t i = 0; i < max_read_length; ++i) {
      for (size_t j = 0; j < kNumMaxTiles; ++j) {
        if(tile_count[j] > 0){
          if (i < kNumBases) {
            tile_ind = (i << kBitShiftTile) | j;
            tile_position_quality[tile_ind] =
            tile_position_quality[tile_ind] / tile_count[j] - mean[i];

            if(pass_per_tile_sequence_quality != "fail") {
              if (tile_position_quality[tile_ind] <= config.limits["tile"]["error"])
                pass_per_tile_sequence_quality = "fail";

              else if (tile_position_quality[tile_ind] <= config.limits["tile"]["warn"])
                pass_per_tile_sequence_quality = "warn";
            }
          }
          else {
            tile_ind = ((i - kNumBases) << kBitShiftTile) | j;
            long_tile_position_quality[tile_ind] =
            (long_tile_position_quality[tile_ind] / tile_count[j]) - long_mean[i - kNumBases];

            if(pass_per_tile_sequence_quality != "fail") {
              if (long_tile_position_quality[tile_ind] <=
                  -config.limits["tile"]["error"])
                pass_per_tile_sequence_quality = "fail";
              else if (long_tile_position_quality[tile_ind] <=
                  -config.limits["tile"]["warn"])
                pass_per_tile_sequence_quality = "warn";
            }
          }
        }
      }
    }
  }
}

/****************** WRITE STATS ***********************/
void
FastqStats::write(ostream &os, const Config &config) {
  // Header
  os << "##FastQC\t0.11.8\n";

  // Basic statistics
  os << ">>Basic Statistics\t" << pass_basic_statistics << "\n";
  os << "#Measure\tValue\n";
  os << "Filename\t" << strip_path (config.filename) << "\n";
  os << "Total Sequences\t" << num_reads << "\n";
  os << "Sequences flagged as poor quality \t" << num_poor << "\n";

  if (config.do_sequence)
    os << "%GC \t" << avg_gc << "\n";
  os << ">>END_MODULE\n";

  // Per base quality
  if (config.do_quality_base) {
    os << ">>Per base sequence quality\t" <<
           pass_per_base_sequence_quality << "\n";

    os << "#Base\tMean\tMedian\tLower Quartile\tUpper Quartile" <<
          "\t10th Percentile 90th Percentile\n";
    for (size_t i = 0; i < max_read_length; ++i) {
      if (i < kNumBases) {
        // Write distribution to new line
        os << i + 1 << "\t"
           << mean[i] << "\t"
           << median[i] << "\t"
           << lquartile[i] << "\t"
           << uquartile[i] << "\t"
           << ldecile[i] << "\t"
           << udecile[i] << "\n";
      } else {
        os << i + 1 << "\t"
           << long_mean[i - kNumBases] << "\t"
           << long_median[i - kNumBases] << "\t"
           << long_lquartile[i - kNumBases] << "\t"
           << long_uquartile[i - kNumBases] << "\t"
           << long_ldecile[i - kNumBases] << "\t"
           << long_udecile[i - kNumBases] << "\n";
      }
    }
    os << ">>END_MODULE\n";
  }

  // Per sequence quality scores
  if (config.do_quality_sequence) {
    os << ">>Per sequence quality scores\t" <<
          pass_per_sequence_quality_scores << "\n";

    os << "#Quality\tCount\n";

    for (size_t i = 0; i < kNumQualityValues; ++i) {
      if (quality_count[i] > 0)
        os << i << "\t" << quality_count[i] << "\n";
    }
    os << ">>END_MODULE\n";

    // Per base sequence content
    os << ">>Per base sequence content\t" <<
          pass_per_base_sequence_content << "\n";

    os << "#Base\tG\tA\tT\tC\n";

    for (size_t i = 0; i < max_read_length; ++i) {
      if (i < kNumBases) {
        os << i+1 << "\t" <<
              g_pct[i] << "\t" <<
              a_pct[i] << "\t" <<
              t_pct[i] << "\t" <<
              c_pct[i] << "\n";
      } else {
        os << i+1 << "\t" <<
              long_g_pct[i - kNumBases] << "\t" <<
              long_a_pct[i - kNumBases] << "\t" <<
              long_t_pct[i - kNumBases] << "\t" <<
              long_c_pct[i - kNumBases] << "\n";
      }
    }
    os << ">>END_MODULE\n";
  }

  // Per tile sequence quality
  if (config.do_tile) {
    os << ">>Per tile sequence quality\t" <<
          pass_per_tile_sequence_quality << "\n";

    size_t tile_ind;
    for (size_t i = 0; i < kNumMaxTiles; ++i) {
      if (tile_count[i] > 0) {
        for (size_t j = 0; j < max_read_length; ++j) {
          if (j < kNumBases) {
            tile_ind = (j << kBitShiftTile) | i;
            os << i + kMinTileNumber << "\t" << j + 1 << "\t" <<
              tile_position_quality[tile_ind];
          } else {
            tile_ind = ((j - kNumBases) << kBitShiftTile) | i;
            os << i + kMinTileNumber << "\t" << j + 1 << "\t" << 
              long_tile_position_quality[tile_ind];
          }
          os << "\n";
        }
      }
    }

    os << ">>END_MODULE\n";
  }

  // Per sequence gc content
  if (config.do_gc_sequence) {
    os << ">>Per sequence gc content\t" << pass_per_sequence_gc_content << "\n";
    os << "#GC Content\tCount\n";
    for (size_t i = 0; i <= 100; ++i) {
      if (gc_count[i] > 0) {
        os << i << "\t" << gc_count[i] << "\n";
      }
    }
    os << ">>END_MODULE\n";

    // Per base N content
    os << ">>Per base N concent\t" << pass_per_base_n_content << "\n";
    os << "#Base\tN-Count\n";

    for (size_t i = 0; i < max_read_length; ++i) {
      if (i < kNumBases)
        os << i+1 << "\t" << n_pct[i] << "\n";
      else
        os << i+1 << "\t" << long_n_pct[i - kNumBases] << "\n";
    }

    os << ">>END_MODULE\n";
  }

  // Sequence length distribution
  if (config.do_sequence_length) {
    os << "Sequence Length Distribution\t" <<
          pass_sequence_length_distribution << "\n";

    os << "Length\tCount\n";
    for (size_t i = 0; i < max_read_length; ++i) {
      if(i < kNumBases) {
        if (read_length_freq[i] > 0)
          os << i+1 << "\t" << read_length_freq[i] << "\n";
      } else {
        if (long_read_length_freq[i - kNumBases] > 0)
          os << i + 1  << "\t" << long_read_length_freq[i - kNumBases] << "\n";
      }
    }
    os << ">>END_MODULE\n";
  }

  // Sequence duplication levels
  if (config.do_duplication) {
    os << ">>Sequence Duplication Levels\t" <<
           pass_duplicate_sequences << "\n";

    os << ">>Total Deduplicated Percentage\t" <<
           total_deduplicated_pct << "\n";

    os << "#Duplication Level  Percentage of deduplicated  Percentage of total\n";
    for(size_t i = 0; i < 9; ++i)
      os << i+1 << "\t" << percentage_deduplicated[i] << "\t"
         << percentage_total[i] << "\n";

    os << ">10\t" << percentage_deduplicated[9]
       << "\t" << percentage_total[9] << "\n";
    os << ">50\t" << percentage_deduplicated[10]
       << "\t" << percentage_total[10] << "\n";
    os << ">100\t" << percentage_deduplicated[11]
       << "\t" << percentage_total[11] << "\n";
    os << ">500\t" << percentage_deduplicated[12]
       << "\t" << percentage_total[12] << "\n";
    os << ">1k\t" << percentage_deduplicated[13]
       << "\t" << percentage_total[13] << "\n";
    os << ">5k\t" << percentage_deduplicated[14]
       << "\t" << percentage_total[14] << "\n";
    os << ">10k+\t" << percentage_deduplicated[15]
       << "\t" << percentage_total[15] << "\n";
    os << ">>END_MOUDLE\n";
  }

  // Overrepresented sequences
  if (config.do_overrepresented) {
    os << ">>Overrepresented sequences\t" <<
          pass_overrepresented_sequences << "\n";
    os << "#Sequence\tCount\tPercentage\tPossible Source\n";

    for (auto seq : overrep_sequences)
        os << seq.first << "\t" << seq.second <<  "\t" <<
          100.0 * seq.second / num_reads << "\t"
          << config.get_matching_contaminant(seq.first) << "\n";
    os << ">>END_MODULE\n";
    os << ">>Adapter Content\t" << pass_adapter_content << "\n";

    os << "#Position\t";
    for (auto v : config.adapters)
      os << v.first << "\t";

    os << "\n";

    // Number of kmers counted
    size_t jj;
    for (size_t i = 0; i < min(kNumBases, kKmerMaxBases); ++i) {
      if (cumulative_read_length_freq[i] > 0) {
        os << i + 1 << "\t";
        jj = 0;
        for (auto v : config.adapters) {
          os << 100.0 * kmer_by_base[i][jj];
          ++jj;

          if(jj != config.adapters.size())
            os << " ";
        }
          os << "\n";
      }
    }
    os << ">>END_MODULE\n";
  }
}

/*************************************************************
 ******************** STREAM READER **************************
 *************************************************************/

// Generic class that does as much as possible without assuming file format
class StreamReader{
 public:
  // config on how to handle reads
  const bool do_duplication,
             do_kmer,
             do_n_content,
             do_overrepresented,
             do_quality_base,
             do_sequence,
             do_gc_sequence,
             do_quality_sequence,
             do_tile,
             do_adapter,
             do_sequence_length,
             do_nogroup;

  // buffer size to store line 2 of each read statically
  const size_t buffer_size;

  // This will tell me which character to look for to go to the next field
  const char field_separator;

  // This will tell me which character to look for to finish processing a record
  const char line_separator;

  // keep track of reads for which to do kmer and tile count
  const size_t num_reads_for_kmer = 50;
  const size_t num_reads_for_tile = 10;
  size_t next_kmer_read;
  size_t next_tile_read;

  // Whether or not to get bases from buffer when reading quality line
  bool read_from_buffer;

  // Whether or not to write bases to buffer when reading sequence
  bool write_to_buffer;

  // Whether to just ignore per tile sequence quality if tiles are not in the
  // name
  bool tile_ignore; 

  // Get a base from the sequence line
  char base_from_buffer;

  // the static buffer that takes the sequence data
  char *buffer;

  // pointer to current read character. Each file type should figure out how to
  // handle this as the file is processed
  char *cur_char;

  // Number of bases that have overflown the buffer
  size_t leftover_ind;

  /********* TILE PARSING ********/
  // tile value parsed from line 1 of each record
  size_t tile_cur;

  // the number of colons (:) needed to be seen until we know we are in a tile
  size_t tile_split_point;

  // Temp variables to be updated as you pass through the file
  size_t base_ind;  // 0,1,2 or 3
  size_t read_pos;  // which base we are at in the read
  size_t quality_value;  // to convert from ascii to number
  size_t cur_gc_count;  // Number of gc bases in read
  size_t cur_quality;  // Sum of quality values in read
  size_t num_bases_after_n;  // count of k-mers that reset at every N
  size_t cur_kmer;  // 32-mer hash as you pass through the sequence line

  // Temporarily store line 2 out of 4 to know the base to which
  // quality characters are associated
  string leftover_buffer;
  string sequence_to_hash;  // the sequence that will be marked for duplication
  string filename;


  /************ FUNCTIONS TO PROCESS READS AND BASES ***********/
  // gets and puts bases from and to buffer
  inline void put_base_in_buffer();  // puts base in buffer or leftover
  inline void get_base_from_buffer();  // gets base from buffer or leftover

  // on the first tile-processed read, we will try to figure out how tiles
  // should be parsed
  inline void get_tile_split_position();
  inline void get_tile_value();

  inline void process_sequence_base_from_buffer(FastqStats &stats);
  inline void process_sequence_base_from_leftover(FastqStats &stats);
  inline void postprocess_sequence_line(FastqStats &stats);

  inline void process_quality_base_from_buffer(FastqStats &stats);
  inline void process_quality_base_from_leftover(FastqStats &stats);

  inline void postprocess_fastq_record(FastqStats &stats);

  /************ FUNCTIONS TO READ LINES IN DIFFERENT WAYS ***********/
  inline void read_fast_forward_line();  // run this to ignore a line
  inline void skip_separator();  // keep going forward while = separator
  inline void read_tile_line(FastqStats &stats);  // get tile from read name
  inline void read_sequence_line(FastqStats &stats);  // parse sequence
  inline void read_quality_line(FastqStats &stats);  // parse quality

  StreamReader (Config &config,
                const size_t buffer_size,
                const char _field_separator,
                const char _line_separator);

  /************ FUNCTIONS TO IMPLEMENT BASED ON FILE FORMAT  ***********/
  virtual void load() = 0;
  virtual bool operator >> (FastqStats &stats) = 0;
  virtual bool is_eof() = 0; // whether file has ended
  virtual ~StreamReader() = 0;
};

StreamReader::StreamReader(Config &config,
                           const size_t _buffer_size,
                           const char _field_separator,
                           const char _line_separator) : 
  // I have to pass the config skips as const to read them fast
  do_duplication(config.do_duplication),
  do_kmer(config.do_kmer),
  do_n_content(config.do_n_content),
  do_overrepresented(config.do_overrepresented),
  do_quality_base(config.do_quality_base),
  do_sequence(config.do_sequence),
  do_gc_sequence(config.do_gc_sequence),
  do_quality_sequence(config.do_quality_sequence),
  do_tile(config.do_tile),
  do_adapter(config.do_adapter),
  do_sequence_length(config.do_sequence_length),
  do_nogroup (config.nogroup), 

  // Here are the usual stream reader configs
  buffer_size (_buffer_size),
  field_separator (_field_separator),
  line_separator (_line_separator) {

  // Allocates buffer to temporarily store reads
  buffer = new char[buffer_size + 1];
  buffer[buffer_size] = '\0';

  // Tile init
  tile_ignore = !do_tile;  //early ignore tile if asked to skip it
  tile_cur = 0;
  tile_split_point = 0;
  next_tile_read = 0;

  // kmer init
  next_kmer_read = 0;
 
  // Subclasses will use this to deflate if necessary
  filename = config.filename;
}

// Makes sure that any subclass deletes the buffer
StreamReader::~StreamReader () {
  delete buffer;
}

/*******************************************************/
/*************** BUFFER MANAGEMENT *********************/
/*******************************************************/
// puts base either on buffer or leftover
inline void
StreamReader::put_base_in_buffer() {
  base_from_buffer = *cur_char;
  if (write_to_buffer) {
    buffer[read_pos] = base_from_buffer;
  } else {
    if (leftover_ind == leftover_buffer.size())
      leftover_buffer.push_back(base_from_buffer);
    else
      leftover_buffer[leftover_ind] = base_from_buffer;
  } 
}

// Gets base from either buffer or leftover
inline void
StreamReader::get_base_from_buffer() {
  if (read_from_buffer) {
    base_from_buffer = buffer[read_pos];
  } else {
    base_from_buffer = leftover_buffer[leftover_ind];
  }
}

/*******************************************************/
/*************** FAST FOWARD ***************************/
/*******************************************************/

// Keeps going forward while the current character is a separator
inline void
StreamReader::skip_separator() {
  for(; *cur_char == field_separator; ++cur_char) {}
}

// Skips lines that are not relevant
inline void
StreamReader::read_fast_forward_line(){
  for (; *cur_char != field_separator; ++cur_char) {}
}

/*******************************************************/
/*************** TILE PROCESSING ***********************/
/*******************************************************/
// Parse the comment 
inline void 
StreamReader::get_tile_split_position(){
  size_t num_colon = 0;

  // Count colons to know the formatting pattern
  for (; *cur_char != field_separator; ++cur_char) {
    if (*cur_char == ':') ++num_colon;
  }

  // Copied from fastqc
  if (num_colon >= 6) tile_split_point = 4;
  else if (num_colon >=4) tile_split_point = 2;

  // We will not get a tile out of this
  else tile_ignore = true;

}

inline void
StreamReader::get_tile_value() {
  tile_cur = 0;
  size_t num_colon = 0;
  for(; *cur_char != field_separator; ++cur_char) {
    if (*cur_char == ':') ++num_colon;
    if (num_colon == tile_split_point) {
      ++cur_char;

      // parse till next colon or \n
      for(; (*cur_char != ':') && (*cur_char != field_separator); ++cur_char) {
        tile_cur = tile_cur*10 + (*cur_char - '0');
      }
      ++num_colon;
    }
  }
}

// Gets the tile from the sequence name (if applicable)
inline void
StreamReader::read_tile_line(FastqStats &stats){
  // if there is no tile information in the fastq header, fast
  // forward this line
  if (tile_ignore) {
    read_fast_forward_line();
    return;
  }

  // Fast forward if this is not a tile line
  if (stats.num_reads != next_tile_read) {
    read_fast_forward_line();
    return;
  }

  // We haven't parsed the first line to know the split point
  if (tile_split_point == 0)  {
    get_tile_split_position();
  }
  else {
    get_tile_value();
    tile_cur -= stats.kMinTileNumber;
  }
}


/*******************************************************/
/*************** SEQUENCE PROCESSING *******************/
/*******************************************************/

// This is probably the most important function for speed, so it must be really
// optimized at all times
inline void
StreamReader::process_sequence_base_from_buffer(FastqStats &stats) {
  // I will count the Ns even if asked to ignore, as checking ifs take time
  if(base_from_buffer == 'N') {
    stats.n_base_count[read_pos]++;
    num_bases_after_n = 1;  // start over the current kmer
  }

  // ATGC bases
  else {
    // two bit base index
    base_ind = actg_to_2bit(base_from_buffer);

    // increments basic statistic counts
    cur_gc_count += (base_ind & 1);
    stats.base_count[(read_pos << stats.kBitShiftNucleotide) | base_ind]++;

    if (do_adapter) {
      if (stats.num_reads == next_kmer_read) {
        if (read_pos < stats.kKmerMaxBases) {
          cur_kmer = ((cur_kmer << stats.kBitShiftNucleotide) | base_ind);

          // registers k-mer if we've seen at least k nucleotides since the last n
          if (num_bases_after_n == stats.kmer_size) {
            stats.kmer_count[(read_pos << stats.kBitShiftKmer) 
                           | (cur_kmer & stats.kmer_mask)]++;
          }

          else
            num_bases_after_n++;
        }
      }
    }
  }
}

// slower version of process_sequence_base_from_buffer that dynamically
// allocates
inline void
StreamReader::process_sequence_base_from_leftover(FastqStats &stats) {
  if(base_from_buffer == 'N') {

    stats.long_n_base_count[leftover_ind]++;
    num_bases_after_n = 1;  // start over the current kmer
  }

  // ATGC bases
  else {
    // two bit base index
    base_ind = actg_to_2bit(base_from_buffer);

    // increments basic statistic counts
    cur_gc_count += (base_ind & 1);
    stats.long_base_count[(leftover_ind << stats.kBitShiftNucleotide)
                         | base_ind]++;

    // WE WILL NOT DO KMER STATS OUTSIDE OF BUFFER
  }
}

// Gets statistics after reading the entire sequence line
inline void
StreamReader::postprocess_sequence_line(FastqStats &stats) {
  if (do_sequence_length) {
    // read length frequency histogram
    if((read_pos != 0) && (read_pos <= stats.kNumBases))
      stats.read_length_freq[read_pos - 1]++;
    else
      stats.long_read_length_freq[leftover_ind - 1]++;
  }

  // Updates maximum read length if applicable
  if (read_pos > stats.max_read_length)
    stats.max_read_length = read_pos;

  // Registers GC % in the bin truncated to the nearest integer
  if (do_gc_sequence)
    stats.gc_count[round(100 * cur_gc_count / static_cast<double>(read_pos))]++;

}

// Reads the line that has the biological sequence
inline void
StreamReader::read_sequence_line(FastqStats &stats){
  // restart line counters
  read_pos = 0;
  cur_gc_count = 0;
  num_bases_after_n = 1;
  write_to_buffer = true;
  leftover_ind = 0;

  /*********************************************************/
  /********** THIS LOOP MUST BE ALWAYS OPTIMIZED ***********/
  /*********************************************************/
  for (; *cur_char != field_separator; ++cur_char) {
    // puts base either on buffer or leftover 
    put_base_in_buffer();
    // Make sure we have memory space to process new base
    if (!write_to_buffer) {
      if (leftover_ind == stats.num_extra_bases) {
        stats.allocate_new_base(tile_ignore);
      }
    }

    // statistics updated base by base
    // use buffer
    if (write_to_buffer) {
      process_sequence_base_from_buffer(stats);
    }

    // use dynamic allocation
    else {
      process_sequence_base_from_leftover(stats);
    }
    // increase leftover position if not writing to buffer anymore
    if (!write_to_buffer) {
      leftover_ind++;
    }

    // either way increase read position
    ++read_pos;

    // if we reached the buffer size, stop
    if(read_pos == buffer_size) {
      write_to_buffer = false;
    }
  }

  // statistics summarized after the read
  postprocess_sequence_line(stats);
}
/*******************************************************/
/*************** QUALITY PROCESSING ********************/
/*******************************************************/

// Process quality value the fast way from buffer
inline void
StreamReader::process_quality_base_from_buffer(FastqStats &stats) {
  // Average quality in position
  stats.position_quality_count[
        (read_pos << stats.kBitShiftQuality) | quality_value]++;

  // Tile processing
  if (!tile_ignore)
    if ((stats.num_reads == next_tile_read) && tile_cur > 0) {
      stats.tile_position_quality[(read_pos << stats.kBitShiftTile) | tile_cur]
          += quality_value;

  }
}

// Slow version of function above
inline void
StreamReader::process_quality_base_from_leftover(FastqStats &stats) {
  // Average quality in position
  stats.long_position_quality_count[
        (leftover_ind << stats.kBitShiftQuality) | quality_value]++;

  // Tile processing
  if (!tile_ignore) {
    if ((stats.num_reads == next_tile_read) && tile_cur > 0) {
      stats.long_tile_position_quality[(leftover_ind << stats.kBitShiftTile) | tile_cur]
          += quality_value;
    }
  }
}

// Reads the quality line of each base.
inline void
StreamReader::read_quality_line(FastqStats &stats){
  // reset quality counts
  read_pos = 0;
  cur_quality = 0;
  read_from_buffer = true;
  leftover_ind = 0;

  // For quality, we do not look for the separator, but rather for an explicit
  // newline or EOF in case the file does not end with newline or we are getting
  // decompressed strings from a stream
  for (; (*cur_char != line_separator) && !is_eof(); ++cur_char) {
    get_base_from_buffer();

    // Converts quality ascii to zero-based
    quality_value = *cur_char - stats.kBaseQuality;

    // Fast bases from buffer 
    if (read_from_buffer) {
      process_quality_base_from_buffer(stats);
    }

    // Slow bases from dynamic allocation
    else {
      process_quality_base_from_leftover(stats);
    }

    // Sums quality value so we can bin the average at the end
    cur_quality += quality_value;

    if (!read_from_buffer) {
      ++leftover_ind;
    }

    // Flag to start reading and writing outside of buffer
    ++read_pos;
    if (read_pos == buffer_size) {
      read_from_buffer = false;
    }
  }

  // Average quality approximated to the nearest integer. Used to make a
  // histogram in the end of the summary.
  if (do_quality_sequence) {
    stats.quality_count[cur_quality / read_pos]++;  // avg quality histogram
  }
}

/*******************************************************/
/*************** POST LINE PROCESSING ******************/
/*******************************************************/

/*************** THIS IS VERY SLOW ********************/
inline void
StreamReader::postprocess_fastq_record(FastqStats &stats) {

  if (do_overrepresented || do_duplication) {
    // if reads are >75pb, truncate to 50
    if(read_pos <= stats.kDupReadMaxSize || do_nogroup)
      buffer[read_pos] = '\0';
    else
      buffer[stats.kDupReadTruncateSize] = '\0';
    sequence_to_hash = string(buffer);

    // New sequence found 
    if(stats.sequence_count.count(sequence_to_hash) == 0) {
      if (stats.num_unique_seen != stats.kDupUniqueCutoff) {
        stats.sequence_count.insert({{sequence_to_hash, 1}});
        stats.count_at_limit = stats.num_reads;
        ++stats.num_unique_seen;
      }
    } else {
      stats.sequence_count[sequence_to_hash]++;
      if (stats.num_unique_seen < stats.kDupUniqueCutoff)
        stats.count_at_limit = stats.num_reads;
    }
  }

  // counts tile if applicable
  if (!tile_ignore) {
    if (stats.num_reads == next_tile_read) {
      next_tile_read += num_reads_for_tile;
      stats.num_reads_tile++;
      if (tile_cur > 0)
        stats.tile_count[tile_cur]++;
    }
  }

  // I counted kmers here so register that I did so
  if (do_adapter) {
    if (stats.num_reads == next_kmer_read) {
      next_kmer_read += num_reads_for_kmer;
      stats.num_reads_kmer++;
    }
  }
}

/*******************************************************/
/*************** READ FASTQ RECORD *********************/
/*******************************************************/
class FastqReader : public StreamReader {
 private:
  // for uncompressed
  struct stat st;
  char *last;
  void *mmap_data;

 public:
  FastqReader(Config &_config,
              const size_t _buffer_size);

  void load();
  inline bool is_eof();
  inline bool operator >> (FastqStats &stats);
  ~FastqReader();
};

// Set fastq separator as \n
FastqReader::FastqReader (Config &_config,
                          const size_t _buffer_size) :
StreamReader(_config, _buffer_size, '\n', '\n') {
}

inline bool
FastqReader::is_eof () {
  return cur_char == (last + 1);
}

// Load uncompressed fastq through memory map
void
FastqReader::load() {
  // uncompressed fastq = memorymap
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

  // Initialize position pointer
  cur_char = static_cast<char*>(mmap_data);
  last = cur_char + st.st_size - 1;
}

// Parses the particular fastq format
inline bool
FastqReader::operator >> (FastqStats &stats) {
  read_tile_line(stats);
  skip_separator();

  read_sequence_line(stats);
  skip_separator();

  read_fast_forward_line();
  skip_separator();

  read_quality_line(stats);
  skip_separator();

  postprocess_fastq_record(stats);

  // Successful read, increment number in stats
  stats.num_reads++;

  // Returns if file should keep being checked
  return !is_eof();
}

FastqReader::~FastqReader()  {
  munmap(mmap_data, st.st_size);
}
/*******************************************************/
/*************** READ FASTQ GZ RCORD *******************/
/*******************************************************/
class GzFastqReader : public FastqReader {
 private:
  static const size_t chunk_size = 16384;
  char gzbuf[chunk_size];
  gzFile fileobj;

 public:
  GzFastqReader(Config &_config,
                const size_t _buffer_size);
  void load();
  inline bool is_eof();
  inline bool operator >> (FastqStats &stats);
  ~GzFastqReader();
};

// the gz fastq constructor is the same as the fastq
GzFastqReader::GzFastqReader (Config &_config,
                              const size_t _buffer_size) :
FastqReader(_config, _buffer_size) {
}

// Load fastq with zlib
void
GzFastqReader::load() {
  fileobj = gzopen(filename.c_str(), "r");
  cur_char = new char[1];
  ++cur_char;
}

// straightforward
inline bool
GzFastqReader::is_eof() {
  return gzeof(fileobj);
}

GzFastqReader::~GzFastqReader() {
  gzclose_r (fileobj);
}

// Parses fastq gz by reading line by line into the gzbuf
inline bool
GzFastqReader::operator >> (FastqStats &stats) {
  cur_char = gzgets (fileobj, gzbuf, chunk_size);

  // need to check here if we did not hit eof
  if(is_eof())
    return false;

  read_tile_line(stats);
  skip_separator();

  cur_char = gzgets (fileobj, gzbuf, chunk_size);
  read_sequence_line(stats);
  skip_separator();

  cur_char = gzgets (fileobj, gzbuf, chunk_size);
  read_fast_forward_line();
  skip_separator();

  cur_char = gzgets (fileobj, gzbuf, chunk_size);
  read_quality_line(stats);
  skip_separator();

  postprocess_fastq_record(stats);

  // Successful read, increment number in stats
  stats.num_reads++;

  // Returns if file should keep being checked
  return !is_eof();
}

/*******************************************************/
/*************** READ SAM RECORD ***********************/
/*******************************************************/

class SamReader : public StreamReader {
 private:
  // for uncompressed
  struct stat st;
  void *mmap_data;
  char *last;
 public:
  SamReader (Config &_config,
             const size_t _buffer_size);
  void load();
  inline bool is_eof();
  inline bool operator >> (FastqStats &stats);
  ~SamReader();
};

// set sam separator as tab
SamReader::SamReader (Config &_config,
                      const size_t _buffer_size) :
StreamReader(_config, _buffer_size, '\t', '\n') {

}

void
SamReader::load() {
  // uncompressed fastq = memorymap
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

  // Initialize position pointer
  cur_char = static_cast<char*>(mmap_data);
  last = cur_char + st.st_size - 1;

  // Skip sam header
  while (*cur_char == '@') {
    for(; *cur_char != line_separator; ++cur_char);
    ++cur_char;
  }
}

inline bool
SamReader::is_eof() {
  return (cur_char == last + 1);
}

inline bool
SamReader::operator >> (FastqStats &stats) {
  read_tile_line(stats);
  skip_separator();
  for (size_t i = 0; i < 8; ++i) {
    read_fast_forward_line();
    skip_separator();
  }
  read_sequence_line(stats);
  read_quality_line(stats);

  // skips all tags after quality until newline
  while (*cur_char != line_separator) {
    read_fast_forward_line(); 
    skip_separator();
  }

  // skip \n
  ++cur_char;
  postprocess_fastq_record(stats);
  stats.num_reads++;

  // Returns if file should keep being checked
  return !is_eof();

}

SamReader::~SamReader() {
  munmap(mmap_data, st.st_size);
}
/*******************************************************/
/*************** READ BAM RECORD ***********************/
/*******************************************************/

class BamReader : public StreamReader {
 private:
  htsFile *hts;
  bam_hdr_t *hdr;
  bam1_t *b;
  int rd_ret, fmt_ret;
  char *last;
 public:
  BamReader (Config &_config, 
             const size_t _buffer_size);
  void load();
  bool is_eof();
  bool operator >> (FastqStats &stats);
  ~BamReader();
};

// set sam separator as tab
BamReader::BamReader (Config &_config,
                      const size_t _buffer_size) : 
StreamReader(_config, _buffer_size, '\t', '\n') {
  rd_ret = 0;

}

void
BamReader::load() {
  if (!(hts = hts_open (filename.c_str(), "r")))
    throw runtime_error ("cannot load bam file : " + filename);

  if (!(hdr = sam_hdr_read(hts)))
     throw runtime_error("failed to read header from file: " + filename);

  if (!(b = bam_init1()))
     throw runtime_error("failed to read record from file: " + filename);

}

// We will check eof on the >> operator
bool 
BamReader::is_eof (){
  return cur_char == last - 1;
}

inline bool
BamReader::operator >> (FastqStats &stats) {
  if ((rd_ret = sam_read1(hts, hdr, b)) >= 0) {
    fmt_ret = 0;
    if ((fmt_ret = sam_format1(hdr, b, &hts->line)) > 0) {
      // define char* values for processing lines char by char
      cur_char = hts->line.s;
      last = cur_char + strlen (hts->line.s) - 1;

      // Now read it as regular sam
      read_tile_line(stats);
      skip_separator();
      for (size_t i = 0; i < 8; ++i) {
        read_fast_forward_line();
        skip_separator();
      }

      read_sequence_line(stats);
      read_quality_line(stats);

      postprocess_fastq_record(stats);
      stats.num_reads++;
      return true;
    } else {
      throw runtime_error ("failed reading record from : " + filename);
    }

    return false;
  }

  // If I could not read another line it means it's eof
  return false;
}

BamReader::~BamReader() {
  if (hdr) {
    bam_hdr_destroy(hdr);
    hdr = 0;
  }

  if (b) {
    bam_destroy1(b);
    b = 0;
  }

  if (hts) {
    hts_close(hts);
    hts = 0;
  }
}

/*******************************************************/
/*************** HTML FACTORY **************************/
/*******************************************************/
class HTMLFactory {
 private:
  void replace_placeholder_with_data (const string &placeholder, 
                                      const string &data);

  void comment_if_skipped (string ph_begin, string ph_end, bool skip);
  // First thing to do: file details
  void put_file_details (const Config &config);

  // Comment out analyses that were skipped
  void comment_out (const Config &config);

  // Functions to replace pass warn fail
  void put_pass_warn_fail (const FastqStats &stats);

  // Function to replace template placeholders with data
  void make_basic_statistics(const FastqStats &stats,
                             Config &config);

  void make_position_quality_data(const FastqStats &stats, 
                                  const Config &config);

  void make_tile_quality_data(const FastqStats &stats,
                              const Config &config);

  void make_sequence_quality_data(const FastqStats &stats,
                                  const Config &config);
  void make_base_sequence_content_data(const FastqStats &stats,
                                       const Config &config);

  void make_sequence_gc_content_data(const FastqStats &stats,
                                     const Config &config);

  void make_base_n_content_data(const FastqStats &stats,
                                const Config &config);

  void make_sequence_length_data(const FastqStats &stats,
                                 const Config &config);

  void make_sequence_duplication_data(const FastqStats &stats,
                                      const Config &config);

  void make_overrepresented_sequences_data(const FastqStats &stats,
                                           const Config &config);

  void make_adapter_content_data(FastqStats &stats,
                                 Config &config);

 public:
  string sourcecode;
  explicit HTMLFactory (string filepath);
  void write(FastqStats &stats, Config &config);
};

HTMLFactory::HTMLFactory (string filepath) {
  sourcecode = "";
  ifstream in(filepath);
  if(!in)
    throw runtime_error ("HTML layout not found: " + filepath);
  string line;

  // pass the whole source code template to a string
  while (getline(in, line))
    sourcecode += line + "\n";
}

void
HTMLFactory::replace_placeholder_with_data (const string &placeholder, 
                                            const string &data) {
  auto pos = sourcecode.find(placeholder);

  // Placeholder not found
  if (pos == string::npos)
    throw runtime_error ("placeholder not found: " + placeholder);

  // at least one placeholder found
  while (pos != string::npos) {
    sourcecode.replace(pos, placeholder.size(), data);
    pos = sourcecode.find(placeholder, pos + 1);
  }
}

// Comments out html pieces if analyses were skipped
void
HTMLFactory::comment_if_skipped (string ph_begin,
                                 string ph_end,
                                 bool done) {
  // put html comments if analysis was skipped
  if (!done) {
    replace_placeholder_with_data(ph_begin, "<!--");
    replace_placeholder_with_data(ph_end, "-->");
  }

  // otherwise delete placeholder
  else {
    replace_placeholder_with_data(ph_begin, "");
    replace_placeholder_with_data(ph_end, "");
  }
}

void
HTMLFactory::comment_out (const Config &config) {
  comment_if_skipped ("{{skipbasesequence_s}}",
                      "{{skipbasesequence_e}}",
                      config.do_quality_base);
  comment_if_skipped ("{{skiptile_s}}",
                      "{{skiptile_e}}",
                      config.do_tile);

  comment_if_skipped ("{{skipsequencequal_s}}",
                      "{{skipsequencequal_e}}",
                      config.do_quality_sequence);

  comment_if_skipped ("{{skipbasecontent_s}}",
                      "{{skipbasecontent_e}}",
                      config.do_sequence);

  comment_if_skipped ("{{skipsequencegc_s}}",
                      "{{skipsequencegc_e}}",
                      config.do_gc_sequence);

  comment_if_skipped ("{{skipbasencontent_s}}",
                      "{{skipbasencontent_e}}",
                      config.do_n_content);

  comment_if_skipped ("{{skipseqlength_s}}",
                      "{{skipseqlength_e}}",
                      config.do_n_content);

  comment_if_skipped ("{{skipseqdup_s}}",
                      "{{skipseqdup_e}}",
                      config.do_duplication);

  comment_if_skipped ("{{skipoverrep_s}}",
                      "{{skipoverrep_e}}",
                      config.do_overrepresented);

  comment_if_skipped ("{{skipadapter_s}}",
                      "{{skipadapter_e}}",
                      config.do_adapter);
}

void
HTMLFactory::put_file_details (const Config &config) {
  // Put filename in filename placeholder
  replace_placeholder_with_data("{{filename}}", 
                                strip_path(config.filename));

  // Put date on date placeholder
  auto tmp = system_clock::to_time_t(system_clock::now());
  string time_fmt = string(ctime(&tmp));
  replace_placeholder_with_data("{{date}}", time_fmt);
}

void
HTMLFactory::put_pass_warn_fail (const FastqStats &stats) {
  replace_placeholder_with_data ("{{passbasic}}",
                                 stats.pass_basic_statistics);
  replace_placeholder_with_data ("{{passbasesequence}}",
                                 stats.pass_per_base_sequence_quality);
  replace_placeholder_with_data ("{{passpertile}}",
                                 stats.pass_per_tile_sequence_quality);
  replace_placeholder_with_data ("{{passpersequencequal}}",
                                 stats.pass_per_sequence_quality_scores);
  replace_placeholder_with_data ("{{passperbasecontent}}",
                                 stats.pass_per_base_sequence_content);
  replace_placeholder_with_data ("{{passpersequencegc}}",
                                 stats.pass_per_sequence_gc_content);
  replace_placeholder_with_data ("{{passperbasencontent}}",
                                 stats.pass_per_base_n_content);
  replace_placeholder_with_data ("{{passseqlength}}",
                                 stats.pass_sequence_length_distribution);
  replace_placeholder_with_data ("{{passseqdup}}",
                                 stats.pass_duplicate_sequences);
  replace_placeholder_with_data ("{{passoverrep}}",
                                 stats.pass_overrepresented_sequences);
  replace_placeholder_with_data ("{{passadapter}}", 
                                 stats.pass_adapter_content);
}
void 
HTMLFactory::make_basic_statistics(const FastqStats &stats,
                                   Config &config) {
  string placeholder = "{{BASICSTATSDATA}}";
  ostringstream data;
  data << "<table><thead><tr><th>Measure</th><th>Value</th></tr></thead><tbody>";
  data << "<tr><td>Filename</td><td>" << strip_path(config.filename)
       << "</td></tr>";
  data << "<tr><td>Total Sequences</td><td>" << stats.num_reads << "</td></tr>";
  data << "<tr><td>Sequences Flagged As Poor Quality</td><td>" 
       << stats.num_poor << "</td></tr>";
  data << "<tr><td>Sequence length</td><td>" ;
  if(stats.min_read_length != stats.max_read_length)
    data << stats.min_read_length << " - " << stats.max_read_length;
  else
    data << stats.max_read_length;
  data << "</td></tr>";

  if (config.do_sequence)
    data << "<tr><td>%GC:</td><td>" << stats.avg_gc << "</td></tr>";
  data << "</tbody></table>";

  replace_placeholder_with_data (placeholder, data.str());
}

void 
HTMLFactory::make_position_quality_data (const FastqStats &stats,
                                         const Config &config) {
  ostringstream data;
  const string placeholder = "{{SEQBASEQUALITYDATA}}";
  if (config.do_quality_base) {
    size_t cur_median;
    for (size_t i = 0; i < stats.max_read_length; ++i) {
      data << "{y : [";

     if (i < stats.kNumBases) {
       cur_median = stats.median[i];
       data << stats.ldecile[i] << ", "
            << stats.lquartile[i] << ", "
            << stats.median[i] << ", "
            << stats.uquartile[i] << ", "
            << stats.udecile[i] << "], ";
     } else {
       cur_median = stats.long_median[i - stats.kNumBases];
       data << stats.long_ldecile[i - stats.kNumBases] << ", "
            << stats.long_lquartile[i - stats.kNumBases] << ", "
            << stats.long_median[i - stats.kNumBases] << ", "
            << stats.long_uquartile[i - stats.kNumBases] << ", "
            << stats.long_udecile[i - stats.kNumBases] << "], ";
     }
     data << "type : 'box', name : ' " << i << "', ";
     data << "marker : {color : '";
     if(cur_median > 30)
       data << "green";
     else if (cur_median > 20)
       data << "yellow";
     else
       data << "red";
     data << "'}}";
     if (i < stats.max_read_length - 1)
       data << ", ";

    }
  }
  replace_placeholder_with_data (placeholder, data.str());
}

void
HTMLFactory::make_tile_quality_data (const FastqStats &stats,
                                     const Config &config) {
  ostringstream data;
  const string placeholder = "{{TILEQUALITYDATA}}";

  if (config.do_tile) {
    // X: position
    data << "{x : [";
    for (size_t i = 0; i < stats.max_read_length; ++i) {
      data << i+1;
      if (i < stats.max_read_length - 1)
        data << ",";
    }

    // Y : Tile
    data << "], y: [";
    bool first_seen = false;
    for (size_t i = 0; i < stats.kNumMaxTiles; ++i) {
      if(stats.tile_count[i] > 0) {
        if(!first_seen)
          first_seen = true;
        else
          data << ",";
        data << i;
      }
    }

    // Z: quality z score
    data << "], z: [";
    first_seen = false;
    for (size_t i = 0; i < stats.kNumMaxTiles; ++i) {
      if(stats.tile_count[i] > 0) {
        if(!first_seen)
          first_seen = true;
        else
          data << ", ";

        // datart new array with all counts
        data << "[";
        for (size_t j = 0; j < stats.max_read_length; ++j) {
          if (j < stats.kNumBases)
            data << stats.tile_position_quality[(j << stats.kBitShiftTile) | i];
          else
            data << stats.long_tile_position_quality[
              ((j - stats.kNumBases) << stats.kBitShiftTile) | i];

          if (j < stats.max_read_length - 1)
            data << ",";
        }
        data << "]";
      }
    }
    data << "]";
    data << ", type : 'heatmap' }";
  }
  replace_placeholder_with_data (placeholder, data.str());
}

void 
HTMLFactory::make_sequence_quality_data (const FastqStats &stats,
                                         const Config &config) {
  ostringstream data;
  const string placeholder = "{{SEQQUALITYDATA}}";

  if (config.do_quality_sequence) {
    // X values : avg quality phred scores
    data << "{x : [";
    for (size_t i = 0; i < 41; ++i) {
      data << i + stats.kBaseQuality;
     if (i < 40)
       data << ", ";
    }

    // Y values: frequency with which they were seen
    data << "], y : [";
    for (size_t i = 0; i < 41; ++i) {
      data << stats.quality_count[i];
     if (i < 40)
       data << ", ";
    }
    data << "], type: 'line', line : {color : 'red'}, "
         << "name : 'Sequence quality distribution'}";
  }
  replace_placeholder_with_data (placeholder, data.str());
}

void 
HTMLFactory::make_base_sequence_content_data (const FastqStats &stats,
                                              const Config &config) {
  ostringstream data;
  const string placeholder = "{{BASESEQCONTENTDATA}}";

  if (config.do_sequence) {
    // ATGC
    for (size_t base = 0; base < stats.kNumNucleotides; ++base){
      // start line
      data << "{";

      // X values : base position
      data << "x : [";
      for (size_t i = 0; i < stats.max_read_length; ++i) {
        data << i+1;
        if (i < stats.max_read_length - 1)
          data << ", ";
      }

      // Y values: frequency with which they were seen
      data << "], y : [";
      for (size_t i = 0; i < stats.max_read_length; ++i) {
        if (base == 0) {
          if (i < stats.kNumBases) data << stats.a_pct[i];
          else data << stats.long_a_pct[i - stats.kNumBases];
        }
        if (base == 1) {
          if (i < stats.kNumBases) data << stats.c_pct[i];
          else data << stats.long_c_pct[i - stats.kNumBases];
        }
        if (base == 2) {
          if (i < stats.kNumBases) data << stats.t_pct[i];
          else data << stats.long_t_pct[i - stats.kNumBases];
        }
        if (base == 3) {
          if (i < stats.kNumBases) data << stats.g_pct[i];
          else data << stats.long_g_pct[i - stats.kNumBases];
        }
        if (i < stats.max_read_length - 1)
          data << ", ";
      }
      data << "], mode : 'lines', name : '" + size_t_to_seq (base, 1) + "', ";

      // color
      data << "line :{ color : '";
      if(base == 0) 
        data << "green";
      else if (base == 1)
        data << "blue";
      else if (base == 2)
        data << "red";
      else 
        data << "black";
      data << "'}";
      // end color

      // end line
      data << "}";
      if (base < stats.kNumNucleotides - 1)
        data << ", ";
    }
  }
  replace_placeholder_with_data (placeholder, data.str());
}

void 
HTMLFactory::make_sequence_gc_content_data (const FastqStats &stats,
                                            const Config &config) {
  ostringstream data;
  const string placeholder = "{{SEQGCCONTENTDATA}}";

  if (config.do_gc_sequence) {
    // Actual count
    data << "{x : [";
    for (size_t i = 0; i < 101; ++i) {
      data << i + 1;
      if (i < 101)
        data << ", ";
    }

    // Y values: frequency with which they were seen
    data << "], y : [";
    for (size_t i = 0; i < 101; ++i) {
      data << stats.gc_count[i];
     if (i < 101)
       data << ", ";
    }
    data << "], type: 'line', line : {color : 'red'}, name : 'GC distribution'}";

    // Theoretical count
    data << ", {x : [";
    for (size_t i = 0; i < 101; ++i) {
      data << i + 1;
      if (i < 101)
        data << ", ";
    }

    // Y values: frequency with which they were seen
    data << "], y : [";
    for (size_t i = 0; i < 101; ++i) {
      data << stats.theoretical_gc_count[i];
     if (i < 101)
       data << ", ";
    }
    data << "], type: 'line', line : {color : 'blue'},"
         << "name : 'Theoretical distribution'}";
  }
  replace_placeholder_with_data (placeholder, data.str());
}

void 
HTMLFactory::make_base_n_content_data (const FastqStats &stats,
                                       const Config &config) {
  ostringstream data;
  const string placeholder = "{{BASENCONTENTDATA}}";

  if (config.do_n_content) {
    // base position
    data << "{x : [";
    for (size_t i = 0; i < stats.max_read_length; ++i) {
      data << i + 1;
      if (i < stats.max_read_length - 1)
        data << ", ";
    }

    // Y values: frequency with which they were seen
    data << "], y : [";
    for (size_t i = 0; i < stats.max_read_length; ++i) {
      if (i < stats.kNumBases)
        data << stats.n_pct[i];
      else
        data << stats.long_n_pct[i - stats.kNumBases];

      if (i < stats.max_read_length - 1)
        data << ", ";
    }
    data << "], type: 'line', line : {color : 'red'}, "
         << "name : 'Fraction of N reads per base'}";
  }
  replace_placeholder_with_data (placeholder, data.str());
}

void 
HTMLFactory::make_sequence_length_data (const FastqStats &stats,
                                        const Config &config) {
  ostringstream data;
  const string placeholder = "{{SEQLENDATA}}";

  if (config.do_sequence_length) {
    // X values : avg quality phred scores
    data << "{x : [";
    bool first_seen = false;
    size_t val;
    for (size_t i = 0; i < stats.max_read_length; ++i) {
      val = 0 ;
      if (i < stats.kNumBases){
        if (stats.read_length_freq[i] > 0)
          val = stats.read_length_freq[i];
      } else {
        if (stats.long_read_length_freq[i - stats.kNumBases] > 0)
          val = stats.long_read_length_freq[i - stats.kNumBases];
      }

      if (val > 0) {
        if (first_seen)
          data << ",";
        data << "\"" << i+1 << " bp\"";
        first_seen = true;
      }
    }

    // Y values: frequency with which they were seen
    data << "], y : [";
    first_seen = false;
    for (size_t i = 0; i < stats.max_read_length; ++i) {
      val = 0;
      if (i < stats.kNumBases){
        if (stats.read_length_freq[i] > 0)
          val = stats.read_length_freq[i];
      } else {
        if (stats.long_read_length_freq[i - stats.kNumBases] > 0)
          val = stats.long_read_length_freq[i - stats.kNumBases];
      }

      if (val > 0) {
        if (first_seen)
          data << ",";
        data << val;
        first_seen = true;
      }
    }

    // Put the sequence value in the text
    data << "], text : [";
    first_seen = false;
    for (size_t i = 0; i < stats.max_read_length; ++i) {
      val = 0;
      if (i < stats.kNumBases){

        if (stats.read_length_freq[i] > 0)
          val = stats.read_length_freq[i];
      } else {
        if (stats.long_read_length_freq[i - stats.kNumBases] > 0)
          val = stats.long_read_length_freq[i - stats.kNumBases];
      }

      if (val > 0) {
        if (first_seen)
          data << ",";
        data << i+1;
        first_seen = true;
      }
    }


    data << "], type: 'bar', marker : {color : 'rgba(55,128,191,1.0)', line : {width : 2}}, "
         << "name : 'Sequence length distribution'}";
  }
  replace_placeholder_with_data (placeholder, data.str());
}


void
HTMLFactory::make_sequence_duplication_data (const FastqStats &stats,
                                             const Config &config) {
  ostringstream data;
  const string placeholder = "{{SEQDUPDATA}}";

  if (config.do_duplication) {
    // non-deduplicated
    data << "{x : [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]";

    // total percentage in each bin
    data << ", y : [";
    for (size_t i = 0; i < 16; ++i) {
      data << stats.percentage_total[i];

      if (i < 15)
        data << ", ";
    }
    data << "], type: 'line', line : {color : 'blue'}, "
         << "name : 'total sequences'}";

    // deduplicated
    data << ", {x : [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]";


    // total percentage in deduplicated
    data << ", y : [";
    for (size_t i = 0; i < 16; ++i) {
      data << stats.percentage_deduplicated[i];

      if (i < 15)
        data << ", ";
    }
    data << "], type: 'line', line : {color : 'red'}, "
         << "name : 'total sequences'}";
  }
  replace_placeholder_with_data (placeholder, data.str());
}


void 
HTMLFactory::make_overrepresented_sequences_data(const FastqStats &stats,
                                                 const Config &config) {
  string placeholder = "{{OVERREPSEQDATA}}";
  ostringstream data;

  if (config.do_overrepresented) {
    // Header
    data << "<table><thead><tr>";
    data << "<th>Sequence</th>";
    data << "<th>Count</th>";
    data << "<th>Percentage</th>";
    data << "<th>Possible Source</th>";
    data << "</tr></thead><tbody>";

    // All overrep sequences
    for (auto v : stats.overrep_sequences) {
      data << "<tr><td>" << v.first << "</td>";
      data << "<td>" << v.second << "</td>";
      data << "<td>" << 100.0 * v.second / stats.num_reads << "</td>";
      data << "<td>" << config.get_matching_contaminant(v.first)
          << "</td>";

      data << "</tr>";
    }
    data << "</tbody></table>";
  }
  replace_placeholder_with_data (placeholder, data.str());
}

void
HTMLFactory::make_adapter_content_data (FastqStats &stats,
                                        Config &config) {
  string placeholder = "{{ADAPTERDATA}}";
  ostringstream data;

  if (config.do_adapter) {
    // Number of bases to make adapter content
    size_t num_bases =  min(stats.kNumBases, stats.kKmerMaxBases);
    bool seen_first = false;

    size_t jj = 0;
    for (auto v : config.adapters) {
      if (!seen_first)
        seen_first = true;
      else
        data << ",";
      data << "{";

      // X values : read position
      data << "x : [";
      for (size_t i = 0; i < num_bases; ++i) {
        if (stats.cumulative_read_length_freq[i] > 0){
          data << i+1;
          if (i < num_bases - 1)
            data << ",";
        }
      }
      data << "]";

      // Y values : cumulative adapter frequency
      data << ", y : [";
      for (size_t i = 0; i < num_bases; ++i) {
        if (stats.cumulative_read_length_freq[i] > 0){
          data << stats.kmer_by_base[i][jj];
          if (i < num_bases - 1)
            data << ",";
        }
      }

      data << "]";
      data << ", type : 'line', ";
      data << "name : '" << v.first << "'}";
      ++jj;
    }
  }
  replace_placeholder_with_data (placeholder, data.str());
}

void
HTMLFactory::write (FastqStats &stats, Config &config) {
  // Filename and date
  put_file_details(config);

  // Put html comments around things to skip
  comment_out(config);

  // Put colors in pass, warn and fail
  put_pass_warn_fail (stats);

  // Put data on tables and plotly javascripts
  make_basic_statistics (stats, config);
  make_position_quality_data (stats, config);
  make_tile_quality_data (stats, config);
  make_sequence_quality_data (stats, config);
  make_base_sequence_content_data (stats, config);
  make_sequence_gc_content_data (stats, config);
  make_base_n_content_data (stats, config);
  make_sequence_length_data (stats, config);
  make_sequence_duplication_data (stats, config);
  make_overrepresented_sequences_data (stats, config);
  make_adapter_content_data (stats, config);
}

/******************************************************
 ********************* MAIN ***************************
 ******************************************************/

// Read any file type until the end and logs progress
template <typename T> void
read_stream_into_stats(T &in,
                       FastqStats &stats,
                       Config &config) {
  in.load();

  // Fast way to report # of reads without modulo arithmetics
  const size_t num_reads_to_log = 1000000;
  size_t next_read = num_reads_to_log;

  // Read record by record
  while (in >> stats) {
    if (!config.quiet) {
      if (stats.num_reads == next_read) {
        log_process ("Processed " +
                    to_string(stats.num_reads / num_reads_to_log) +
                    "M reads");
        next_read += num_reads_to_log;
      }
    }
  }

  // if I could not get tile information from read names, I need to tell this to
  // config so it does not output tile data on the summary or html
  if (in.tile_ignore)
    config.do_tile = false;
}

int main(int argc, const char **argv) {
  /****************** COMMAND LINE OPTIONS ********************/
  bool help = false;
  bool version = false;
  Config config;
  string forced_file_format;
  string outdir = "";

  OptionParser opt_parse(strip_path(argv[0]),
                         "A high throughput sequence QC analysis tool",
                         "seqfile1 seqfile2 ... seqfileN");

  opt_parse.add_opt("-help", 'h', "print this help file adn exit", 
                     false, help);

  opt_parse.add_opt("-version", 'v', "print the version of the program and exit", 
                     false, version);

  opt_parse.add_opt("-outdir", 'o', "Create all output files in the specified "
                    "output directory. If not provided, files will be created "
                    "in the same directory as the input file.",
                    false, outdir);

  opt_parse.add_opt("-casava", 'C',
                    "Files come from raw casava output (currently ignored)",
                    false, config.casava);

  opt_parse.add_opt("-nano", 'n',
                    "Files come from fast5 nanopore sequences",
                    false, config.nanopore);

  opt_parse.add_opt("-nofilter", 'F',
                    "If running with --casava do not remove poor quality "
                    "sequences (currently ignored)",
                    false, config.nofilter);

  opt_parse.add_opt("-noextract", 'e',
                    "If running with --casava do not remove poor quality "
                    " sequences (currently ignored)",
                    false, config.casava);

  opt_parse.add_opt("-nogroup", 'g',
                    "Disable grouping of bases for reads >50bp",
                    false, config.nogroup);

  opt_parse.add_opt("-format", 'f',
                    "Force file format",
                    false, forced_file_format);

  opt_parse.add_opt("-threads", 't',
                    "Specifies number of simultaneous files ",
                    false, config.threads );

  opt_parse.add_opt("-contaminants", 'c',
                    "Non-default filer with a list of contaminants",
                    false, config.contaminants_file);

  opt_parse.add_opt("-adapters", 'a',
                    "Non-default file with a list of adapters",
                    false, config.contaminants_file);

  opt_parse.add_opt("-limits", 'l',
                    "Non-default file with limits and warn/fail criteria",
                    false, config.contaminants_file);

  opt_parse.add_opt("-kmer", 'k',
                    "k-mer size (default = 7, max = 10)", false, 
                    config.kmer_size);


  opt_parse.add_opt("-quiet", 'q', "print more run info", false, config.quiet);
  opt_parse.add_opt("-dir", 'd', "directory in which to create temp files", 
                     false, config.quiet);

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

  if (leftover_args.size() == 0) {
    cerr << opt_parse.help_message() << endl;
    return EXIT_SUCCESS;
  }

  for (auto file : leftover_args) {
    auto file_start = system_clock::now();
    config.filename = file;

    // if format was not provided, we have to guess it by the filename
    if (!forced_file_format.empty())
      config.format = forced_file_format;
    else
      config.format = "";

    /****************** BEGIN PROCESSING CONFIG ******************/
    config.setup();  // define file type, read limits, adapters, contaminants

    /****************** END PROCESSING CONFIG *******************/
    if (!config.quiet)
      log_process("Started reading file " +config.filename);

    // Allocates vectors to summarize data
    FastqStats stats = FastqStats();

    // Initializes a reader given the file format
    if (config.format == "sam") {
      log_process("reading file as sam format");
      SamReader in(config, stats.kNumBases);
      read_stream_into_stats(in,stats,config);
    }
    else if (config.format == "bam") {
      log_process("reading file as bam format");
      BamReader in (config, stats.kNumBases);
      read_stream_into_stats(in,stats,config);
    }

    else if (config.compressed) {
      log_process("reading file as gzipped fastq format");
      GzFastqReader in (config, stats.kNumBases);
      read_stream_into_stats(in,stats,config);
    }
    else {
      log_process("reading file as uncompressed fastq format");
      FastqReader in (config, stats.kNumBases);
      read_stream_into_stats(in,stats,config);
    }

    if (!config.quiet)
      log_process("Finished reading file");

    if (!config.quiet)
      log_process("Summarizing data");

    // This function has to be called before writing to output. This is where we
    // calculate all the summary statistics that will be written to output. 
    stats.summarize(config);
    /************************ WRITE TO OUTPUT *****************************/
    string outfile = file;
    if (!outdir.empty())
      outfile = outdir + "/" + strip_path(file);
    outfile += "_qc_data.txt";

    if (!config.quiet)
      log_process("Writing text data to " + outfile);

    // define output
    ofstream of (outfile);

    // Write
    stats.write(of, config);
    of.close();

    /************************ WRITE TO HTML *****************************/
    // Take the html template and populate it with stats data. 
    // Use config to know which parts not to write
    HTMLFactory factory ("Configuration/template.html");
    factory.write (stats, config);

    // Decide html filename based on input
    string htmlfile = file;
    if (!outdir.empty())
      htmlfile = outdir + "/" + strip_path(file);
    htmlfile += "_report.html";

    if (!config.quiet)
      log_process("Writing HTML report to " + htmlfile);

    // Start writing
    ofstream html (htmlfile);
    html << factory.sourcecode;
    html.close();

    /************** TIME SUMMARY *********************************/
    if (!config.quiet) 
      cerr << "Elapsed time for file " << file << ": " <<
      std::chrono::duration_cast<std::chrono::seconds>(
          system_clock::now() - file_start
      ).count() 
      << "s\n";
  }

  return EXIT_SUCCESS;
}


