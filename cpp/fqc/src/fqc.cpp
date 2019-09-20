/* fqc: quality control for fastq files
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
using std::min;

/*************************************************************
 ******************** AUX FUNCTIONS **************************
 *************************************************************/

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
size_t
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

/******************* IMPLEMENTATION OF FASTQC FUNCTIONS **********************/
// FastQC extrapolation of counts to the full file size
double get_corrected_count (size_t count_at_limit,
                            size_t num_reads,
                            size_t dup_level,
                            size_t num_obs) {

  // See if we can bail out early
  if (count_at_limit == num_reads)
    return num_obs;

  // If there aren't enough sequences left to hide another sequence with this count then
  // we can also skip the calculation
  if (num_reads - num_obs < count_at_limit) 
    return num_obs;

  // If not then we need to see what the likelihood is that we had 
  // another sequence with this number of observations which we would 
  // have missed. We'll start by working out the probability of NOT seeing a 
  // sequence with this duplication level within the first count_at_limit 
  // sequences of num_obs.  This is easier than calculating
  // the probability of seeing it.
  double pNotSeeingAtLimit = 1.0;

  // To save doing long calculations which are never going to produce anything meaningful
  // we'll set a limit to our p-value calculation.  This is the probability below which we
  // won't increase our count by 0.01 of an observation.  Once we're below this we stop caring
  // about the corrected value since it's going to be so close to the observed value that
  // we can just return that instead.
  double limitOfCaring = 1.0 - (num_obs/(num_obs + 0.01));
  for (size_t i = 0; i < count_at_limit; ++i) {
    pNotSeeingAtLimit *= static_cast<double>((num_reads-i)-dup_level) /
                         static_cast<double>(num_reads-i);

    if (pNotSeeingAtLimit < limitOfCaring) {
      pNotSeeingAtLimit = 0;
      break;
    }
  }

  // Now we can assume that the number we observed can be 
  // scaled up by this proportion
  return num_obs/(1 - pNotSeeingAtLimit);
}

// Function to calculate the deviation of a histogram with 100 bins from a
// theoretical normal distribution with same mean and standard deviation
double 
sum_deviation_from_normal (const array <size_t, 101> &gc_content) {
  double mode = 0.0, freq = 0.0, num_reads = 0.0;

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
  array <double, 101> theoretical;
  for (size_t i = 0; i < 101; ++i) {
    z = i - mode;
    theoretical[i] = exp(- (z*z)/ (2 * stdev *stdev));
    theoretical_sum += theoretical[i];
  }

  // Normalize theoretical so it sums to the total of reads 
  for (size_t i = 0; i < 101; ++i)
    theoretical[i] = theoretical[i] * num_reads / theoretical_sum;

  for (size_t i = 0; i < 101; ++i)
    ans += fabs(gc_content[i] - theoretical[i]);
  
  // Fractional deviation
  return ans / num_reads;
}

/*************************************************************
 ******************** FASTQ STATS ****************************
 *************************************************************/

struct FastqStats {
 public:
  /************************************************************
   *************** UNIVERSAL CONSTANTS ************************
   ************************************************************/

  // number of bases for which time optimization will be performed at the cost
  // of memory
  static const size_t kNumBases = 60;

  // Value to subtract quality characters to get the actual quality value
  static const size_t kBaseQuality = 33;  // The ascii for the lowest quality

  // Smallest power of two that comprises all possible Illumina quality values.
  // Illumina gives qualities from 0 to 40, therefore we set it as 64. Power of
  // is to avoid double pointer jumps and to get indices with bit shifts.
  static const size_t kNumQualityValues = 64;

  // How many possible nucleotides (must be power of 2!)
  static const size_t kNumNucleotides = 4;  // A = 00,C = 01,T = 10,G = 11

  // maximum tile value
  static const size_t kNumMaxTiles = 65536;

  // threshold for a sequence to be considered  poor quality
  static const size_t kPoorQualityThreshold = 20;

  /************* DUPLICATION ESTIMATES *************/
  // Number of unique sequences to see before stopping counting sequences
  static const size_t kDupUniqueCutoff = 1e5;

  // Maximum read length to store the entire read in memory
  static const size_t kDupReadMaxSize = 75;

  // Prefix size to cut if read length exceeds the value above
  static const size_t kDupReadTruncateSize = 50;

  /************ KMER **********/
  const size_t kKmerMaxBases = 500;  // fastqc's maximum number of bases
 
  /*********** SINGLE NUMBERS FROM THE ENTIRE FASTQ ****************/
  size_t kBitShiftNucleotide;  // log 2 of value above
  size_t kBitShiftTile;
  size_t kBitShiftQuality;  // log 2 of value above

  // Number of unique sequences seen thus far
  size_t num_unique_seen;

  // How many reads were processed before num_unique_seen = kDupUniqueCutoff
  size_t count_at_limit;

  size_t total_bases;  // sum of all bases in all reads
  size_t avg_read_length;  // average of all read lengths
  size_t num_reads;  // total number of lines read
  size_t max_read_length;  // total number of lines read
  size_t num_poor;  // reads whose average quality was <= poor
  size_t num_extra_bases;  // number of bases outside of buffer
  // Kmer size given as input
  size_t kmer_size;

  // mask to get only the first 2*k bits of the sliding window
  size_t kmer_mask;

  double avg_gc;  // (sum of g bases + c bases) / (num_reads)
  double total_deduplicated_pct;  // number of reads left if deduplicated

  /************ OVERREPRESENTATION ESTIMTES **********/
  // fraction of the number of slow reads a sequence needs to be seen to be
  // considered a candiate for overrepresentation
  constexpr static double kOverrepMinFrac = 0.001;


  /*********************************************************
   *********** METRICS COLLECTED DURING IO *****************
   *********************************************************/
  /*********** PER BASE METRICS ****************/

  // counts the number of bases in every read position
  array<size_t, kNumNucleotides * kNumBases> base_count;  // ATGC
  array<size_t, kNumNucleotides * kNumBases> n_base_count; // N

  // Sum of base qualities in every read position
  array<size_t, kNumNucleotides* kNumBases> base_quality;  // ATGC
  array<size_t, kNumBases> n_base_quality;  // N

  /*********** PER QUALITY VALUE METRICS ****************/
  // Counts of quality in each base position
  array<size_t, kNumQualityValues * kNumBases> position_quality_count;

  // Counts of average quality (truncated) per sequence
  array<size_t, kNumQualityValues> quality_count;

  /*********** PER GC VALUE METRICS ****************/
  // histogram of GC fraction in each read from 0 to 100%
  array<size_t, 101> gc_count;

    /*********** PER READ METRICS ***************/
  // Distribution of read lengths
  array<size_t, kNumBases> read_length_freq;
  array<size_t, kNumBases> cumulative_read_length_freq;


  /*********** PER TILE SEQUENCE QUALITY OVERSERQUENCES ********/
  array <size_t, kNumMaxTiles * kNumBases> tile_position_quality;
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
  vector<size_t> long_base_quality;
  vector<size_t> long_n_base_quality;
  vector<size_t> long_position_quality_count;
  vector<size_t> long_read_length_freq;
  vector<size_t> long_cumulative_read_length_freq;
  vector<size_t> long_tile_position_quality;
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
         pass_kmer_content,
         pass_adapter_content;

    /**************** FUNCTIONS ****************************/

  // Default constructor that zeros everything
  explicit FastqStats(const size_t _kmer_size);

  // Allocation of more read positions
  inline void allocate_new_base(const bool ignore_tile);

  /******* DUPLICATION AND OVERREPRESENTATION *******/
  // Makes a hash map with keys as 32-bit suffixes and values as all the
  // candidate frequent sequences with that given suffix

  // Summarize all statistics we need before writing
  void summarize();

  // Writes to outpute fastqc-style
  void write(ostream &os, string filename);
};

// Default constructor
FastqStats::FastqStats(const size_t _kmer_size) {
  total_bases = 0;
  num_extra_bases = 0;
  avg_read_length = 0;
  avg_gc = 0;
  num_reads = 0;
  max_read_length = 0;
  num_poor = 0;

  num_unique_seen = 0;
  count_at_limit = 0;

  // Initialize IO arrays
  base_count.fill(0);
  n_base_count.fill(0);
  base_quality.fill(0);
  n_base_quality.fill(0);
  read_length_freq.fill(0);
  quality_count.fill(0);
  gc_count.fill(0);
  position_quality_count.fill(0);
  tile_position_quality.fill(0);
  tile_count.fill(0);

  // Defines bit shift values
  kBitShiftNucleotide = log2exact(kNumNucleotides);
  kBitShiftQuality = log2exact(kNumQualityValues);
  kBitShiftTile = log2exact(kNumMaxTiles);

  // Defines k-mer mask, length and allocates vector
  kmer_size = _kmer_size;
  kmer_mask = (1ll << (2*kmer_size)) - 1;
  kmer_count = vector<size_t>(min(kNumBases, kKmerMaxBases) * (kmer_mask + 1), 0);
}

// When we read new bases, dynamically allocate new space for their statistics
inline void
FastqStats::allocate_new_base(const bool ignore_tile) {
  
  for (size_t i = 0; i < kNumNucleotides; ++i)
    long_base_count.push_back(0);
  long_n_base_count.push_back (0);

  // space to allocate quality average
  for (size_t i = 0; i < kNumNucleotides; ++i)
    long_base_quality.push_back(0);
  long_n_base_quality.push_back(0);

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
FastqStats::summarize() {
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
  for (size_t i = 0; i < kPoorQualityThreshold; ++i)
    num_poor += quality_count[i];

  // Cumulative read length frequency
  size_t cumulative_sum = 0;
  for (size_t i = 0; i < max_read_length; ++i) {
    if (i < kNumBases)
      cumulative_sum += read_length_freq[i];
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
      if (cur_lquartile < 5)
        pass_per_base_sequence_quality = "fail";
      else if (cur_lquartile < 10)
        pass_per_base_sequence_quality = "warn";

      if (cur_median < 20)
        pass_per_base_sequence_quality = "fail";
      else if (cur_median < 25)
        pass_per_base_sequence_quality = "warn";
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
    if (i < kNumBases) {
      g_pct[i] = 100.0*g / total;
      a_pct[i] = 100.0*a / total;
      t_pct[i] = 100.0*t / total;
      c_pct[i] = 100.0*c / total;
      n_pct[i] = 100.0*n / total;
    } else {
      long_g_pct.push_back(100.0*g / total);
      long_a_pct.push_back(100.0*a / total);
      long_t_pct.push_back(100.0*t / total);
      long_c_pct.push_back(100.0*c / total);
      long_n_pct.push_back(100.0*n / total);
    }

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

  /******************* PER SEQUENCE GC CONTENT *****************/
  pass_per_sequence_gc_content = "pass";
  // fix the zero gcs by averaging around the adjacent values
  for (size_t i = 1; i < 99; ++i) 
    if (gc_count[i] == 0)
      gc_count[i] = (gc_count[i+1] + gc_count[i-1])/2;

  // Calculate pass warn fail statistics
  cerr << "gc_deviation: " << gc_deviation << "\n";
  if (gc_deviation >= 0.3)
    pass_per_sequence_gc_content = "fail";
  else if (gc_deviation >= 0.15)
    pass_per_sequence_gc_content = "warn";

  /******************* PER BASE N CONTENT **********************/

  pass_per_base_n_content = "pass";
  double cur_n_pct;
  for (size_t i = 0; i < max_read_length; ++i) {
    if (pass_per_base_n_content != "fail") {
      if (i < kNumBases)
        cur_n_pct = n_pct[i];
      else
        cur_n_pct = long_n_pct[i - kNumBases];

      if (cur_n_pct > 20.0) {
        pass_per_base_n_content = "fail";
      } else if (cur_n_pct > 10.0)
        pass_per_base_n_content = "warn";
    }
  }

  /************** SEQUENCE LENGTH DISTRIBUTION *****************/
  pass_sequence_length_distribution = "pass";
  size_t freq_of_avg;
  
  if (avg_read_length < kNumBases)
    freq_of_avg = read_length_freq[avg_read_length];
  else
    freq_of_avg = long_read_length_freq[avg_read_length - kNumBases];

  if (freq_of_avg != num_reads)
    pass_sequence_length_distribution = "warn";
  if (read_length_freq[0] > 0)
    pass_sequence_length_distribution = "fail";

  /************** DUPLICATE SEQUENCES **************************/
  pass_duplicate_sequences = "pass";

  double seq_total = 0.0;
  double seq_dedup = 0.0;
  unordered_map <size_t, size_t> counts_by_freq;

  for (auto v : sequence_count) {
    if(counts_by_freq.count(v.second) == 0)
      counts_by_freq[v.second] = 0;
    counts_by_freq[v.second]++;
  }

  for (auto v: counts_by_freq) {
    counts_by_freq[v.first] = get_corrected_count (count_at_limit, 
                                                             num_reads, 
                                                             v.first, 
                                                             v.second);
  }

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

  total_deduplicated_pct = 100.0 * seq_dedup / seq_total;

  // Convert to percentage
  for (auto &v : percentage_deduplicated)
    v = 100.0 * v / seq_dedup;  // Percentage of unique sequences in bin

   // Convert to percentage
  for (auto &v : percentage_total)
    v = 100.0 * v / seq_total;  // Percentage of sequences in bin

  // pass warn fail criteria : unique reads must be >80% 
  // (otherwise warn) or >50% (otherwisefail)
  if (percentage_total[0] <= 50)
    pass_duplicate_sequences = "fail";
  else if (percentage_total[0] <= 80)
    pass_duplicate_sequences = "warn";

  /************** OVERREPRESENTED SEQUENCES ********************/
  pass_overrepresented_sequences = "pass";

  // Keep only sequences that pass the input cutoff
  for (auto it = sequence_count.begin(); it != sequence_count.end(); ++it) {
    if (it->second > num_reads * kOverrepMinFrac)
      overrep_sequences.push_back (make_pair(it->first, 
                                              it->second));
  }

  // Sort strings by frequency
  sort(overrep_sequences.begin(), overrep_sequences.end(), 
       [](auto &a, auto &b){
          return a.second > b.second;
          });
  /************** ADAPTER CONTENT ******************************/
  pass_adapter_content = "pass";


  /************** KMER CONTENT *********************************/
  pass_kmer_content = "pass";

  /************** PER TILE SEQUENCE QUALITY ********************/
}

/****************** WRITE STATS ***********************/
void
FastqStats::write(ostream &os, string filename) {
  // Header
  os << "##FastQC\t0.11.8\n";

  // Basic statistics
  os << ">>Basic Statistics\t" << pass_basic_statistics << "\n";
  os << "#Measure\tValue\n";
  os << "Filename\t" << filename << "\n";
  os << "File type\tConventional base calls\n";
  os << "Total Sequences\t" << num_reads << "\n";
  os << "Sequences flagged as poor quality \t" << num_poor << "\n";
  os << "%GC \t" << avg_gc << "\n";
  os << ">>END_MODULE\n";

  // Per base quality
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

  // Per sequence quality scores
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

  // Per tile sequence quality
  os << ">>Per tile sequence quality\t" <<
        pass_per_tile_sequence_quality << "\n";

  size_t tile_pos_ind;
  for (size_t i = 0; i < kNumMaxTiles; ++i) {
    if (tile_count[i] > 0) {
      for (size_t j = 0; j < kNumBases; ++j) {
        tile_pos_ind = (j << kBitShiftTile) | i;
        if(tile_position_quality[tile_pos_ind] > 0) {
          os << i << "\t" << j + 1 << "\t" << 
            tile_position_quality[tile_pos_ind] / 
            static_cast<double>(tile_count[i]) << "\n";
        }
      }
    }
  }

  os << ">>END_MODULE\n";
  // Per sequence gc content
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

  // Sequence length distribution
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

  // Sequence duplication levels
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

  // Overrepresented sequences
  os << ">>Overrepresented sequences\t" <<
        pass_overrepresented_sequences << "\n";
  os << "#Sequence\tCount\tPercentage\tPossible Source\n";

  for (auto seq : overrep_sequences)
      os << seq.first << "\t" << seq.second <<  "\t" <<
        100.0 * seq.second / num_reads << "\t???\n";
  os << ">>END_MODULE\n";
  os << ">>Adapter Content\t" << pass_adapter_content << "\n";
  os << ">>END_MODULE\n";
}


/*************************************************************
 ******************** FASTQ READER ***************************
 *************************************************************/

struct FastqReader{
 public:
  // Memory map variables
  char *curr;  // current position in file
  char *last;  // last position in file
  struct stat st;

  // buffer size to store line 2 of each read
  size_t buffer_size;

  // Number of bases that have overflown the buffer
  size_t leftover_ind;

  // Whether or not to get bases from buffer when reading quality line
  bool read_from_buffer;

  // Whether or not to write bases to buffer when reading sequence
  bool write_to_buffer;

  // Get a base from the sequence line
  char base_from_buffer;

  /********* TILE PARSING ********/
  bool tile_ignore;  // Whether to just ignore per tile sequence quality
  size_t tile_cur;  // tile value parsed from line 1 of each record
  size_t tile_split_point;

  // Temp variables to be updated as you pass through the file
  size_t base_ind;  // 0,1,2 or 3
  size_t read_pos;  // which base we are at in the read
  size_t quality_value;  // to convert from ascii to number
  size_t cur_gc_count;  // Number of gc bases in read
  size_t cur_quality;  // Sum of quality values in read
  size_t num_bases_after_n;  // count of k-mers that reset at every N
  void *mmap_data;

  // Temporarily store line 2 out of 4 to know the base to which
  // quality characters are associated
  string buffer;
  string leftover_buffer;
  string sequence_to_hash;  // the sequence that will be marked for duplication

  size_t cur_kmer;  // 32-mer hash as you pass through the sequence line

  /************ FUNCTIONS TO PROCESS READS AND BASES ***********/
  inline void put_base();  // puts base in buffer or leftover
  inline void get_base();  // gets base from buffer or leftover

  // Check if it is a line that is done only occasionally
  inline bool is_tile_line(const FastqStats &stats);
  inline bool is_kmer_line(const FastqStats &stats);

  // gets and puts bases from and to buffer

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
  inline void read_tile_line(FastqStats &stats);  // get tile from read name
  inline void read_sequence_line(FastqStats &stats);  // parse sequence
  inline void read_quality_line(FastqStats &stats);  // parse quality

  explicit FastqReader(const size_t buffer_size);
  void memorymap(const string &filename,
                 const bool VERBOSE);
  void memoryunmap();
  inline bool operator >> (FastqStats &stats);
};

FastqReader::FastqReader(const size_t _buffer_size) {
  buffer_size = _buffer_size;

  // Allocates buffer to temporarily store reads
  buffer.resize(buffer_size + 1);
  buffer[buffer_size] = '\0';

  // Tile init
  tile_ignore= false;
  tile_cur= 0;
  tile_split_point = 0;
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

void
FastqReader::memoryunmap() {
  munmap(mmap_data, st.st_size);
  buffer.clear();
  leftover_buffer.clear();
}

// Fastqc only counts kmers once every 50 reads. We will do it once every 32
inline bool
FastqReader::is_kmer_line(const FastqStats &stats) {
  return (!(stats.num_reads & 31));
}

/*******************************************************/
/*************** BUFFER MANAGEMENT *********************/
/*******************************************************/
// puts base either on buffer or leftover
inline void
FastqReader::put_base() {
  base_from_buffer = *curr;
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
FastqReader::get_base() {
  if (read_from_buffer) {
    base_from_buffer = buffer[read_pos];
  } else {
    base_from_buffer = leftover_buffer[leftover_ind];
  }
}

/*******************************************************/
/*************** TILE PROCESSING ***********************/
/*******************************************************/

// Fastqc does tile statistics once every 10 lines, 
// we will do once every 8 so we can use bits instead of module arithmetic
inline bool
FastqReader::is_tile_line(const FastqStats &stats) {
  return (!(stats.num_reads & 7));
}

// Parse the comment 
inline void 
FastqReader::get_tile_split_position(){
  size_t num_colon = 0;

  // Count colons to know the formatting pattern
  for (; *curr != '\n'; ++curr) {
    if (*curr == ':')
      ++num_colon;
  }

  // Copied from fastqc
  if (num_colon >= 6)
    tile_split_point = 4;
  else if (num_colon >=4)
    tile_split_point = 2;

  // We will not get a tile out of this
  else
    tile_ignore= true;

}

inline void
FastqReader::get_tile_value() {
  tile_cur = 0;
  size_t num_colon = 0;
  for(; *curr != '\n'; ++curr) {
    if (*curr == ':')
      ++num_colon;

    if (num_colon == tile_split_point) {
      ++curr;  // pass the colon

      // parse till next colon or \n
      for(; (*curr != ':') && (*curr != '\n'); ++curr)
        tile_cur = tile_cur*10 + (*curr - '0');

      ++num_colon;
    }
  }
}

// Gets the tile from the sequence name (if applicable)
inline void
FastqReader::read_tile_line(FastqStats &stats){
  // if there is no tile information in the fastq header, fast
  // forward this line
  if (tile_ignore) {
    read_fast_forward_line();
    return;
  }

  // Fast forward if this is not a tile line
  if (!is_tile_line(stats)) {
    read_fast_forward_line();
    return;
  }

  // We haven't parsed the first line to know the split point
  if (tile_split_point == 0) 
    get_tile_split_position();
  else
    get_tile_value();

  ++curr;  // skip newline
}

// Skips lines that are not relevant
inline void
FastqReader::read_fast_forward_line(){
  for (; *curr != '\n'; ++curr) {}
  ++curr;  // skips \n from line 1
}

/*******************************************************/
/*************** SEQUENCE PROCESSING *******************/
/*******************************************************/

// This is probably the most important function for speed, so it must be really
// optimized at all times
inline void
FastqReader::process_sequence_base_from_buffer(FastqStats &stats) {
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

    if (is_kmer_line(stats)) {
      if (read_pos < stats.kKmerMaxBases) {
        cur_kmer = ((cur_kmer << stats.kBitShiftNucleotide) | base_ind);

        // registers k-mer if we've seen at least k nucleotides since the last n
        if (num_bases_after_n == stats.kmer_size) {
          stats.kmer_count[(read_pos << stats.kmer_size) | (cur_kmer & stats.kmer_mask)]++;
        }

        else
          num_bases_after_n++;
      }
    }
  }
}

// slower version of process_sequence_base_from_buffer that dynamically
// allocates
inline void
FastqReader::process_sequence_base_from_leftover(FastqStats &stats) {
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

    size_t test = (leftover_ind << stats.kBitShiftNucleotide) | base_ind;
    stats.long_base_count[test]++;

    // WE WILL NOT DO KMER STATS OUTSIDE OF BUFFER
  }
}

// Gets statistics after reading the entire sequence line
inline void
FastqReader::postprocess_sequence_line(FastqStats &stats) {
  // read length frequency histogram
  if((read_pos != 0) && (read_pos <= stats.kNumBases))
    stats.read_length_freq[read_pos - 1]++;
  else
    stats.long_read_length_freq[leftover_ind - 1]++;

  // Updates maximum read length if applicable
  if (read_pos > stats.max_read_length)
    stats.max_read_length = read_pos;

  // Registers GC % in the bin truncated to the nearest integer
  stats.gc_count[round(100 * cur_gc_count / static_cast<double>(read_pos))]++;

}

// Reads the line that has the biological sequence
inline void
FastqReader::read_sequence_line(FastqStats &stats){
  // restart line counters
  read_pos = 0;
  cur_gc_count = 0;
  num_bases_after_n = 1;
  write_to_buffer = true;
  leftover_ind = 0;

  /*********************************************************/
  /********** THIS LOOP MUST BE ALWAYS OPTIMIZED ***********/
  /*********************************************************/
  for (; *curr != '\n'; ++curr) {
    // puts base either on buffer or leftover 
    put_base();
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
    if (!write_to_buffer)
      leftover_ind++;
      
    // either way increase read position
    ++read_pos;
    
    // if we reached the buffer size, stop
    if(read_pos == buffer_size) {
      write_to_buffer = false;
    }
  }
  ++curr;  // skip \n

  // statistics summarized after the read
  postprocess_sequence_line(stats);
}


/*******************************************************/
/*************** QUALITY PROCESSING ********************/
/*******************************************************/

// Process quality value the fast way from buffer
inline void
FastqReader::process_quality_base_from_buffer(FastqStats &stats) {
  // N quality stats
  if (base_from_buffer == 'N') {
    stats.n_base_quality[read_pos] += quality_value;

 // ATGC QUALITY STATS
  } else {
    base_ind = actg_to_2bit(base_from_buffer);
    stats.base_quality[(read_pos << stats.kBitShiftNucleotide) | base_ind] 
         += quality_value;
  }

  // Average quality in position
  stats.position_quality_count[
        (read_pos << stats.kBitShiftQuality) | quality_value]++;

  // Tile processing
  if (!tile_ignore)
    if (is_tile_line(stats) && tile_cur > 0) {
      stats.tile_position_quality[(read_pos << stats.kBitShiftTile) | tile_cur]
          += quality_value;

  }
}

// Slow version of function above
inline void
FastqReader::process_quality_base_from_leftover(FastqStats &stats) {

  // N quality stats
  if (base_from_buffer == 'N') {
    stats.long_n_base_quality[leftover_ind] += quality_value;

 // ATGC QUALITY STATS
  } else {
    base_ind = actg_to_2bit(base_from_buffer);

    stats.long_base_quality[(leftover_ind << stats.kBitShiftNucleotide) | base_ind] 
         += quality_value;
  }

  // Average quality in position
  stats.long_position_quality_count[
        (leftover_ind << stats.kBitShiftQuality) | quality_value]++;

  // Tile processing
  if (!tile_ignore)
    if (is_tile_line(stats) && tile_cur > 0) {

      stats.long_tile_position_quality[(leftover_ind << stats.kBitShiftTile) | tile_cur]
          += quality_value;

  }
}

// Reads the quality line of each base.
inline void
FastqReader::read_quality_line(FastqStats &stats){
  // reset quality counts
  read_pos = 0;
  cur_quality = 0;
  read_from_buffer = true;
  leftover_ind = 0;

  for (; (*curr != '\n') && (curr < last); ++curr) {
    get_base();

    // Converts quality ascii to zero-based
    quality_value = *curr - stats.kBaseQuality;

    // Fast bases from buffer 
    if (read_from_buffer) 
      process_quality_base_from_buffer(stats);

    // Slow bases from dynamic allocation
    else
      process_quality_base_from_leftover(stats);

    // Sums quality value so we can bin the average at the end
    cur_quality += quality_value;

    if (!read_from_buffer)
      ++leftover_ind;

    // Flag to start reading and writing outside of buffer
    ++read_pos;
    if(read_pos == buffer_size) {
      read_from_buffer = false;
    }
  }
  ++curr;  // skip \n

  // Average quality approximated to the nearest integer. Used to make a
  // histogram in the end of the summary.
  stats.quality_count[cur_quality / read_pos]++;  // avg quality histogram
}

/*******************************************************/
/*************** POST LINE PROCESSING ******************/
/*******************************************************/

/*************** THIS IS VERY SLOW ********************/
inline void
FastqReader::postprocess_fastq_record(FastqStats &stats) {
  if(read_pos <= stats.kDupReadMaxSize)
    sequence_to_hash = buffer.substr(0, read_pos);
  else
    sequence_to_hash = buffer.substr(0, stats.kDupReadTruncateSize);

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

  // counts tile if applicable
  if (!tile_ignore)
    if (is_tile_line(stats) && tile_cur > 0)
      stats.tile_count[tile_cur]++;
  
}

/*******************************************************/
/*************** READ FASTQ RECORD *********************/
/*******************************************************/

inline bool
FastqReader::operator >> (FastqStats &stats) {
  read_tile_line(stats);
  read_sequence_line(stats);
  read_fast_forward_line();
  read_quality_line(stats);
  postprocess_fastq_record(stats);

  // Successful read, increment number in stats
  stats.num_reads++;

  // Returns if file should keep being checked
  return (curr < last - 1);

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
    cerr << "K-mer size should be larger than 2\n";
    cerr << opt_parse.help_message() << endl;
    return EXIT_SUCCESS;
  }
  filename = leftover_args.front();

  /****************** END COMMAND LINE OPTIONS *****************/
  if (VERBOSE)
    cerr << "Started reading file " << filename << ".\n";

  // Allocates vectors to summarize data
  FastqStats stats(kmer_size);

  // Initializes a reader given the maximum number of bases to summarie
  FastqReader in (stats.kNumBases);
  in.memorymap(filename, VERBOSE);

  // Read record by record
  const size_t num_reads_to_log = 1000000;
  size_t next_read = num_reads_to_log;
  while (in >> stats) {
    if(VERBOSE)
      // Equality is faster than modular arithmetics
      if (stats.num_reads == next_read) {
        cerr << "Processed " << stats.num_reads / num_reads_to_log
             << "M reads.\n";
        next_read += num_reads_to_log;
      }
  }

  // Free memory
  in.memoryunmap();
  if (VERBOSE)
    cerr << "Finished reading file.\n";

  if (VERBOSE)
    cerr << "Summarizing data.\n";

  // This function has to be called before writing to output
   stats.summarize();

  /************************ WRITE TO OUTPUT *****************************/
  if (VERBOSE)
    cerr << "Writing data.\n";

  // define output
  ofstream of;
  if (!outfile.empty())
    of.open(outfile.c_str(), ofstream::binary);

  ostream os(outfile.empty() ? cout.rdbuf() : of.rdbuf());

  // Write
  stats.write(os, filename);
  /************** TIME SUMMARY *********************************/
  if (VERBOSE)
    cerr << "Elapsed time: "
         << (clock() - begin) / CLOCKS_PER_SEC
         << " seconds\n";
}

