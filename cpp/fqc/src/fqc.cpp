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
using std::pair;
using std::make_pair;
using std::ifstream;

/*************************************************************
 ******************** AUX FUNCTIONS **************************
 *************************************************************/

// converts 64 bit integer to a sequence string by reading 3 bits at a time and
// converting back to ACTGN
static string
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

/*************************************************************
 ******************** FASTQ STATS ****************************
 *************************************************************/

struct FastqStats {
 public:
  /*************** UNIVERSAL CONSTANTS ************************/

  // This can be a very large number, the maximum illumina read size. Does not
  // affect memory very much
  static const size_t kNumBases = 1000;

  // Value to subtract quality characters to get the actual quality value
  static const size_t kBaseQuality = 33;
  static const size_t kNumChars = 4;  // A = 000,C = 001,T = 010,G = 011,N = 111
  static const size_t kNumAsciiChars = 256;

  // threshold for a sequence to be considered  poor quality
  static const size_t kPoorQualityThreshold = 20;

  // Kmer size given as input
  size_t kmer_size;

  // 4^kmer size
  size_t kNumKmers;

  // mask to get only the first 2*k bits of the sliding window
  size_t kmer_mask;

  /*********** CONSTANT NUMBERS FROM THE ENTIRE FASTQ ****************/
  size_t total_bases;  // sum of all bases in all reads
  size_t avg_read_length;  // average of all read lengths
  size_t avg_gc;  // sum of g bases + c bases
  size_t avg_n;  // n bases
  size_t num_reads;
  size_t num_poor;

  /*********** PER BASE METRICS ****************/

  // counts the number of bases in every read position
  array<size_t, kNumChars * kNumBases> base_count;

  // counts the number of N bases in every read position
  array<size_t, kNumChars * kNumBases> n_base_count;

  // Sum of base qualities in every read position
  array<size_t, kNumChars * kNumBases> base_quality;

  // Sum of N nucleotide base qualities in every position
  array<size_t, kNumChars * kNumBases> n_base_quality;

  /*********** PER QUALITY VALUE METRICS ****************/
  // Counts of quality in each base position
  array<size_t, kNumAsciiChars * kNumBases> position_quality_count;

  // Counts of average quality (truncated) per sequence 
  array<size_t, kNumAsciiChars> quality_count;

  /*********** PER GC VALUE METRICS ****************/
  array<size_t, 101> gc_count;

  /********** KMER FREQUENCY ****************/
  // A 2^K + 1 vector to count all possible kmers
  vector<size_t> kmer_count;

  /*********** PER READ METRICS ***************/
  // count reads based on 21-base prefix
  vector<size_t> read_freq;

  // Distribution of read lengths
  array<size_t, kNumBases> read_length_freq;

  // I need this to know what to divide each base by
  array<size_t, kNumBases> cumulative_read_length_freq;

///////////////// FUNCTIONS /////////////////////////////////
  // Default constructor that zeros everything
  explicit FastqStats(const size_t _kmer_size);

  // Summarize all statistics we need before writing
  void calc_summary();

  // Writing functions
  // HEADER
  void write_header(ostream &os);

  // BASIC STATISTICS: columns = measure, value
  void write_basic_statistics(ostream &os);

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

FastqStats::FastqStats(const size_t _kmer_size) {
  total_bases = 0;
  avg_read_length = 0;
  avg_gc = 0;
  avg_n = 0;
  num_reads = 0;
  num_poor = 0;

  // Initialize arrays
  base_count.fill(0);
  n_base_count.fill(0);
  base_quality.fill(0);
  n_base_quality.fill(0);
  read_length_freq.fill(0);

  // Defines k-mer mask, length and allocates vector
  kmer_size = _kmer_size;
  kNumKmers = (1ll << (2*kmer_size));
  kmer_mask = kNumKmers - 1;
  kmer_count = vector<size_t>(kmer_mask + 1, 0);
  read_freq = vector<size_t>(kmer_mask + 1, 0);
}

void
FastqStats::calc_summary() {
  // 1) Average read length
  avg_read_length = 0;
  for (size_t i = 0; i < kNumBases; ++i)
    avg_read_length += i * read_length_freq[i];
  avg_read_length /= num_reads;

  // 2) GC %
  avg_gc = 0;
  for (size_t i = 0; i < kNumBases; ++i)
    avg_gc += gc_count[i];
  avg_gc /= num_reads;

  // 3) Poor quality reads
  num_poor = 0;
  for (size_t i = 0; i < kPoorQualityThreshold + kBaseQuality; ++i)
    num_poor += quality_count[i];

  // 4) Cumulative read length frequency
  size_t cumulative_sum = 0;
  for (int i = kNumBases - 1; i >= 0; --i) {
    cumulative_sum += read_length_freq[i];
    cumulative_read_length_freq[i] = cumulative_sum;
  }
}

/*************************************************************
 ******************** FASTQ READER ***************************
 *************************************************************/

struct FastqReader{
 private:
  // this is begging for a bug but it's repeated in two classes
  static const size_t kNumBases = 1000;

  // Memory map variables
  char *curr;  // current position in file
  char *last;  // last position in file
  struct stat st;

  // Temp variables to be updated as you pass through the file
  size_t base_ind;  // 0,1,2 or 3
  size_t read_pos;  // which base we are at in the read
  size_t cur_gc_count;
  size_t cur_quality;
  size_t kmer_pos;
  void *mmap_data;

  // Temporarily store the second read of the 4 lines to know the base to which
  // quality characters are associated
  array<size_t, kNumBases> buff;

  size_t cur_kmer;  // kmer hash

 public:
  FastqReader();
  void memorymap(const string &filename,
                 const bool VERBOSE);
  void memoryunmap();
  bool operator >> (FastqStats &stats);
};

FastqReader::FastqReader() {
  base_ind = 0;  // 0,1,2 or 3
  read_pos = 0;  // which base we are at in the read
  cur_gc_count = 0;
  cur_quality = 0;
  kmer_pos = 0;
  cur_kmer = 0;
  buff.fill(0);
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

bool
FastqReader::operator >> (FastqStats &stats) {
  // *************READ NAME LINE****************************/
  // fast forward first line
  for (; *curr != '\n'; ++curr) {}
    ++curr;  // skips \n from line 1

    // *************NUCLEOTIDE LINE***************************/
  for (cur_gc_count = 0, read_pos = 0; *curr != '\n'; ++curr, ++read_pos) {
    // Gets data for ATGC
    if (*curr != 'N') {
      // Transforms base into 3-bit index
      // Bits 2,3 and 4 of charcters A,C,G,T and N are distinct so we can just
      // use them to get the index instead of doing if-elses.
      base_ind = ((*curr) >> 1) & 3;

      // Increments count of base
      stats.base_count[(read_pos << 2) | base_ind]++;

      // Hash the current kmer
      cur_kmer = ((cur_kmer << 2) | base_ind);

      // number of gcs in the read, G and C have odd bit
      cur_gc_count += (base_ind & 1);

      // If we already read >= k bases, increment the k-mer count
      if (kmer_pos >= stats.kmer_size - 1)
        stats.kmer_count[cur_kmer & stats.kmer_mask]++;

      kmer_pos++;
    } else {
      // Gets data for N
      stats.n_base_count[read_pos]++;
      kmer_pos = 0;  // start over the current kmer
    }

    // Need this to know what base the quality is associated to, therefore
    // we will store each read in memory
    buff[read_pos] = base_ind;
  }

  ++curr;  // skips \n from line 2

  // Store read statistics
  stats.read_length_freq[read_pos - 1]++;  // registers the read length

  // Store GC value read (0 to 100%)
  stats.gc_count[100 *cur_gc_count / read_pos]++;

  // Adds suffix to read frequency
  stats.read_freq[cur_kmer & stats.kmer_mask]++;

  // *************OTHER READ NAME LINE**********************/
  // fast forward third line
  for (; *curr != '\n'; ++curr) {}
  ++curr;  // skips \n from line 3

  // *************QUALITY LINE******************************/
  cur_quality = 0;
  for (read_pos = 0; (*curr != '\n') && (curr < last); ++curr, ++read_pos) {
    // Gets the base from buffer, finds index, adds quality value
    if (buff[read_pos] == 'N')
      stats.n_base_quality[read_pos] += *curr;
    else
      stats.base_quality[(read_pos << 2) | buff[read_pos]] += *curr;
    cur_quality += *curr;

    // We store for every position x quality value pair the count
    stats.position_quality_count[(read_pos << 8) | *curr]++;
  }
  // Slow but idk how else to do it
  stats.quality_count [cur_quality / read_pos]++;

  // Successful read, increment number in stats
  stats.num_reads++;

  // Returns if file should keep being checked
  if (curr < last - 1) {
    ++curr;  // skips the \n from the 4th line
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
FastqStats::write_basic_statistics(ostream &os) {
  // TODO(gui): Check what defines pass or fail
  os << ">>Basic Statistics\tpass\n";
  os << "#Measure\tValue\n";
  os << "Total Sequences\t" << num_reads << "\n";
  os << "Sequences flagged as poor quality \t" << num_poor << "\n";
  os << ">>END_MODULE\n";
}

void
FastqStats::write_per_base_quality(ostream &os) {
  // TODO(gui): Check what defines pass or fail
  // TODO(gui): Make median, deciles and quartiles
  os << ">>Per base sequence quality pass\n";
  os << "#Base\tMean\tMedian\tLower Quartile\tUpper Quartile" << 
        "\t10th Percentile 90th Percentile\n";
  size_t ldecile, lquartile, median, uquartile, udecile;
  double mean,cur;

  for (size_t i = 0; i < kNumBases; ++i) {
    if (cumulative_read_length_freq[i] > 0) {
      mean = 0;
      size_t counts = 0;
      size_t nreads = cumulative_read_length_freq[i];

      // Number of counts I need to see to know in which bin each *ile is
      double ldecile_thresh = 0.1*nreads;
      double lquartile_thresh = 0.25*nreads;
      double median_thresh = 0.5*nreads;
      double uquartile_thresh = 0.75*nreads;
      double udecile_thresh = 0.9*nreads;

      // Iterate through quality values
      for (size_t j = 0; j < kNumAsciiChars; ++j) {
        cur = position_quality_count[(i << 8) + j];
        if(counts < ldecile_thresh && counts + cur >= ldecile_thresh)
          ldecile = j - kBaseQuality;
        if(counts < lquartile_thresh && counts + cur >= lquartile_thresh)
          lquartile = j - kBaseQuality;
        if(counts < median_thresh && counts + cur >= median_thresh)
          median = j - kBaseQuality;
        if(counts < uquartile_thresh && counts + cur >= uquartile_thresh)
          uquartile = j - kBaseQuality;
        if(counts < udecile_thresh && counts + cur >= udecile_thresh)
          udecile = j - kBaseQuality;
        mean += cur*j;
        counts += cur;
      }

      // Normalize mean 
      mean = mean / nreads - kBaseQuality;
      // Write distribution to new line
      os << i + 1 << "\t" << mean << "\t" << median << "\t" << lquartile << "\t"
         << uquartile << "\t" << ldecile << "\t" << udecile << "\n";
    }
  }
  os << ">>END_MODULE\n";
}

void
FastqStats::write_per_sequence_quality(ostream &os) {
  os << ">>Per sequence quality scores\tpass\n";
  os << "#Quality\tCount\n";

  for (size_t i = kBaseQuality; i < kNumAsciiChars; ++i) {
    if (quality_count[i] > 0)
      os << i - kBaseQuality << "\t" << quality_count[i] << "\n";
  }
  os << ">>END_MODULE\n";
}

void
FastqStats::write_per_base_sequence_content(ostream &os) {
  os << ">>Per base sequence content\tpass\n";
  os << "#Base\tG\tA\tT\tC\n";

  size_t a, t, g, c, n;
  double total;
  for (size_t i = 0; i < kNumBases; ++i) {
    if (cumulative_read_length_freq[i] > 0) {
      a = base_count[(i << 2)];
      c = base_count[(i << 2) + 1];
      t = base_count[(i << 2) + 2];
      g = base_count[(i << 2) + 3];
      n = n_base_count[i];
      total = static_cast<double>(a+c+t+g+n);
      os << i+1 << "\t" <<
            100.0*(g + 0.25*n) / total << "\t" <<
            100.0*(a + 0.25*n) / total << "\t" <<
            100.0*(t + 0.25*n) / total << "\t" <<
            100.0*(c + 0.25*n) / total << "\n";
    }
  }
  os << ">>END_MODULE\n";
}

void
FastqStats::write_per_sequence_gc_content(ostream &os) {
  os << ">>Per sequence gc content\tpass\n";
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
  os << ">>Per base N concent\tpass\n";
  os << "#Base\tN-Count\n";

  size_t a, t, g, c, n;
  double total;
  for (size_t i = 0; i < kNumBases; ++i) {
    if (cumulative_read_length_freq[i] > 0) {
      a = base_count[(i << 2)];
      c = base_count[(i << 2) + 1];
      t = base_count[(i << 2) + 2];
      g = base_count[(i << 2) + 3];
      n = n_base_count[i];
      total = static_cast<double>(a+c+t+g+n);
      os << i+1 << "\t" << 100.0*(n / total) << "\n";
    }
  }

  os << ">>END_MODULE\n";
}

void
FastqStats::write_sequence_length_distribution(ostream &os) {
  os << "Sequence Length Distribution\tpass\n";
  os << "Length\tCount\n";
  for (size_t i = 0; i < kNumBases; ++i)
    if (read_length_freq[i] > 0)
      os << i+1 << "\t" << read_length_freq[i] << "\n";
  os << ">>END_MODULE\n";
}

void
FastqStats::write_sequence_duplication_levels(ostream &os) {
  os << ">>Sequence Duplication Levels\tpass\n";
  os << "#Duplication Level  Percentage of deduplicated  Percentage of total\n";
  os << ">>END_MOUDLE\n";
}

void
FastqStats::write_overrepresented_sequences(ostream &os) {
  os << ">>Overrepresented sequences warn\n";
  os << "#Sequence Count Percentage  Possible Sourc\n";
  os << ">>END_MODULE\n";
}

void
FastqStats::write_adapter_content(ostream &os) {
  os << ">>Adapter Content\tpass\n";
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
  size_t kmer_size = 8;
  const size_t MAX_KMER_SIZE = 16;

  OptionParser opt_parse(strip_path(argv[0]),
                         "Quality control metrics for fastq files",
                         "<fastq-file>");

  opt_parse.add_opt("kmer", 'k',
                    "k-mer size (default = 8, max = 10)", false, kmer_size);

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
  if (VERBOSE)
    cerr << "Started reading file " << filename << ".\n";
  FastqReader in;
  in.memorymap(filename, VERBOSE);

  if (VERBOSE)
    cerr << "Started processing file " << filename << ".\n";

  FastqStats stats(kmer_size);
  while (in >> stats) {}

  in.memoryunmap();
  if (VERBOSE)
    cerr << "Finished reading file.\n";

  if (VERBOSE)
    cerr << "Summarizing data.\n";

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
  stats.write_basic_statistics(os);
  stats.write_per_base_quality(os);
  stats.write_per_sequence_quality(os);
  stats.write_per_base_sequence_content(os);
  stats.write_per_sequence_gc_content(os);
  stats.write_per_base_n_content(os);
  stats.write_sequence_length_distribution(os);
  stats.write_sequence_duplication_levels(os);
  stats.write_adapter_content(os);

  /************** TIME SUMMARY *********************************/
  // TODO(gui) : find adapters with significant kmer enrichment
  if (VERBOSE)
    cerr << "Elapsed time: "
         << (clock() - begin) / CLOCKS_PER_SEC
         << " seconds\n";
}

