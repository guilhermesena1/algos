#include "Config.hpp"
#include <fstream>
#include <sstream>
using std::string;
using std::vector;
using std::unordered_map;
using std::pair;
using std::make_pair;
using std::ifstream;
using std::runtime_error;
using std::istringstream;

// Check if a std::string ends with another, to be use to figure out the file format
static inline bool
endswith(std::string const & value, std::string const & ending) {
  if (ending.size() > value.size()) {
    return false;
  }
  return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
}

string
Config::strip_path(string full_path) const {
  size_t start = full_path.find_last_of('/');
  if (start == string::npos)
    start = 0;
  else
    ++start;
  return full_path.substr(start);
}

const vector<string> Config::values_to_check({
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
Config::Config() {
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
  contaminants_file = string(__MYFILE__) + "/Configuration/contaminant_list.txt";
  adapters_file = string(__MYFILE__) + "/Configuration/adapter_list.txt";
  limits_file = string(__MYFILE__) + "/Configuration/limits.txt";
  html_file = string(__MYFILE__) + "/Configuration/template.html";
  kmer_size = 7;
  quiet = false;
  tmpdir = ".";

  is_sam = false;
  is_bam = false;
  is_fastq = false;
  is_fastq_gz = false;

}

void
Config::setup() {
  define_file_format();
  read_limits();
  if (limits["adapter"]["ignore"] == 0.0)
    read_adapters();
  if (limits["adapter"]["ignore"] == 0.0)
    read_contaminants();

  filename_stripped = strip_path(filename);
}

void
Config::define_file_format() {
  if (format == "") {
    if (endswith(filename, "sam")) {
      is_sam = true;
    }

    if (endswith(filename, "bam")) {
      is_bam = true;
    }

    if (endswith(filename, "fastq.gz")) {
      is_fastq_gz = true;
    }

    if (endswith(filename, "fq.gz")) {
      is_fastq_gz = true;
    }

    if (endswith(filename, "fastq")) {
      is_fastq = true;
    }

    if (endswith(filename, "fq")) {
      is_fastq = true;
    }
  }
}

void
Config::read_limits() {
  ifstream in(limits_file);
  if (!in)
    throw runtime_error("limits file does not exist: " + limits_file);

  // Variables to parse lines
  string line, limit, instruction;
  double value;
  while (getline(in, line)) {
    // Lines with # are comments and should be skipped
    if (line[0] != '#') {
      istringstream iss(line);

      // Every line is a limit, warn/error/ignore and the value
      iss >> limit >> instruction >> value;

      if (find(values_to_check.begin(), values_to_check.end(), limit)
          == values_to_check.end())
        throw runtime_error("unknown limit option: " + limit);

      if (instruction != "warn" &&
          instruction != "error" &&
          instruction != "ignore")
        throw runtime_error("unknown instruction for limit " + limit +
                            ": " + instruction);

      limits[limit][instruction] = value;
    }
  }

  for (auto v : values_to_check) {
    if (limits.count(v) == 0)
      throw runtime_error("instruction for limit " + v +
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
Config::read_adapters() {
  ifstream in(adapters_file);
  if (!in)
    throw runtime_error("adapter file not found: " + adapters_file);

  string line, _tmp;
  vector<string> line_by_space;
  string adapter_name, adapter_seq;
  size_t adapter_hash;

  // The contaminants file has a space separated name, and the last instance is
  // the biological sequence
  while (getline(in, line)) {
    if (line[0] != '#') {
      adapter_name = "";
      adapter_seq = "";
      istringstream iss(line);
      while (iss >> _tmp) {
        line_by_space.push_back(_tmp);
      }

      if (line_by_space.size() > 1) {
        for (size_t i = 0; i < line_by_space.size() - 1; ++i)
          adapter_name += line_by_space[i] + " ";
        adapter_seq = line_by_space.back();

        if (adapter_seq.size() > kmer_size) {
          adapter_seq = adapter_seq.substr(0, kmer_size);
        }

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
Config::read_contaminants() {
  ifstream in(contaminants_file);
  if (!in)
    throw runtime_error("contaminants file not found: " + contaminants_file);

  string line, _tmp;
  vector<string> line_by_space;
  string contaminant_name, contaminant_seq;

  // The contaminants file has a space separated name, and the last instance is
  // the biological sequence
  while (getline(in, line)) {
    if (line[0] != '#') {
      contaminant_name = "";
      contaminant_seq = "";
      istringstream iss(line);
      while (iss >> _tmp) {
        line_by_space.push_back(_tmp);
      }

      if (line_by_space.size() > 1) {
        for (size_t i = 0; i < line_by_space.size() - 1; ++i) {
          contaminant_name += line_by_space[i] + " ";
        }
        contaminant_seq = line_by_space.back();
        contaminants.push_back(make_pair(contaminant_name, contaminant_seq));
      }

      line_by_space.clear();
    }
  }
  in.close();
}

// Find contaminant with highest overlap with sequence or return "No Hit" if
// there is none
string
Config::get_matching_contaminant(string seq) const {
  size_t best = 0;
  string ret;
  for (auto v : contaminants) {
    if (seq.size() > v.second.size()) {
      // contaminant contained in sequence
      if (seq.find(v.second) != string::npos) {
        if (v.second.size() > best) {
          best = v.second.size();
          ret = v.first;
        }
      }
    } else {
      // sequence contained in contaminant
      if (v.second.find(seq) != string::npos) {
        // In this case this is the best possible match so return it
        return v.first;
      }
    }
  }

  // If any sequence is a match, return the best one
  if (best > 0) {
    return ret;
  }
  return "No Hit";
}

