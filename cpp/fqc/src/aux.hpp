#ifndef _AUX_HPP
#define _AUX_HPP
#include <string>
#include <chrono>
#include <ctime>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <array>

/*************************************************************
 ******************** AUX FUNCTIONS **************************
 *************************************************************/
// converts 64 bit integer to a sequence std::string by reading 2 bits at a time and
// converting back to ACTG
static inline std::string
size_t_to_seq(size_t v, const size_t seq_length) {
  std::string ans;
  for (size_t i = 0; i < seq_length; ++i) {
    switch (v & 3) {
      case 0: ans.push_back('A'); break;
      case 1: ans.push_back('C'); break;
      case 2: ans.push_back('T'); break;
      case 3: ans.push_back('G'); break;
    }   
    v >>= 2;
  }

  std::reverse(ans.begin(), ans.end());
  return ans;
}

// Converts A,T,G,C to 2-bit values
static inline size_t
actg_to_2bit(const char &c) {
  return ((c >> 1) & 3);
}

// custom class to hash std::strings fast
class SmallSequence {
 public:
  size_t seq[3];
  size_t size;
  size_t i;
  size_t cur_ind;
  SmallSequence() {
    seq[0] = seq[1] = seq[2] = 0;
    size = 0;
  }

  inline void rewrite(const char *in, const size_t &num_chars) {
    cur_ind = 0;
    for (i = 0; i != num_chars; ++i) {
      seq[cur_ind] = (seq[cur_ind] << 2) | actg_to_2bit(in[i]);
      if (i == 31 || i == 63)
        ++cur_ind;
    }
  }

  std::string to_string() const {
    // we take all 96 bases even if we don't use them
    std::string ans = size_t_to_seq(seq[0], 32) +
                 size_t_to_seq(seq[1], 32) +
                 size_t_to_seq(seq[2], 32);

    // get just the prefix of the given std::string size
    return ans.substr(0, size);
  }

  inline bool operator ==(const SmallSequence &rhs) const {
    return (seq[0] == rhs.seq[0] &&
            seq[1] == rhs.seq[1] &&
            seq[2] == rhs.seq[2]);
  }
};

// The hash will be the first 32 bases
template <>
struct std::hash<SmallSequence> {
  inline size_t operator ()(const SmallSequence &k) const {
    return k.seq[0];
  }
};
#endif

