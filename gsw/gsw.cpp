#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>
using std::string;
using std::cin;
using std::cout;
using std::cerr;
using std::endl;
using std::max;
using std::vector;
using std::istream;
using std::ifstream;
using std::ostream;
using std::getline;
using std::runtime_error;

struct Sequence {
  string name;
  string seq;

  Sequence() {
    seq = "";
  }

  string tostring() const {
    return ">" + name + "\n" + seq;
  }
};

istream & operator
>> (istream &lhs, Sequence &s) {
  string line;
  if (lhs.peek() != '>')
    throw runtime_error("malformatted fasta file");

  while(getline(lhs, s.name) &&
        s.name.size() > 0 &&
        s.name[0] != '>');

  // remove >
  s.name.erase(s.name.begin());

  // read sequence
  s.seq.clear();
  while (lhs.peek() != '>' && getline(lhs, line))
    s.seq.append(line);

  return lhs;
}

namespace alignment {
  static const int sa = 1; //match
  static const int sb = 3; //mismatch
  static const int sg = 5; //gap
}

// match-mismatch score
inline int
s(char a, char b) {
  return (alignment::sa * (a == b)) -
         (alignment::sb * (a != b));
}

template <typename N>
void
get_alignment_positions(const N &H,
                        size_t &ts,
                        size_t &target_end,
                        size_t &qs,
                        size_t &query_end,
                        const string &target,
                        const string &query,
                        int    &max_score) {
  max_score = 0;
  size_t i,j;

  // find largest score
  for (i = 0; i != H.size(); ++i)
    for (j = 0; j != H[0].size(); ++j)
      if (H[i][j] > max_score) {
        max_score = H[i][j];
        target_end = i;
        query_end = j;
      }

  ts = target_end;
  qs = query_end;
  while (H[ts][qs] != 0) {
    if (H[ts][qs] == H[ts-1][qs-1] + s(target[ts-1], query[qs-1])) {
      ts--;
      qs--;
    } else if (H[ts][qs] == H[ts-1][qs] - alignment::sg) {
      ts--;
    } else if (H[ts][qs] == H[ts][qs-1] - alignment::sg) {
      qs--;
    }
    // alignment ends here
    else return;
  }
}

// aligns two string, allowing a maximum of one indel during mapping
template <typename S,
          typename T,
          typename N> void
smith_waterman (const S &target,
                       const T &query,
                       N &H) {
  size_t i,j;
  // first row
  for (i = 1; i != target.size() + 1; ++i) {
    for (j = 1; j != query.size() + 1; ++j) {
      H[i][j] = 0;
      H[i][j] = max(H[i][j], H[i-1][j-1] + s(target[i-1], query[j-1]));
      H[i][j] = max(H[i][j], H[i-1][j] - alignment::sg);
      H[i][j] = max(H[i][j], H[i][j-1] - alignment::sg);
    }
  }
}

void
read_fasta(istream &is, vector<Sequence> &v) {
  Sequence s;
  while (is >> s)
    v.push_back(s);
  v.push_back(s);
}

int
main(int argc, char **argv) {
  if (argc != 3) {
    cerr << "run: gsw <target.fa> <query.fa>" << endl;
    return EXIT_SUCCESS;
  }
  vector <Sequence> targets, queries;
  ifstream target_file(argv[1]),
           query_file(argv[2]);
  read_fasta(target_file, targets);
  target_file.close();

  read_fasta(query_file, queries);
  query_file.close();

  // alignment matrix and results
  vector<vector<int>> H;
  size_t target_start, target_end, query_start, query_end;
  int max_score;

  for (size_t j = 0; j < queries.size(); ++j) {
    for (size_t i = 0; i < targets.size(); ++i) {
      // alignment matrix
      H = vector<vector<int>>(
          targets[i].seq.size()+1, vector<int>(
          queries[j].seq.size()+1, 0)
      );

      // run alignment
      smith_waterman(targets[i].seq,
                     queries[j].seq,
                     H);

      // get traceback
      get_alignment_positions(H, target_start, target_end,
                                 query_start, query_end,
                                 targets[i].seq,
                                 queries[j].seq,
                                 max_score);
      cout << targets[i].name << "\t" <<
              target_start << "\t" <<
              target_end << "\t" <<
              queries[j].name << "\t" <<
              query_start << "\t" <<
              query_end << "\t" <<
              max_score << "\n";

      /*
      cout << targets[i].seq << "\n";
      cout << queries[j].seq << "\n";

      cout << " ";
      for (size_t k = 0; k < H[0].size(); ++k)
        cout << queries[j].seq[k] << " ";
      cout << "\n";
      for (size_t k = 0; k < H.size(); ++k) {
        cout << targets[i].seq[k] << " ";
        for (size_t l = 0; l < H[0].size(); ++l) {
          cout << H[k][l] << " ";
        }
        cout << "\n";
      }
      */
    }
  }
}
