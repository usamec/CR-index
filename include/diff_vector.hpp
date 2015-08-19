#ifndef DIFF_VECTOR_HPP__
#define DIFF_VECTOR_HPP__

#include <string>
#include <vector>
#include <utility>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <limits>
#include <map>
#include <iostream>
#include <fstream>
#include <tuple>
#include <time.h>
#include <exception>
#include <stdexcept>
#include <unordered_map>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <sdsl/bit_vectors.hpp>

using namespace std;

typedef tuple<int, int, char> t_diff;

inline int BaseToInt(char c) {
  if (c == 'A') return 0;
  if (c == 'C') return 1;
  if (c == 'G') return 2;
  if (c == 'T') return 3;
  assert(false);
}

inline char IntToBase(int a) {
  char alp[] = "ACGT";
  return alp[a];
}

class DiffVector {
  typedef sd_vector<> t_vector;
  typedef t_vector::rank_1_type rank_type;
  typedef t_vector::select_1_type select_type;

  public:
    DiffVector() {
    }

    DiffVector operator=(const DiffVector& d) {
      data = d.data;
      rank = rank_type(&data);
      select = select_type(&data);
      diffs = d.diffs;
      read_length = d.read_length;
      return *this;
    }

    DiffVector(size_t read_count, size_t read_length_, 
               vector<t_diff> diffs_) :
        data(read_count*read_length_), read_length(read_length_), diffs(diffs_.size()) {
      bit_vector data_(read_count*read_length_);
      for (size_t i = 0; i < diffs_.size(); i++) {
        int read_id, read_pos;
        char change;
        tie(read_id, read_pos, change) = diffs_[i];
        int pos = read_id * read_length_ + read_pos;
        data_[pos] = 1;
        diffs[i] = BaseToInt(change);
      }
      data = t_vector(data_);
      rank = rank_type(&data);
      select = select_type(&data);
      sdsl::util::bit_compress(diffs);
    }

    int memory_size() {
      cout << "data   " << sdsl::size_in_bytes(data) << endl;
      cout << "rank   " << sdsl::size_in_bytes(rank) << endl;
      cout << "select " << sdsl::size_in_bytes(select) << endl;
      cout << "diffs  " << sdsl::size_in_bytes(diffs) << endl;
      return sdsl::size_in_bytes(data) +
             sdsl::size_in_bytes(rank) +
             sdsl::size_in_bytes(select) +
             sdsl::size_in_bytes(diffs);
    }

    vector<t_diff> GetDiffs(int read_id) {
      int start_pos = read_id * read_length;
      int end_pos = (read_id + 1) * read_length;
      int start_rank = rank(start_pos);
//      int end_rank = rank(end_pos);
      vector<t_diff> ret;
/*      for (int i = start_rank; i < end_rank; i++) {
        int pos = select(i+1);
        ret.push_back(make_tuple(read_id, pos - start_pos, IntToBase(diffs[i])));
      }*/
      for (int i = start_rank; i < diffs.size(); i++) {
        int pos = select(i+1);
        if (pos >= end_pos) break;
        ret.push_back(make_tuple(read_id, pos - start_pos, IntToBase(diffs[i])));
      }

      return ret;
    }

    int size() const {
      return diffs.size();
    }

    t_vector data;
    rank_type rank;
    select_type select;
    size_t read_length;
    int_vector<2> diffs;
};

#endif
