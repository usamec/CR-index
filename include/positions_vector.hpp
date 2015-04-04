#ifndef POSITIONS_VECTOR_HPP__
#define POSITIONS_VECTOR_HPP__

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

typedef tuple<int, int, bool> t_pos;

template <class t_bit_vector = bit_vector,
          class t_bit_vector2 = bit_vector,
          class t_bit_vector3 = bit_vector,
          class t_int_vector = int_vector<>>
class PositionsVector {
  typedef typename t_bit_vector::rank_1_type rank_type;
  typedef typename t_bit_vector::select_1_type select_type;
  typedef typename t_bit_vector3::rank_1_type rank_type3;
  public:
    PositionsVector() {}
    PositionsVector(const PositionsVector& v) :
      hits(v.hits), reverses(v.reverses), read_ids(v.read_ids), dups_bitmap(v.dups_bitmap),
      dups_ids(v.dups_ids), dups_reverses(v.dups_reverses) {
      hits_rank = rank_type(&hits);
      hits_select = select_type(&hits);
      dups_rank = rank_type3(&dups_bitmap);
    }
    PositionsVector& operator=(const PositionsVector& v) {
      hits = v.hits;
      reverses = v.reverses;
      read_ids = v.read_ids;
      hits_rank = rank_type(&hits);
      hits_select = select_type(&hits);
      dups_bitmap = v.dups_bitmap;
      dups_ids = v.dups_ids;
      dups_reverses = v.dups_reverses;
      dups_rank = rank_type3(&dups_bitmap);
      return *this;
    }

    PositionsVector(int superstring_size, vector<t_pos> positions) {
      bit_vector hits_(superstring_size, 0);
      bit_vector reverses_(positions.size());
      t_int_vector read_ids_(positions.size());
      int ind = 0;
      map<int, pair<int, bool>> dups;
      for (int i = 0; i < positions.size(); i++) {
        if (i > 0) {
          assert(get<0>(positions[i]) >= get<0>(positions[i-1]));
          if (get<0>(positions[i])== get<0>(positions[i-1])) {
            dups[get<1>(positions[i-1])] = make_pair(get<1>(positions[i]), get<2>(positions[i]));
            continue;
          }
        }
        hits_[get<0>(positions[i])] = 1;
        reverses_[ind] = get<2>(positions[i]);
        read_ids_[ind] = get<1>(positions[i]);
        ind++;
      }
      reverses_.resize(ind);
      read_ids_.resize(ind);
      hits = t_bit_vector(hits_);
      reverses = t_bit_vector2(reverses_);
      read_ids = t_int_vector(read_ids_);
      sdsl::util::bit_compress(read_ids);
      hits_rank = rank_type(&hits);
      hits_select = select_type(&hits);
      cout << "sizes " << hits.size() << " " << hits_rank.size() << " " << dups.size() << endl;    
      dups_ids.resize(dups.size());
      dups_reverses.resize(dups.size());
      bit_vector dups_bitmap_(positions.size());
      ind = 0;
      for (auto &e: dups) {
        dups_bitmap_[e.first] = 1;
        dups_reverses[ind] = e.second.second;
        dups_ids[ind] = e.second.first;
        ind++;
      }
      dups_bitmap = t_bit_vector3(dups_bitmap_);
      sdsl::util::bit_compress(dups_ids);
      dups_rank = rank_type3(&dups_bitmap);
    }

    vector<t_pos> GetRange(int start_index, int end_index) {
      int rank_start = hits_rank(start_index);
      int rank_end = hits_rank(end_index);
      vector<t_pos> ret;
      for (int i = rank_start; i < rank_end; i++) {
        int pos = hits_select(i+1);
        ret.push_back(make_tuple(pos, read_ids[i], (bool)reverses[i]));
        int cur_read = read_ids[i];
        while (dups_bitmap[cur_read] == 1) {
          int rank = dups_rank(cur_read);
          int next_id = dups_ids[rank];
          bool rev = dups_reverses[rank];
          ret.push_back(make_tuple(pos, next_id, rev));
          cur_read = next_id;
        }
      }
      return ret;
    }

    int memory_size() {
      cout << "pv hits      " << sdsl::size_in_bytes(hits) << endl;
      cout << "pv hits rank " << sdsl::size_in_bytes(hits_rank) << endl;
      cout << "pv hits sel  " << sdsl::size_in_bytes(hits_select) << endl;
      cout << "pv reverses  " << sdsl::size_in_bytes(reverses) << endl;
      cout << "pv read_ids  " << sdsl::size_in_bytes(read_ids) << endl;
      cout << "pv dups bit  " << sdsl::size_in_bytes(dups_bitmap) << endl;
      cout << "pv dups rank " << sdsl::size_in_bytes(dups_rank) << endl;
      cout << "pv dups ids  " << sdsl::size_in_bytes(dups_ids) << endl;
      cout << "pv dups revs " << sdsl::size_in_bytes(dups_reverses) << endl;
      return sdsl::size_in_bytes(hits) + 
             sdsl::size_in_bytes(hits_rank) + 
             sdsl::size_in_bytes(hits_select) + 
             sdsl::size_in_bytes(reverses) + 
             sdsl::size_in_bytes(read_ids) +
             sdsl::size_in_bytes(dups_ids) +
             sdsl::size_in_bytes(dups_bitmap) +
             sdsl::size_in_bytes(dups_rank) +
             sdsl::size_in_bytes(dups_reverses);
    }

    t_bit_vector hits;
    t_bit_vector2 reverses;
    t_int_vector read_ids;
    rank_type hits_rank;
    select_type hits_select;

    t_bit_vector3 dups_bitmap;
    t_int_vector dups_ids;
    t_bit_vector2 dups_reverses;
    rank_type3 dups_rank;
};

#endif
