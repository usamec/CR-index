#ifndef BLOOM_FILTER_HPP__
#define BLOOM_FILTER_HPP__

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
#include <cmath>
#include <boost/algorithm/string.hpp>
#include <boost/functional/hash.hpp>
#include <boost/filesystem.hpp>
#include <sdsl/bit_vectors.hpp>

using namespace std;
using namespace sdsl;

class BloomFilter {
  public:
    BloomFilter() {}
    BloomFilter(int items, double p = 0.1) {
      int m1 = - items * log(p) / log(2) / log(2);
      k = log(2) * m1 / items;
      int m = 1;
      sh = 0;
      while (m < m1) {m *= 2; sh += 1;}
      cout << "bloom size " << m << " " << m / 8 + 1 << endl;
      data.resize(m);
      util::set_to_value(data, 0);
    }

    void Add(const string& s1, const string& s2) {
      auto h1 = hh(s1);
      auto h2 = hh(s2);
      for (int i = 0; i < k; i++) {
        size_t h = k;
        boost::hash_combine(h, h1);
        boost::hash_combine(h, h2);
        h &= ((1 << sh) - 1);
        data[h] = 1;
      }
    }

    bool Test(const string& s1, const string& s2) const {
      auto h1 = hh(s1);
      auto h2 = hh(s2);
      for (int i = 0; i < k; i++) {
        size_t h = k;
        boost::hash_combine(h, h1);
        boost::hash_combine(h, h2);
        h &= ((1 << sh) - 1);
        if (data[h] == 0) return false;
      }
      return true;
    }

    int k;
    bit_vector data;
    hash<string> hh;
    int sh;
};

#endif
