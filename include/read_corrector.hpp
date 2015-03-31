#include <unordered_map>
#include <string>
#include "util.hpp"

using cr_util::rev_compl; 

class ReadCorrector {
  public:
    ReadCorrector(int kmer=31, int bad_kmer_limit=3) :
        kmer_(kmer), bad_kmer_limit_(bad_kmer_limit) {
    }

    void AddRead(string read) {
      for (int i = 0; i + kmer_ <= read.size(); i++) {
        kmer_counts_[read.substr(i, kmer_)]++;
        kmer_counts_[rev_compl(read.substr(i, kmer_))]++;
      }
    }

    string CorrectRead(string read) {
      char alpha[] = "ACGT";
      for (int i = 0; i + kmer_ <= read.size(); i++) {
        if (kmer_counts_[read.substr(i, kmer_)] < bad_kmer_limit_) {
          string bad_kmer = read.substr(i, kmer_);
          string best_correction = bad_kmer;
          string cur_kmer;
          for (int j = 0; j < kmer_; j++) {
            cur_kmer = bad_kmer;
            for (int k = 0; k < 4; k++) {
              cur_kmer[j] = alpha[k];
              if (kmer_counts_[cur_kmer] > kmer_counts_[best_correction]) {
                best_correction = cur_kmer;
              }
            }
          }
          for (int j = 0; j < kmer_; j++) {
            read[i+j] = best_correction[j];
          }
          i += kmer_ + 15;
        }
      }
      return read;
    }

    int kmer_;
    int bad_kmer_limit_;
    unordered_map<string, int> kmer_counts_;
};

