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
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include "fm_wrapper.hpp"
#include "util.hpp"
#include "positions_vector.hpp"
#include "diff_vector.hpp"
#include "bloom_filter.hpp"

using namespace std;

typedef PositionsVector<sd_vector<>, bit_vector, sd_vector<>> t_pv_vector;

class CRIndex {
    public:
        static const bool DEFAULT_VERBOSITY;
        static const int DEFAULT_READ_LENGTH;
        static bool verbose;

        CRIndex(string path, int read_length = DEFAULT_READ_LENGTH,
                bool verbose = DEFAULT_VERBOSITY, bool crappy_preprocess = false);
        CRIndex(string superstring, vector<t_pos> positions, vector<t_diff> diff,
                int read_length = DEFAULT_READ_LENGTH,
                bool verbose = DEFAULT_VERBOSITY);
        vector<int> find_indexes(const string& s);
        vector<string> find_reads(const string& s);
        ~CRIndex();

        static tuple<string, vector<t_pos>, vector<t_diff>>
                     preprocess(string path, bool verbose = DEFAULT_VERBOSITY);
        static tuple<string, vector<t_pos>, vector<t_diff>>
                     crappy_preprocess(string path, bool verbose = DEFAULT_VERBOSITY);

        void SaveToFile(string filename);
        static CRIndex* LoadFromFile(string filename, bool verbose = DEFAULT_VERBOSITY);
    private:
        int read_length;
/*        vector<t_pos> positions;
        vector<t_diff> diff;*/
        t_pv_vector pv_vector;
        DiffVector diff_vector;
        BloomFilter bloom_filter;
        FMWrapper<> fm_index;

        void FinishInit(string supersting, vector<t_pos> positions,
                        vector<t_diff> diff);
        vector<t_pos> locate_positions(const string& s);
        vector<t_pos> locate_positions2(const string& s, const string& check_s);
        string extract_original_read(t_pos read);

        static void debug(string msg);
        static void debug(vector<string> msg);
        static void debug(vector<int> msg);
        static void info(string msg);

        CRIndex() {}
};
