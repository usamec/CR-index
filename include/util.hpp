#ifndef UTIL_HPP__
#define UTIL_HPP__

#include <string>
#include <vector>
#include <fstream>
#include <stdexcept>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <pstreams/pstream.h>

using namespace std;
class BloomFilter;

namespace cr_util {
    bool check_read(const string& read);
    bool check_contig(const string& contig);
    string rev_compl(const string& s);
    string load_contigs(const string& contigs_path);
    string execute_command(const string &command);
    boost::filesystem::path create_tmpdir();
    void check_path_existence(boost::filesystem::path path);
    vector<int> diff_indexes(const string& s1, const string& s2);
    bool indexes_close(vector<int> indexes, int k);
    vector<string> strings_with_edt1(const string& s);
    vector<string> strings_with_edt1(const string& s, const BloomFilter &filter);
}

#endif
