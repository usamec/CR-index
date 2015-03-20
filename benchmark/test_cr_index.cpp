#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <iostream>
#include <vector>
#include <string>
#include <utility>
#include <chrono>
#include <boost/algorithm/string.hpp>
#include <malloc.h>
#include "cr_index.hpp"

using namespace std;

bool is_memory_mapping(const string &line) {
    istringstream s(line);
    s.exceptions(istream::failbit | istream::badbit);
    try {
        // Try reading the fields to verify format.
        string buf;
        s >> buf; // address range
        s >> buf; // permissions
        s >> buf; // offset
        s >> buf; // device major:minor
        size_t inode;
        s >> inode;
        if (inode == 0) {
            return true;
        }
    } catch (std::ios_base::failure &e) { }

    return false;
}

size_t get_field_value(const string &line) {
    istringstream s(line);
    size_t result;
    string buf;
    s >> buf >> result;

    return result;
}

size_t get_referenced_memory_size() {
    ifstream smaps("/proc/self/smaps");
    smaps.exceptions(istream::failbit | istream::badbit);
    string line;
    size_t result = 0;

    try {
        while (1) {
            while (!is_memory_mapping(line)) {
                getline(smaps, line);
            }
            while (!(line.substr(0, 11) == "Referenced:")) {
                getline(smaps, line);
            }

            result += get_field_value(line);
        }
    } catch (std::ios_base::failure &e) { }

    return result;
}

void go(string index_filename) {
    cout << "Building CRIndex ... " << endl;
    CRIndex *rm = CRIndex::LoadFromFile(index_filename);
    cout << "Referenced memory: " << get_referenced_memory_size() << "kB" << endl;
    malloc_trim(42);
    cout << "Referenced memory: " << get_referenced_memory_size() << "kB" << endl;
    vector<int> res = rm->find_indexes("CCATCTT");
    cout << "Referenced memory: " << get_referenced_memory_size() << "kB" << endl;
    cout << res.size() << endl;
    cout << "Referenced memory: " << get_referenced_memory_size() << "kB" << endl;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <filename>" << endl;
        cerr << "Examples: " << endl;
        cerr << "  Read index file and construct CR-index" << endl;
        cerr << "  " << argv[0] << " index.data" << endl;
        exit(1);
    }

    try {
        go(argv[1]);
    } catch(exception &e) {
        cerr << "Error: " << e.what() << endl;
        exit(1);
    }

    return 0;
}
