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

void go(string index_filename, string bac_filename) {
    cout << "Building CRIndex ... " << endl;
    CRIndex *rm = CRIndex::LoadFromFile(index_filename);
    cout << "Referenced memory: " << get_referenced_memory_size() << "kB" << endl;
//    rm->SaveToFile("saving-test.index");
    malloc_trim(42);
    cout << "Referenced memory: " << get_referenced_memory_size() << "kB" << endl;
    cout << "Querying ..." << endl;

    ifstream bac(bac_filename);
    string line;
    string buf;
    while (getline(bac, line)) {
      if (line[0] == '>') continue;
      boost::to_upper(line);
      buf += line;
    }
    vector<string> queries;
    for (int i = 0; i < 100000; i++) {
      string q;
      while (true) {
        q = buf.substr(rand()%(buf.size() - 20), 15);
        bool good = true;
        for (int j = 0; j < q.size(); j++) {
          if (q[j] == 'A' || q[j] == 'C' || q[j] == 'G' || q[j] == 'T') continue;
          good = false;
          break;
        }
        if (good) break;
      }
      queries.push_back(q);
    }
    chrono::time_point<std::chrono::system_clock> t1, t2;
    chrono::duration<double> elapsed;
    t1 = std::chrono::system_clock::now();
    int total_found = 0;
    for (size_t i = 0; i < queries.size(); i++) {
      total_found += rm->find_indexes(queries[i]).size();
    }
    cout << "total found" << total_found << endl;
    t2 = std::chrono::system_clock::now();
    elapsed = t2 - t1;
    cout << "Querying took " << elapsed.count() << "s" << endl;
}

int main(int argc, char** argv) {
    if (argc < 3) {
        cerr << "Usage: " << argv[0] << " <filename>" << endl;
        cerr << "Examples: " << endl;
        cerr << "  Read index file and construct CR-index" << endl;
        cerr << "  " << argv[0] << " index.data" << endl;
        exit(1);
    }

    try {
        go(argv[1], argv[2]);
    } catch(exception &e) {
        cerr << "Error: " << e.what() << endl;
        exit(1);
    }

    return 0;
}
