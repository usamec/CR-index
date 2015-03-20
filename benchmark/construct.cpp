#include <iostream>
#include <string>
#include <chrono>
#include <malloc.h>
#include "cr_index.hpp"
#include "hash_index.hpp"
#include <libGkArrays/gkArrays.h>

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

bool test(string index_type, string reads_filename, int read_length,
        int query_length) {
    chrono::time_point<std::chrono::system_clock> t1, t2;
    chrono::duration<double> elapsed;

    cout << "index_type: " << index_type << endl;
    cout << "read_filename: " << reads_filename << endl;
    cout << "read_length: " << read_length << endl;
    cout << "query_length: " << query_length << endl;
    t1 = std::chrono::system_clock::now();
    if (index_type == "cr") {
        CRIndex cr = CRIndex(reads_filename, read_length, true);
        malloc_trim(42);
        cout << "Referenced memory: " << get_referenced_memory_size() << "kB ";
        cout << cr.find_indexes("CCATCCT").size() << " hits for CCATCCT" << endl;
    } else if (index_type == "hash") {
        HashIndex h = HashIndex(reads_filename, query_length, true);
    } else if (index_type == "gk") {
        char* f = const_cast<char*>(reads_filename.c_str());
        gkarrays::gkArrays *reads = new gkarrays::gkArrays(f, query_length,
                true, 0, true, 4);
    } else {
        cerr << "Unknown index_type '" << index_type << "'" << endl;
        return false;
    }
    t2 = std::chrono::system_clock::now();
    elapsed = t2 - t1;
    cout << fixed << "Construction took " << elapsed.count() << "s" << endl;
    cout << "Referenced memory: " << get_referenced_memory_size() << "kB" << endl;

    return true;
}

int main(int argc, char** argv) {
    if (argc < 5) {
        cerr << "Construct index " << endl;
        cerr << "Usage: " << argv[0] << " <index_type> <filename> ";
        cerr << "<read_length> <query_length>" << endl;
        cerr << "index_type: cr, hash, gk" << endl;
        cerr << "Example: " << endl;
        cerr << "  " << argv[0] << "cr bacteria.fastq 100 13" << endl;
        exit(1);
    }

    try {
        if (test(argv[1], argv[2], stoi(argv[3]), stoi(argv[4]))) {
            cout << "OK." << endl;
            return 0;
        } else {
            cerr << "Failed to pass." << endl;
            exit(1);
        }
    } catch (exception &e) {
        cerr << "Error: " << e.what() << endl;
        exit(1);
    }
}
