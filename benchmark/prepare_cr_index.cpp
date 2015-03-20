#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <iostream>
#include <vector>
#include <string>
#include <utility>
#include <chrono>
#include <boost/algorithm/string.hpp>
#include "cr_index.hpp"

using namespace std;

void go(string genome_filename, string index_filename) {
    ifstream f(genome_filename);
    string l;
    getline(f, l);
    getline(f, l);
    f.close();
    int read_length = l.size();

    cout << "Building CRIndex ... " << endl;
    CRIndex rm = CRIndex(genome_filename, read_length, true, true);
    rm.SaveToFile(index_filename);
}

int main(int argc, char** argv) {
    if (argc < 3) {
        cerr << "Usage: " << argv[0] << " <filename> <filename>" << endl;
        cerr << "Examples: " << endl;
        cerr << "  Read .fastq file and construct CR-index and save to file" << endl;
        cerr << "  " << argv[0] << " bacteria.fastq index.data" << endl;
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
