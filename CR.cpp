#include "CR.hpp"

using namespace std;

const bool CR::DEFAULT_VERBOSITY = false;
const int CR::DEFAULT_READ_LENGTH = 100;

CR::CR(string p, int rl, bool v) {
    this->positions = vector<pair<int, int>>();
    this->read_length = rl;
    verbose = v;

    string superstring;

    tie(superstring, this->positions) = preprocess(p);
    this->fm_index = FMWrapper(superstring);
    sort(this->positions.begin(), this->positions.end());
    superstring.clear();
}

CR::CR(string superstring, vector<pair<int, int>> p, int rl, bool v) {
    this->positions = p;
    this->read_length = rl;
    verbose = v;

    this->fm_index = FMWrapper(superstring);
    sort(this->positions.begin(), this->positions.end());
    superstring.clear();
}

pair<string, vector<pair<int, int>>> CR::preprocess(string p) {
    vector<pair<int, int>> _positions = vector<pair<int, int>>();

    boost::filesystem::path genome_path = boost::filesystem::path(p);

    ifstream f(genome_path.string().c_str());

    if (!f) {
        throw runtime_error("Error opening file '" +
                genome_path.string()+ "'");
    }

    // create temporary directory
    boost::filesystem::path tmpdir = boost::filesystem::temp_directory_path();
    srand (time(NULL));
    tmpdir += "/" + to_string(rand() % 1000000) + "/";
    string prefix = to_string(rand() % 1000000);
    debug("creating temporary directory " + tmpdir.string() + "\n");
    if (!boost::filesystem::create_directory(tmpdir)) {
        throw runtime_error("failed to create temporary directory");
    }

    // sga index
    string cmd = "sga index -v -a ropebwt -c -t 4 -p " + tmpdir.string() +
            "/" + prefix + " " + genome_path.string();
    redi::ipstream ips(cmd);
    string s;
    while (getline(ips, s)) {
        debug(s);
    }
    ips.close();
    int exit_code = ips.rdbuf()->status();
    if (exit_code != 0) {
        throw runtime_error("Error while executing '" + cmd +
                "'. Exit code: " + to_string(exit_code));
    }

    // sga overlap
    string cmd2 = "sga overlap -v -t 4 -p " + tmpdir.string() + "/" + prefix +
            " " + genome_path.string();
    redi::ipstream ips2(cmd2);
    string s2;
    while (getline(ips2, s2)) {
        debug(s2);
    }
    ips2.close();
    int exit_code2 = ips2.rdbuf()->status();
    if (exit_code2 != 0) {
        throw runtime_error("Error while executing '" + cmd2 +
                "'. Exit code: " + to_string(exit_code));
    }

    boost::filesystem::path overlap_path = boost::filesystem::path(
            boost::filesystem::current_path() / (genome_path.stem().string() + ".asqg.gz"));
    if (!boost::filesystem::exists(overlap_path)) {
        throw runtime_error("Error. Overlap file '" + overlap_path.string() +
                "' not found");
    }

    // move overlap file to tmpdir
    info(overlap_path.string() + " created.");
    boost::filesystem::path aux = boost::filesystem::path(tmpdir / (prefix +
            ".asqg.gz"));
    boost::filesystem::rename(overlap_path, aux);
    overlap_path = aux;

    // sga assemble
    string cmd3 = "sga assemble -v -o " + tmpdir.string() +
            "/" + prefix + " " + overlap_path.string();
    redi::ipstream ips3(cmd3);
    string s3;
    while (getline(ips3, s3)) {
        debug(s3);
    }
    ips3.close();
    int exit_code3 = ips3.rdbuf()->status();
    if (exit_code3 != 0) {
        throw runtime_error("Error while executing '" + cmd3 +
                "'. Exit code: " + to_string(exit_code3));
    }

    boost::filesystem::path contigs_path = boost::filesystem::path(
            tmpdir / (prefix + "-contigs.fa"));
    if (!boost::filesystem::exists(contigs_path)) {
        throw runtime_error("Error. Contigs file '" + overlap_path.string() +
                "' not found");
    }

    info(contigs_path.string() + " created.");

    // read contigs
    string superstring = cr_util::load_contigs(contigs_path.string());

    // construct FM-index, find each read in it and add missing reads
    FMWrapper _fm_index = FMWrapper(superstring);

    int missing_read_count = 0;
    int total_reads_size = 0;

    debug("Querying FM index for reads");
    string l;
    int read_count = 0;
    int line_count = 0;
    int skipped = 0;
    while (getline(f, l)) {
        if (line_count++ % 4 != 1) {
            continue;
        }
        boost::algorithm::trim(l);
        if (!cr_util::check_read(l)) {
            skipped++;
            continue;
        }

        total_reads_size += l.length();

        if (read_count % 10000 == 0) {
            debug("Processed " + to_string(read_count) + " reads");
        }

        vector<int> matches = _fm_index.locate(l);
        vector<int> matches2 = _fm_index.locate(cr_util::rev_compl(l));

        if (matches.size() == 0 && matches2.size() == 0) {
            missing_read_count++;
            _positions.push_back(pair<int, int>(superstring.length(),
                    read_count));
            superstring += l;
        } else {
            for (int m : matches) {
                _positions.push_back(pair<int, int>(m, read_count));
            }
            for (int m : matches2) {
                _positions.push_back(pair<int, int>(m, read_count));
            }
        }

        read_count++;
    }

    debug("Missing " + to_string(missing_read_count) + " reads");
    debug("Positions size: " + to_string(_positions.size()));
    info("Total reads size: " + to_string(total_reads_size) + "\n" +
            "Superstring size:  " + to_string(superstring.size()) + "\n" +
            "Compress ratio: " + to_string((float)total_reads_size / (float)superstring.size()));

    return pair<string, vector<pair<int, int>>>(superstring, _positions);
}

vector<pair<int, int>> CR::locate_positions(string s) {
    vector<pair<int, int>> retval;
    vector<int> indexes = this->fm_index.locate(s);

    for (auto i : indexes) {
        pair<int, int> start_index(i + s.length() - this->read_length, -1);
        pair<int, int> end_index(i, numeric_limits<int>::max());
        auto low = lower_bound(this->positions.begin(), this->positions.end(), start_index);
        auto up = upper_bound(this->positions.begin(), this->positions.end(), end_index);
        for (auto it = low; it != up; it++) {
            retval.push_back(*it);
        }
    }

    return retval;
}

vector<int> CR::find_indexes(const string& s) {
    boost::algorithm::trim_copy(s);
    vector<int> retval;

    for (auto i : this->locate_positions(s)) {
        retval.push_back(i.second);
    }
    for (auto i : this->locate_positions(cr_util::rev_compl(s))) {
        retval.push_back(i.second);
    }

    sort(retval.begin(), retval.end());

    string debug_string = "";
    for (auto i : retval) {
        debug_string += to_string(i) + ", ";
    }
    debug(debug_string);

    return retval;
}

vector<string> CR::find_reads(const string& s) {
    boost::algorithm::trim_copy(s);
    vector<string> retval;

    for (auto i : this->locate_positions(s)) {
        retval.push_back(this->fm_index.extract(i.first, this->read_length));
    }
    for (auto i : this->locate_positions(cr_util::rev_compl(s))) {
        retval.push_back(this->fm_index.extract(i.first, this->read_length));
    }

    //sort(retval.begin(), retval.end());

    string debug_string = "";
    for (auto i : retval) {
        debug_string += i + '\n';
    }
    debug(debug_string);

    return retval;
}

void CR::debug(string msg) {
    if (verbose) {
        cout << msg << endl;
    }
    return;
}

void CR::info(string msg) {
    cout << msg << endl;
    return;
}

CR::~CR() {
    //
}
