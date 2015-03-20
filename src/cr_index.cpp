#include "cr_index.hpp"

using namespace std;

const bool CRIndex::DEFAULT_VERBOSITY = false;
const int CRIndex::DEFAULT_READ_LENGTH = 100;
bool CRIndex::verbose = CRIndex::DEFAULT_VERBOSITY;

CRIndex::CRIndex(string p, int rl, bool v, bool crappy) {
    this->positions = vector<t_pos>();
    this->diff = vector<t_diff>();
    this->read_length = rl;
    CRIndex::verbose = v;

    debug("Constructing CRIndex...");
    string superstring;

    if (crappy) {
      tie(superstring, this->positions, this->diff) = crappy_preprocess(p, v);
    } else {
      tie(superstring, this->positions, this->diff) = preprocess(p, v);
    }
    this->fm_index = FMWrapper(superstring);
    debug("FM index done\n");
    sort(this->positions.begin(), this->positions.end());
    superstring.clear();
}

CRIndex::CRIndex(string superstring, vector<t_pos> p, vector<t_diff> d, int rl, bool v) {
    this->positions = p;
    this->diff = d;
    this->read_length = rl;
    CRIndex::verbose = v;

    this->fm_index = FMWrapper(superstring);
    sort(this->positions.begin(), this->positions.end());
    superstring.clear();
}

void CRIndex::SaveToFile(string filename) {
    ofstream of(filename);
    of << read_length << endl;
    of << fm_index.extract(0, fm_index.index_size()-1) << endl;
    of << positions.size() << endl;
    for (auto &e: positions) {
        of << get<0>(e) << " " << get<1>(e) << " " << get<2>(e) << endl;
    }
    of << diff.size() << endl;
    for (auto &e: diff) {
        of << get<0>(e) << " " << get<1>(e) << " " << get<2>(e) << endl;
    }
    of.close();
}

CRIndex* CRIndex::LoadFromFile(string filename) {
    ifstream f(filename);
    CRIndex* index = new CRIndex();
    f >> index->read_length;
    string superstring;
    f >> superstring;
    index->fm_index = FMWrapper(superstring);
    superstring.clear();
    int positions_length;
    f >> positions_length;
    index->positions.resize(positions_length);
    for (int i = 0; i < positions_length; i++) {
        f >> get<0>(index->positions[i]);
        f >> get<1>(index->positions[i]);
        f >> get<2>(index->positions[i]);
    }
    int diff_length;
    f >> diff_length;
    index->diff.resize(diff_length);
    for (int i = 0; i < diff_length; i++) {
        f >> get<0>(index->diff[i]);
        f >> get<1>(index->diff[i]);
        f >> get<2>(index->diff[i]);
    }
    f.close();
    return index;
}

tuple<string, vector<t_pos>, vector<t_diff>> CRIndex::crappy_preprocess(string p, bool v) {
    CRIndex::verbose = v;
    string superstring;
    vector<t_pos> _positions = vector<t_pos>();
    vector<t_diff> _diff = vector<t_diff>();
    string read_label, read, read_meta, read_q;
    ifstream reads_istream(p);
    int total_reads_size = 0;
    int read_count = 0;
    while (getline(reads_istream, read_label)) {
        getline(reads_istream, read);
        getline(reads_istream, read_meta);
        getline(reads_istream, read_q);

        boost::algorithm::trim(read);
        boost::algorithm::trim(read_label);

        int read_id = read_count;

        total_reads_size += read.length();

        if (read_count % 10000 == 0) {
            debug("Processed " + to_string(read_count) + " reads");
        }

        _positions.push_back(make_tuple(superstring.length(), read_id, 0));
        superstring += read;

        read_count++;
    }
    reads_istream.close();
    debug("Positions size: " + to_string(_positions.size()));
    info("Total reads size: " + to_string(total_reads_size) + "\n" +
            "Superstring size:  " + to_string(superstring.size()) + "\n" +
            "Compress ratio: " + to_string((float)total_reads_size / (float)superstring.size()));
    return make_tuple(superstring, _positions, _diff);
}

tuple<string, vector<t_pos>, vector<t_diff>> CRIndex::preprocess(string p, bool v) {
    CRIndex::verbose = v;
    vector<t_pos> _positions = vector<t_pos>();
    vector<t_diff> _diff = vector<t_diff>();

    boost::filesystem::path orig_reads_path = boost::filesystem::path(p);

    // create temporary directory
    boost::filesystem::path tmpdir = cr_util::create_tmpdir();
    debug("Temporary directory " + tmpdir.string() + " created");
    boost::filesystem::path p_tmpdir = tmpdir / "cr"; // e.g. /tmp/634634/cr

    debug("Running sga index");
    debug(cr_util::execute_command("sga index -a ropebwt -c -t 4 -p " +
            p_tmpdir.string() + " " + orig_reads_path.string()));
    debug("Running sga correct");
    debug(cr_util::execute_command("sga correct -t 4 -p " +
                p_tmpdir.string() + " -o " + p_tmpdir.string() + ".ec.fa " +
                orig_reads_path.string()));
    boost::filesystem::path corr_reads_path = boost::filesystem::path(
            p_tmpdir.string() + ".ec.fa");
    cr_util::check_path_existence(corr_reads_path);

    boost::filesystem::path corr_ncrit_reads_path = boost::filesystem::path(
                p_tmpdir.string() + ".ncrit.corr");
    boost::filesystem::path orig_crit_reads_path = boost::filesystem::path(
                p_tmpdir.string() + ".crit");

    ifstream orig_reads_istream(orig_reads_path.string());
    ifstream corr_reads_istream(corr_reads_path.string());
    ofstream corr_ncrit_reads_ostream(corr_ncrit_reads_path.string());
    ofstream orig_crit_reads_ostream(orig_crit_reads_path.string());

    string orig_read;
    string corr_read;
    string read_label;
    string read_meta;
    string read_q;
    string ll;
    int crit_count = 0;
    int ncrit_count = 0;
    int read_count = 0;
    debug("Searching for critical reads");
    while (getline(orig_reads_istream, read_label)) {
        getline(corr_reads_istream, ll);
        getline(corr_reads_istream, corr_read);
        getline(corr_reads_istream, ll);
        getline(corr_reads_istream, ll);

        getline(orig_reads_istream, orig_read);
        getline(orig_reads_istream, read_meta);
        getline(orig_reads_istream, read_q);

        boost::algorithm::trim(orig_read);
        boost::algorithm::trim(corr_read);
        vector<int> diff_indexes = cr_util::diff_indexes(orig_read, corr_read);
        int read_id = read_count;
        if (diff_indexes.size() >= 2 && cr_util::indexes_close(diff_indexes, 15)) {
            orig_crit_reads_ostream << "@" << read_count << endl;
            orig_crit_reads_ostream << orig_read << endl;
            orig_crit_reads_ostream << read_meta << endl;
            orig_crit_reads_ostream << string(orig_read.size(), '~') << endl;
            crit_count += 1;
        } else {
            corr_ncrit_reads_ostream << "@" << read_count << endl;
            corr_ncrit_reads_ostream << corr_read << endl;
            corr_ncrit_reads_ostream << read_meta << endl;
            corr_ncrit_reads_ostream << string(corr_read.size(), '~') << endl;
            for (int i : diff_indexes) {
                _diff.push_back(make_tuple(read_id, i, orig_read[i]));
            }
            ncrit_count += 1;
        }

        read_count++;
    }
    corr_ncrit_reads_ostream.close();
    orig_crit_reads_ostream.close();
    orig_reads_istream.close();
    corr_reads_istream.close();

    debug("Total " + to_string(read_count) + " reads. critical: " +
            to_string(crit_count) + " noncritical: " + to_string(ncrit_count));

    // rerun sga index with noncritical reads only
    debug("Running sga index");
    debug(cr_util::execute_command("sga index -a ropebwt -c -t 4 -p " +
            p_tmpdir.string() + " " + corr_ncrit_reads_path.string()));
    debug("Running sga overlap");
    debug(cr_util::execute_command("sga overlap -t 4 -p " +
            p_tmpdir.string() + " " + corr_ncrit_reads_path.string()));

    // get path of overlap file because dickish sga overlap saves it to workdir
    boost::filesystem::path overlap_path = boost::filesystem::path(
            boost::filesystem::current_path() / (corr_ncrit_reads_path.stem().string() + ".asqg.gz"));
    cr_util::check_path_existence(overlap_path);

    // move overlap file to tmpdir
    boost::filesystem::path aux = boost::filesystem::path(p_tmpdir.string() + ".asqg.gz");
    boost::filesystem::rename(overlap_path, aux);
    overlap_path = aux;

    debug("Running sga assemble");
    debug(cr_util::execute_command("sga assemble -o " + p_tmpdir.string() +
            " " + overlap_path.string()));

    boost::filesystem::path contigs_path = boost::filesystem::path(
            p_tmpdir.string() + "-contigs.fa");
    cr_util::check_path_existence(contigs_path);

    // read contigs
    string superstring = cr_util::load_contigs(contigs_path.string());

    // construct FM-index, find each read in it and add missing reads
    FMWrapper _fm_index = FMWrapper(superstring);

    debug("Superstring size: " + to_string(superstring.size()) + " (before" +
            " querying)");

    debug("Querying FM index for noncritical corrected reads");
    int missing_read_count = 0;
    int total_reads_size = 0;
    string read;
    read_count = 0;
    int skipped = 0;
    ifstream corr_ncrit_reads_istream(corr_ncrit_reads_path.string());
    while (getline(corr_ncrit_reads_istream, read_label)) {
        getline(corr_ncrit_reads_istream, read);
        getline(corr_ncrit_reads_istream, read_meta);
        getline(corr_ncrit_reads_istream, read_q);

        boost::algorithm::trim(read);
        boost::algorithm::trim(read_label);
        if (!cr_util::check_read(read)) {
            skipped++;
            continue;
        }

        int read_id = stoi(read_label.substr(1));

        total_reads_size += read.length();

        if (read_count % 10000 == 0) {
            debug("Processed " + to_string(read_count) + " reads");
        }

        vector<int> matches = _fm_index.locate(read);
        vector<int> matches2 = _fm_index.locate(cr_util::rev_compl(read));

        if (matches.size() == 0 && matches2.size() == 0) {
            missing_read_count++;
            _positions.push_back(make_tuple(superstring.length(), read_id, 0));
            superstring += read;
        } else {
            for (int m : matches) {
                _positions.push_back(make_tuple(m, read_id, 0));
            }
            for (int m : matches2) {
                _positions.push_back(make_tuple(m, read_id, 1));
            }
        }

        read_count++;
    }
    corr_ncrit_reads_istream.close();

    debug("Querying FM index for original critical reads");
    ifstream orig_crit_reads_istream(orig_crit_reads_path.string());
    while (getline(orig_crit_reads_istream, read_label)) {
        getline(orig_crit_reads_istream, read);
        getline(orig_crit_reads_istream, read_meta);
        getline(orig_crit_reads_istream, read_q);

        boost::algorithm::trim(read);
        boost::algorithm::trim(read_label);
        if (!cr_util::check_read(read)) {
            skipped++;
            continue;
        }

        int read_id = stoi(read_label.substr(1));

        total_reads_size += read.length();

        if (read_count % 10000 == 0) {
            debug("Processed " + to_string(read_count) + " reads");
        }

        _positions.push_back(make_tuple(superstring.length(), read_id, 0));
        superstring += read;

        read_count++;
    }
    orig_crit_reads_istream.close();
    debug("Superstring:" + superstring.substr(0, 20) + "..." + 
          superstring.substr(superstring.size()-20));

    if (skipped > 0) {
        info("WARNING: skipped " + to_string(skipped) + " reads");
    }
    debug("Processed " + to_string(read_count) + " reads");
    debug("Missing " + to_string(missing_read_count) + " reads (noncritical " +
            "dropped by sga assembler + critical)");
    debug("Positions size: " + to_string(_positions.size()));
    info("Total reads size: " + to_string(total_reads_size) + "\n" +
            "Superstring size:  " + to_string(superstring.size()) + "\n" +
            "Compress ratio: " + to_string((float)total_reads_size / (float)superstring.size()) + "\n" +
            "Diff size: " + to_string(_diff.size()));

    return make_tuple(superstring, _positions, _diff);
}

vector<t_pos> CRIndex::locate_positions2(const string& s, const string& s_check) {
    vector<t_pos> retval;
    vector<int> indexes = this->fm_index.locate(s);

    for (auto i : indexes) {
        t_pos start_index(i + s.length() - this->read_length, -1, 0);
        t_pos end_index(i, numeric_limits<int>::max(), 1);
        auto low = lower_bound(this->positions.begin(), this->positions.end(), start_index);
        auto up = upper_bound(this->positions.begin(), this->positions.end(), end_index);
        for (auto it = low; it != up; it++) {
            // drop false positives
            int pos = get<0>(*it);
            int read_id = get<1>(*it);
            int s_pos_in_read = i - pos;
            int s_pos_in_rev_compl_read = abs(i - pos + (int)s.size() - (int)this->read_length);
            bool rev_compl = get<2>(*it);

            t_diff start_index2(read_id, -1, 'A');
            t_diff end_index2(read_id, numeric_limits<int>::max(), 'A');
            auto low2 = lower_bound(this->diff.begin(), this->diff.end(), start_index2);
            auto up2 = upper_bound(this->diff.begin(), this->diff.end(), end_index2);

            string s2 = s;
            if (rev_compl) {
                s2 = cr_util::rev_compl(s2);

                for (auto it2 = low2; it2 != up2; it2++) {
                    int j = get<1>(*it2) - s_pos_in_rev_compl_read;
                    if (j >= 0 && (size_t)j < s2.size()) {
                        s2[j] = get<2>(*it2);
                    }
                }

                if (s2 == cr_util::rev_compl(s_check)) {
                    retval.push_back(*it);
                }
            } else {
                for (auto it2 = low2; it2 != up2; it2++) {
                    int j = get<1>(*it2) - s_pos_in_read;
                    if (j >= 0 && (size_t)j < s2.size()) {
                        s2[j] = get<2>(*it2);
                    }
                }

                if (s2 == s_check) {
                    retval.push_back(*it);
                }
            }
        }
    }

    return retval;
}

vector<t_pos> CRIndex::locate_positions(const string& s) {
    vector<t_pos> retval = locate_positions2(s, s);

    for (string s2 : cr_util::strings_with_edt1(s)) {
        vector<t_pos> temp = locate_positions2(s2, s);
        retval.insert(retval.end(), temp.begin(), temp.end());
    }

    return retval;
}

vector<int> CRIndex::find_indexes(const string& s) {
    boost::algorithm::trim_copy(s);
    vector<int> retval;

    for (auto i : this->locate_positions(s)) {
        retval.push_back(get<1>(i));
    }
    for (auto i : this->locate_positions(cr_util::rev_compl(s))) {
        retval.push_back(get<1>(i));
    }

    sort(retval.begin(), retval.end());
    retval.erase(unique(retval.begin(), retval.end()), retval.end());

    debug("Found " + to_string(retval.size()) + " occurrences.");
    debug(retval);

    return retval;
}

vector<string> CRIndex::find_reads(const string& s) {
    boost::algorithm::trim_copy(s);
    vector<string> retval;

    for (auto i : this->locate_positions(s)) {
        retval.push_back(extract_original_read(i));
    }
    for (auto i : this->locate_positions(cr_util::rev_compl(s))) {
        retval.push_back(extract_original_read(i));
    }

    sort(retval.begin(), retval.end());
    retval.erase(unique(retval.begin(), retval.end()), retval.end());

    debug("Found " + to_string(retval.size()) + " occurrences.");
    debug(retval);

    return retval;
}

string CRIndex::extract_original_read(t_pos read_p) {
    int pos = get<0>(read_p);
    int read_id = get<1>(read_p);
    bool rev_compl = get<2>(read_p);

    t_diff start_index(read_id, -1, 'A');
    t_diff end_index(read_id, numeric_limits<int>::max(), 'A');
    auto low = lower_bound(this->diff.begin(), this->diff.end(), start_index);
    auto up = upper_bound(this->diff.begin(), this->diff.end(), end_index);

    string retval = this->fm_index.extract(pos, this->read_length);
    if (rev_compl) {
        retval = cr_util::rev_compl(retval);
    }
    for (auto it = low; it != up; it++) {
        retval[get<1>(*it)] = get<2>(*it);
    }

    return retval;
}

void CRIndex::debug(string msg) {
    if (!CRIndex::verbose) {
        return;
    }

    cout << msg << endl;

    return;
}

void CRIndex::debug(vector<string> msg) {
    if (!CRIndex::verbose) {
        return;
    }

    string debug_string = "";
    for (auto i : msg) {
        debug_string += i + '\n';
    }
    cout << debug_string << endl;

    return;
}

void CRIndex::debug(vector<int> msg) {
    if (!CRIndex::verbose) {
        return;
    }

    string debug_string = "";
    for (auto i : msg) {
        debug_string += to_string(i) + ", ";
    }
    cout << debug_string << endl;

    return;
}


void CRIndex::info(string msg) {
    cout << msg << endl;
    return;
}

CRIndex::~CRIndex() {
    //
}
