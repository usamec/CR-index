#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <queue>
#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>

using namespace std;

string rev_compl(const string& s) {
    string rc = "";
    for (size_t i = 0; i < s.length(); i++) {
        switch (s[i]) {
        case 'A':
            rc += 'T';
            break;
        case 'T':
            rc += 'A';
            break;
        case 'C':
            rc += 'G';
            break;
        case 'G':
            rc += 'C';
            break;
        default:
            rc += s[i];
            break;
        }
    }
    reverse(rc.begin(), rc.end());
    return rc;
}

void Erase(string &x, unordered_map<string, unordered_set<string>>& gg) {
  for (auto &y: gg[x]) {
    gg[y].erase(x);
  }
  gg.erase(x);
}

int GetBestOverL(const string &a, const string &b) {
  string rc = rev_compl(a);
  int start = min(a.size(), b.size());
  start = min(151, start);
  for (int i = start; i > 0; i--) {
    if (a.substr(a.size() - i, i) == b.substr(0, i)) {
      return i;
    }
    if (b.substr(b.size() - i, i) == a.substr(0, i)) {
      return i;
    }
    if (rc.substr(rc.size() - i, i) == b.substr(0, i)) {
      return i;
    }
    if (b.substr(b.size() - i, i) == rc.substr(0, i)) {
      return i;
    }
  }
  return 0;
}

string GetBestOver(const string &a, const string &b) {
  string rc = rev_compl(a);
  int start = min(a.size(), b.size());
  start = min(151, start);
  for (int i = start; i > 0; i--) {
    if (a.substr(a.size() - i, i) == b.substr(0, i)) {
      return a + b.substr(i);
    }
    if (b.substr(b.size() - i, i) == a.substr(0, i)) {
      return b + a.substr(i);
    }
    if (rc.substr(rc.size() - i, i) == b.substr(0, i)) {
      return rc + b.substr(i);
    }
    if (b.substr(b.size() - i, i) == rc.substr(0, i)) {
      return b + rc.substr(i);
    }
  }
  return a + b;
}

string DoOverlap(vector<string>& cur_g, vector<string>& cur_c, unordered_map<string, unordered_set<string>>& g) {
  unordered_map<string, unordered_set<string>> gg;
  unordered_map<string, string> rr;
  for (int i = 0; i < cur_c.size(); i++) {
    rr[cur_g[i]] = cur_c[i];
    gg[cur_g[i]] = g[cur_g[i]];
  }

  string ret;
  while (gg.size() > 0) {
    string start = "";
    string cur = "";
    for (auto &x: gg) {
      if (x.second.size() <= 1) {
        start = x.first;
        break;
      }
    }
    if (start == "") {
      start = gg.begin()->first;
    }
    cur = rr[start];
    while (true) {
      string next = "";
      int bo = 0;
      for (auto &y: gg[start]) {
        int l = GetBestOverL(cur, rr[y]);
        if (l > bo) {
          next = y;
          bo = l;
        }
      }
      cur = GetBestOver(rr[next], cur);
      Erase(start, gg);
      if (next == "") break;
      start = next;
    }
    ret += cur;
  }

  for (int i = 0; i < cur_c.size(); i++) {
    if (ret.find(cur_c[i]) == string::npos &&
        ret.find(rev_compl(cur_c[i])) == string::npos) {
      printf("Dafug\n");
      exit(457);
    }
  }

  return ret;
}


int main(int argc, char** argv) {
  FILE *outf = fopen(argv[2], "w");
  ifstream file(argv[1], ios_base::in | ios_base::binary);
  boost::iostreams::filtering_streambuf<boost::iostreams::input> in;
  in.push(boost::iostreams::gzip_decompressor());
  in.push(file);
  istream input(&in);
  string header;
  getline(input, header);
  cout << header << endl;
  string line;
  int total_len = 0;
  int over_len = 0;
  unordered_map<string, vector<int>> g_lens;
  unordered_map<string, unordered_set<string>> g;
  unordered_map<string, string> reads;
  vector<pair<pair<string, string>, int>> edges;
  while (getline(input, line)) {
    if (line[0] == 'V') {
      vector<string> items;
      boost::algorithm::split(items, line, boost::algorithm::is_any_of("\t"));
      total_len += items[2].size();
      if (g.count(items[1]) == 0) {
        g[items[1]] = unordered_set<string>();
        g[items[1]] = unordered_set<string>();
      }
      reads[items[1]] = items[2];
    }
    if (line[0] == 'E') {
      vector<string> items;
      boost::algorithm::split(items, line, boost::algorithm::is_any_of("\t "));
      over_len += abs(stoi(items[3]) - stoi(items[4]));
      g[items[1]].insert(items[2]);
      g[items[2]].insert(items[1]);
      int over_len = stoi(items[4]) - stoi(items[3]);
      g_lens[items[1]].push_back(over_len);
      g_lens[items[2]].push_back(over_len);
      edges.push_back(make_pair(make_pair(items[1], items[2]), over_len));
    }
  }
  printf("graph loaded\n");
//  g = g1;
  vector<string> cur_c;
  vector<string> cur_g;
  int before = 0;
  int after = 0;
  int total_reads = 0;
  fprintf(outf, ">r1\n");
  unordered_set<string> visited;
  for (auto &x: g) {
    cur_c.clear();
    cur_g.clear();
    if (visited.count(x.first)) continue;
    int c_size = 0;
    queue<string> fr;
    fr.push(x.first);
    visited.insert(x.first);
    while (!fr.empty()) {
      c_size++;
      string x = fr.front(); fr.pop();
      cur_c.push_back(reads[x]);
      cur_g.push_back(x);
      for (auto &y: g[x]) {
        if (visited.count(y)) continue;
        visited.insert(y);
        fr.push(y);
      }
    }
    for (auto &x: cur_c) {
      before += x.size();
    }
    string over = DoOverlap(cur_g, cur_c, g);
    fprintf(outf, "%s", over.c_str());
    after += over.size();
    total_reads += cur_c.size();
  }
  fprintf(outf, "\n");
  fprintf(stdout, "before %d after %d reads %d\n", before, after, total_reads);
}
