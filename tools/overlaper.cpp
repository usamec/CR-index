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

string DoOverlap(vector<string> r) {
  if (r.size() >= 50) 
    fprintf(stdout, "over %d\n", r.size());
  unordered_set<string> reads(r.begin(), r.end());
  priority_queue<pair<int, pair<string, string>>> overs;
  for (auto &x: reads) {
    for (auto &y: reads) {
      if (x == y) {
        continue;
      }
      string q = GetBestOver(x, y);
      overs.push(make_pair(x.size() + y.size() - q.size(), make_pair(x, y)));
    }
  }
  while (reads.size() > 1) {
    pair<string, string> over;
    while (true) {
      over = overs.top().second;
      overs.pop();
      if (reads.count(over.first) && reads.count(over.second)) {
        break;
      }
    }
    string q = GetBestOver(over.first, over.second);
    reads.erase(over.first);
    reads.erase(over.second);
    reads.insert(q);
    for (auto &x: reads) {
      if (x == q) continue;
      string q2 = GetBestOver(x, q);
      overs.push(make_pair(x.size() + q.size() - q2.size(), make_pair(x, q)));
    }
  }
  for (int i = 0; i < r.size(); i++) {
    if (reads.begin()->find(r[i]) == string::npos &&
        reads.begin()->find(rev_compl(r[i])) == string::npos) {
      printf("Dafug\n");
      exit(457);
    }
  }
  return *reads.begin();
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
  unordered_map<string, vector<string>> g;
  unordered_map<string, string> reads;
  while (getline(input, line)) {
    if (line[0] == 'V') {
      vector<string> items;
      boost::algorithm::split(items, line, boost::algorithm::is_any_of("\t"));
      total_len += items[2].size();
      if (g.count(items[1]) == 0) {
        g[items[1]] = vector<string>();
      }
      reads[items[1]] = items[2];
    }
    if (line[0] == 'E') {
      vector<string> items;
      boost::algorithm::split(items, line, boost::algorithm::is_any_of("\t "));
      over_len += abs(stoi(items[3]) - stoi(items[4]));
      g[items[1]].push_back(items[2]);
      g[items[2]].push_back(items[1]);
    }
  }
  unordered_set<string> visited;
  map<int, int> c_dist;

  vector<string> cur_c;
  int before = 0;
  int after = 0;
  int total_reads = 0;
  fprintf(outf, ">r1\n");
  for (auto &x: g) {
    cur_c.clear();
    if (visited.count(x.first)) continue;
    int c_size = 0;
    queue<string> fr;
    fr.push(x.first);
    visited.insert(x.first);
    while (!fr.empty()) {
      c_size++;
      string x = fr.front(); fr.pop();
      cur_c.push_back(reads[x]);
      for (auto &y: g[x]) {
        if (visited.count(y)) continue;
        visited.insert(y);
        fr.push(y);
      }
    }
    c_dist[c_size]++;
    for (auto &x: cur_c) {
      before += x.size();
    }
    string over = DoOverlap(cur_c);
    fprintf(outf, "%s", over.c_str());
    after += over.size();
    total_reads += cur_c.size();
  }
  fprintf(outf, "\n");
  fprintf(stdout, "before %d after %d reads %d\n", before, after, total_reads);
  for (auto &e: c_dist) {
    cout << e.first << " " << e.second << " " << e.first * e.second << endl;
  }
}
