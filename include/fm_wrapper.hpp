#include <vector>
#include <string>
#include <sdsl/suffix_arrays.hpp>

using namespace sdsl;
using namespace std;

//typedef csa_wt<wt_rlmn<bit_vector, rank_support_v5<>, select_support_scan<>,
//               wt_huff<bit_vector, rank_support_v5<>, select_support_scan<>, select_support_scan<>>>,
//typedef csa_wt<wt_huff<bit_vector, rank_support_v5<>, select_support_scan<>, select_support_scan<>>,
//        16, 32,
//        text_order_sa_sampling<sd_vector<>>> fm_index_type;
//typedef csa_wt<wt_huff<rrr_vector<127> >, 16, 32> fm_index_type;

template<uint32_t n1=8, uint32_t n2=5000000>
class FMWrapper {
    typedef csa_wt<wt_huff<bit_vector>, n1, n2> fm_index_type;
    public:
        FMWrapper() {}
        FMWrapper(const string& data){
          fm_index_type fm;
          sdsl::construct_im(fm, data.c_str(), 1);
          this->fm_index = fm;
        }
 
        vector<int> locate(const string& query) {
          auto retval = sdsl::locate(this->fm_index, query.begin(), query.end());
          return vector<int>(retval.begin(), retval.end());
        }

        string extract(int start, int length) {
          return sdsl::extract(this->fm_index, start, start + length - 1);
        }

        int memory_size() {
          return sdsl::size_in_bytes(this->fm_index);
        }

        size_t index_size() {
          return fm_index.size();
        }

    private:
         fm_index_type fm_index;
};
