#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include <iomanip>
#include <chrono>

#include <sdsl/bit_vectors.hpp>

#include <zombit/zombit_vector.hpp>


std::string g_vec_path = "";
std::string g_block_model = "";
uint64_t g_postings_list_size = 0;
sdsl::bit_vector g_bv;


template <typename T>
void read_input_file(const char *input_file, std::vector<T> &vec_ref) {
    FILE *fp = fopen(input_file, "r");
    uint64_t u, n;

    fseek(fp, 0, SEEK_END);
    n = ftell(fp) / sizeof(T);
    fseek(fp, 0, SEEK_SET);

    std::cerr << ">> Number of integers: " << n << std::endl;
    vec_ref.resize(n);
    std::cerr << ">> Reading input...";
    for (uint64_t i = 0; i < n; i++) {
        fread(&u, sizeof(T), 1, fp);
        vec_ref[i] = u;
    }


    fclose(fp);  // closing inputfile
    std::cerr << " DONE" << std::endl;

}


// creating new vector v from w -> v[i] = w[i] - w[i-1]
void create_diff_vector(const std::vector<uint64_t> &seq_ref2, std::vector<uint64_t> &seq_diff_ref, uint64_t &sumOrig_r) {
    std::cerr << ">> Creating difference vector...";
    for (uint64_t i = 1; i < seq_ref2.size(); i++) {
        if (seq_ref2[i] > seq_ref2[i - 1]) {
            seq_diff_ref[i] = seq_ref2[i] - seq_ref2[i - 1];
            sumOrig_r += (seq_ref2[i] - seq_ref2[i - 1]);
        } else {
            seq_diff_ref[i] = seq_ref2[i];
            sumOrig_r += seq_ref2[i];
        }
    }
    std::cerr << " DONE" << std::endl;
}

void create_bit_vector(sdsl::bit_vector &bv, const std::vector<uint64_t> &vec_ref, const uint64_t &n) {
    bv = sdsl::bit_vector(n + 1);

    for (uint64_t i = 0; i < n; i++) bv[i] = 0;

    uint64_t sum = vec_ref[0];
    bv[sum] = 1;
    for (uint64_t i = 1; i < vec_ref.size(); i++) {
        sum += vec_ref[i];
        bv[sum] = 1;
    }
}

template<class T>
void test_zombit(const uint32_t b,
    const std::string label,
    const std::vector<uint64_t> &r) {

  int32_t max_rec_depth = (int32_t) std::log2(b);
  float prev_bpp = 1000000.0;
  for (int32_t rec_level = 0; rec_level < max_rec_depth; rec_level++) {
    T zombit{};
    std::cerr << ">> building zombit" << label << " b:" << b << " rec_level:" << rec_level << "...";
    zombit.build_zombit(g_bv,b, rec_level, g_block_model);
    std::cerr << "DONE";

    std::cerr << " benchmarking...";
    std::chrono::duration<double, std::micro> ms_double;

    float bpp;
    double time_avg = 0;
    double time_avg_total = 0;

    // rank scan testing
    for (uint32_t j = 0; j < 3; j++) {
      time_avg = 0;
      for (uint64_t i = 0; i < r.size(); i++) {
        auto t1 = std::chrono::high_resolution_clock::now();
        zombit.nextGEQ(r[i]);
        auto t2 = std::chrono::high_resolution_clock::now();
        ms_double = t2 - t1;
        time_avg += ms_double.count();
      }
      time_avg_total += (time_avg / r.size());
    }
    bpp = (float) zombit.size_in_bits() / g_postings_list_size;
    std::cout << label << ";" << b << ";" << bpp;
    std::cout << ";" << ((float)zombit.u_vector_size_in_bits() / g_postings_list_size);
    std::cout << ";" << ((float)zombit.o_vector_size_in_bits() / g_postings_list_size);
    std::cout << ";" << ((float)zombit.m_vector_size_in_bits() / g_postings_list_size);
    std::cout << ";" << ((float)zombit.u_vector_size_in_bits() / zombit.size_in_bits());
    std::cout << ";" << ((float)zombit.o_vector_size_in_bits() / zombit.size_in_bits());
    std::cout << ";" << ((float)zombit.m_vector_size_in_bits() / zombit.size_in_bits());
    std::cout << ";" << zombit.block_n << ";" << zombit.m_blocks << ";" << zombit.runs_n << ";" << rec_level;
    std::cout << ";" << (time_avg_total/3) << ";rank_scan" << "\n";
    std::cerr << "DONE\n";

    if (prev_bpp <= bpp) break;
    prev_bpp = bpp;
  }
}


int main(int argc, char *argv[]) {

	if (argc < 5) {
        std::cout << ">> Program usage: " << argv[0] << " <input_file> <mode> <zombit_models> <vec_path> <rec_model> \n";
        exit(1);
    }

    std::cerr << ">> Program START\n";
    std::string mode = argv[2];
    if (mode == "raw") {

        std::vector<uint64_t> seq;
        read_input_file<uint64_t>(argv[1], seq);
        g_postings_list_size = seq.size();
        std::cerr << ">> inputfile size: " << g_postings_list_size << "\n";

        std::vector<uint64_t> seq_diff(seq.size());
        seq_diff[0] = seq[0];
        uint64_t sumOrig = seq[0];
        create_diff_vector(seq, seq_diff, sumOrig);
        seq.clear();
        seq.shrink_to_fit();

        // creating bit_vector
        std::cerr << ">> creating bitvector of size: " << (sumOrig+1) << "...";
        create_bit_vector(g_bv, seq_diff, sumOrig);
        std::cerr << "Done\n";
        seq_diff.clear();
        seq_diff.shrink_to_fit();
    } else if (mode == "notraw") {
        sdsl::load_from_file(g_bv, argv[1]);
        g_postings_list_size = sdsl::rank_support_v5<1>(&g_bv)(g_bv.size());
    } else {
      std::cerr << "unknown model\n";
      exit(1);
    }

    // setting global variables
    std::string run_mode = argv[3];
    g_vec_path = argv[4];
    g_block_model = argv[5];

    // reading block lengths for test from user input
    uint32_t b_n;
    std::cin >> b_n;
    std::vector<uint32_t> block_sizes(b_n);
    for (int i = 0; i < b_n; i++) std::cin >> block_sizes[i];

    // creating test query values
    uint64_t n = g_bv.size() > 1000000 ? 1000000 : g_bv.size();
    uint64_t u = g_bv.size()-2;
    std::vector<uint64_t> benchmark_quesries(n);
    srand(0);
    for (uint32_t i = 0; i < n; i++) benchmark_quesries[i] = rand() % u + 1;

    std::cerr << ">> start testing\n";
    std::cout << "zombit<U,O,M>;block size;overall size;U size;O size;M size;U%;O%;M%;number of blocks;mixed blocks;runs of 1s;recursio level;nextGEQ avg (μs);query type\n";

    ////////////////////////
    //
    // ACTUAL ZOMBIT TESTING
    //
    ////////////////////////



    // bit_vector, bit_vector, bit_vector
    if (run_mode[0] == '1') {
      for (uint32_t b : block_sizes ) {
        test_zombit<zombit_bv_bv_bv_O2>(b, "<bv,bv,bv>", benchmark_quesries);
      }
    }

    ////  bv,bv, rrr
    //if (run_mode[1] == '1') {
    //  for (uint32_t b : block_sizes ) {
    //    test_zombit<zombit_bv_bv_rrr>(b, "<bv,bv,rrr>", benchmark_quesries);
    //  }
    //}

    //if (run_mode[2] == '1') {
    //  for (uint32_t b : block_sizes ) {
    //    zombit_hyb_rrr_rrr zombit{};
    //    test_zombit(zombit, b, "<hyb,rrr,rrr>", benchmark_quesries);
    //  }
    //}

    ////  bv,bv,sd
    //if (run_mode[3] == '1') {
    //  for (uint32_t b : block_sizes ) {
    //    test_zombit<zombit_bv_bv_sd>(b, "<bv,bv,sd>", benchmark_quesries);
    //  }
    //}

    std::cerr << ">> Program END\n";

  	return 0;
}
