#ifndef _ZOMBIT_VECTOR_HPP_
#define _ZOMBIT_VECTOR_HPP_

#include <iostream>
#include <algorithm>
#include <vector>
#include <string>

#include <sdsl/bit_vectors.hpp>
#include <sdsl/rank_support_v5O2.hpp>

// basic template
template <
    typename T_u_vec, typename T_u_vec_rank,
    typename T_o_vec, typename T_o_vec_rank, typename T_o_vec_slc,
    typename T_m_vec, typename T_m_vec_rank, typename T_m_vec_slc
    >
class Zombit {
    public:

        Zombit();
        ~Zombit();
        void build_zombit(sdsl::bit_vector &X, const uint32_t b, const uint32_t recursio_level, std::string rec_model);
        void build_zombit(sdsl::bit_vector &X, const uint32_t recursio_level, std::string rec_model);
        const uint32_t access(const uint64_t i); // return value of i-th bit
        const uint64_t nextGEQ(const uint64_t x) const; // return next greater or equal than x
        //const uint64_t nextGEQ_scan(const uint64_t x) const; // return next greater or equal than x
        //const uint64_t nextGEQ_rank_scan(const uint64_t x) const; // return next greater or equal than x
        const uint64_t size(); // return lenght of the original bitvector
        const uint64_t size_in_bits();
        const uint64_t u_vector_size_in_bits();
        const uint64_t o_vector_size_in_bits();
        const uint64_t m_vector_size_in_bits();

        uint32_t block_size;
        uint64_t orig_bv_size = 0;
        uint64_t runs_n = 0;
        uint64_t m_blocks = 0;
        uint64_t block_n = 0;
        uint32_t rec_level = 0;
        std::string b_model = "";

        T_u_vec u_vector;
        T_o_vec o_vector;
        T_m_vec m_vector;
        Zombit<T_u_vec, T_u_vec_rank,\
            T_o_vec, T_o_vec_rank, T_o_vec_slc,\
            T_m_vec, T_m_vec_rank, T_m_vec_slc>* zombit_rec;

        T_u_vec_rank u_rank;
        T_o_vec_rank o_rank;
        T_o_vec_slc o_select;
        T_m_vec_rank m_rank;
        T_m_vec_slc m_select;

        size_t o_ones;
        size_t m_ones;

        const uint32_t block_model(const sdsl::bit_vector &M);
};

template <
    typename T_u_vec, typename T_u_vec_rank,
    typename T_o_vec, typename T_o_vec_rank, typename T_o_vec_slc,
    typename T_m_vec, typename T_m_vec_rank, typename T_m_vec_slc
    >
const uint32_t Zombit<T_u_vec,T_u_vec_rank,
      T_o_vec,T_o_vec_rank,T_o_vec_slc,
      T_m_vec,T_m_vec_rank,T_m_vec_slc>::block_model(const sdsl::bit_vector &M) {
    if (b_model == "div2") return block_size/2;
    if (b_model == "div4") return block_size/4;
    if (b_model == "div8") return block_size/8;
    return 1;
}

template <
    typename T_u_vec, typename T_u_vec_rank,
    typename T_o_vec, typename T_o_vec_rank, typename T_o_vec_slc,
    typename T_m_vec, typename T_m_vec_rank, typename T_m_vec_slc
    >
Zombit<T_u_vec,T_u_vec_rank,
    T_o_vec,T_o_vec_rank,T_o_vec_slc,
    T_m_vec,T_m_vec_rank,T_m_vec_slc>::Zombit() {};

template <
    typename T_u_vec, typename T_u_vec_rank,
    typename T_o_vec, typename T_o_vec_rank, typename T_o_vec_slc,
    typename T_m_vec, typename T_m_vec_rank, typename T_m_vec_slc
    >
Zombit<T_u_vec,T_u_vec_rank,
    T_o_vec,T_o_vec_rank,T_o_vec_slc,
    T_m_vec,T_m_vec_rank,T_m_vec_slc>::~Zombit() {
    if (rec_level > 0) delete zombit_rec;
}

/*
 * Note if |X| mod block_size != 0, bitvector X is modified to zero's are added end of X, such that
 * |X| mod block_size = 0 is valid.
 */
template <
    typename T_u_vec, typename T_u_vec_rank,
    typename T_o_vec, typename T_o_vec_rank, typename T_o_vec_slc,
    typename T_m_vec, typename T_m_vec_rank, typename T_m_vec_slc
    >
void Zombit<T_u_vec,T_u_vec_rank,
     T_o_vec,T_o_vec_rank,T_o_vec_slc,
     T_m_vec,T_m_vec_rank,T_m_vec_slc>::build_zombit(
        sdsl::bit_vector &X, const uint32_t b, const uint32_t recursio_level, std::string rec_model ) {

    block_size = b;
    orig_bv_size = X.size();
    rec_level = recursio_level;
    b_model = rec_model;
    if (runs_n == 0) {
        uint32_t prev = 0;
        for (uint64_t i = 0; i < X.size(); i++) {
            if (X[i] && !prev) {
                runs_n++;
                prev = 1;
            } else if (!X[i] && prev) {
                prev = 0;
            }
        }
    }
    // adding 0s end of X, so that all blocks have same length.
    if (X.size() % block_size != 0) {
        X.resize(X.size() + (block_size - (X.size() % block_size)));
    }
    block_n = X.size() / block_size;
    sdsl::bit_vector U_bv = sdsl::bit_vector(block_n);
    sdsl::bit_vector O_bv = sdsl::bit_vector(block_n);
    sdsl::bit_vector M_bv;
    char label = 'O';
    m_blocks = 0;
    std::vector<uint64_t> m_block_idx;

    // building U and O-vector
    for (uint64_t i = 0; i < block_n; i++) {
        //block first val
        if (X[i*block_size] == 0) label = 'Z';
        else label = 'O';
        for (uint64_t j = 1; j < block_size; j++) {

            if (X[i * block_size + j] == 0 && label == 'Z') continue;
            if (X[i * block_size + j] == 1 && label == 'O') continue;
            if (X[i * block_size + j] == 0 && label == 'O') {
                label = 'M';
                continue;
            }
            if (X[ (i * block_size) + j] == 1 && label == 'Z') {
                label = 'M';
                continue;
            }
        }
        // setting value for U-vector at index b_idx
        if (label == 'Z' || label == 'O') U_bv[i] = 1; // all values 0s or 1s
        else U_bv[i] = 0; // mix of 0s and 1s
        // setting value for O-vector at indexx b_idx
        if (label == 'Z') O_bv[i] = 0; // all 0s
        else O_bv[i] = 1; // block contain atleast one 1 bit

        if (label == 'M') {
            m_blocks++;
            m_block_idx.push_back(i);
        }
    }
    // building M-vector
    M_bv = sdsl::bit_vector(block_size * m_blocks);
    uint64_t m_idx = 0;
    for (uint64_t x : m_block_idx) {
        for (uint64_t k = 0; k < block_size; k++) {
            M_bv[m_idx] = X[(x * block_size) + k];
            m_idx++;
        }
    }

    m_block_idx.clear();
    m_block_idx.shrink_to_fit();

    // constructing U and O-vector and it's rank/select support
    u_vector = T_u_vec(U_bv);
    o_vector = T_o_vec(O_bv);
    u_rank = T_u_vec_rank(&u_vector);
    o_rank = T_o_vec_rank(&o_vector);
    o_select = T_o_vec_slc(&o_vector);

    if (rec_level > 0) {
        uint32_t m_b = block_model(M_bv);
        if (m_b <= 1) {
            rec_level = 0;
            rec_level = false;
        } else {
            zombit_rec = new Zombit<T_u_vec, T_u_vec_rank,\
                         T_o_vec, T_o_vec_rank, T_o_vec_slc,\
                         T_m_vec, T_m_vec_rank, T_m_vec_slc>();
            zombit_rec->build_zombit(M_bv, m_b, rec_level-1, rec_model);
        }
    }
    if (rec_level == 0) {
        m_vector = T_m_vec(M_bv);
        m_rank = T_m_vec_rank(&m_vector);
        m_select = T_m_vec_slc(&m_vector);
        m_ones = T_m_vec_rank(&m_vector)(M_bv.size());
    }

    o_ones = T_o_vec_rank(&o_vector)(O_bv.size());
}

// build zombit-vector with optimal block size sqrt(n/runs) according article
template <
    typename T_u_vec, typename T_u_vec_rank,
    typename T_o_vec, typename T_o_vec_rank, typename T_o_vec_slc,
    typename T_m_vec, typename T_m_vec_rank, typename T_m_vec_slc
    >
void Zombit<T_u_vec,T_u_vec_rank,
     T_o_vec,T_o_vec_rank,T_o_vec_slc,
     T_m_vec,T_m_vec_rank,T_m_vec_slc>::build_zombit(
        sdsl::bit_vector &X, const uint32_t recursio_level, std::string rec_model) {
    uint32_t prev = 0;
    rec_level = recursio_level;
    runs_n = 0;
    for (uint64_t i = 0; i < X.size(); i++) {
        if (X[i] && !prev) {
            runs_n++;
            prev = 1;
        } else if (!X[i] && prev) {
            prev = 0;
        }
    }
    if (runs_n == 0) runs_n = 1; // all 0s or 1s
    uint32_t b = (uint32_t) sqrt(X.size()/runs_n);
    //if (b == 1) b++;

    build_zombit(X, b, recursio_level, rec_model);
}


template <
    typename T_u_vec, typename T_u_vec_rank,
    typename T_o_vec, typename T_o_vec_rank, typename T_o_vec_slc,
    typename T_m_vec, typename T_m_vec_rank, typename T_m_vec_slc
    >
const uint32_t Zombit<T_u_vec,T_u_vec_rank,
      T_o_vec,T_o_vec_rank,T_o_vec_slc,
      T_m_vec,T_m_vec_rank,T_m_vec_slc>::access(const uint64_t i) {
    uint64_t j = i / block_size;
    if (u_vector[j] == 1) {
        if (o_vector[j] == 1) return 1;
        else return 0;
    }
    uint64_t q = u_rank(j);
    uint64_t delta_i = i % block_size;
    uint64_t beg_q = q * block_size;
    return rec_level == 0 ? m_vector[beg_q + delta_i] : zombit_rec->access(beg_q + delta_i);
}


/*
 * Idea to scan next 1bit from x.
 * Why this would be good solution is that if bv is random or have more 1s than 0s,
 * the next 1bit is usually in same word as x or its in the next word.
 */
inline const uint64_t succ_scan(const sdsl::bit_vector &bv, const uint64_t x) {
    int next_one_bit = 0;
    uint64_t x_w_idx = x/64;
    uint64_t curr_w_bv = bv.data()[x_w_idx] >> ( x & (64-1) );
    if (curr_w_bv > 0) return x + __builtin_ctzll(curr_w_bv); // next 1bit is in the same word block
    // next possible 1bit is in some of the next word blocks
    uint64_t bv_words = (uint64_t) ceil((bv.size()/64.0));//; + 1;
    for (uint64_t i = 1; i <= ( bv_words - x_w_idx - 1 ); i++) {
        if (bv.data()[i + (x/64)] > 0) {
            return  __builtin_ctzll(
                        bv.data()[ i + x_w_idx ]
                    ) + ((i-1)*64) + x + 1 + ( 64 - (x & 63) - 1);
        }
    }
    return 0; // there is not equeal or bigger val in bv

}

// nextGEQ qeury
// block size is expected to be power of 2
template <
    typename T_u_vec, typename T_u_vec_rank,
    typename T_o_vec, typename T_o_vec_rank, typename T_o_vec_slc,
    typename T_m_vec, typename T_m_vec_rank, typename T_m_vec_slc
    >
const uint64_t Zombit<T_u_vec,T_u_vec_rank,
      T_o_vec,T_o_vec_rank,T_o_vec_slc,
      T_m_vec,T_m_vec_rank,T_m_vec_slc>::nextGEQ(const uint64_t x) const {
    uint64_t j = x / block_size;
    uint64_t q = u_rank(j);
    // x is in uniform block of 1s
    if (u_vector[j] == 1) {
        if (o_vector[j] == 1) {
            return x;
        }
    } else {
        uint64_t beg_q = q * block_size;
        //uint64_t delta_x = x % block_size;
        uint64_t delta_x = x & (block_size - 1);
        if (rec_level > 0) {
            uint64_t s = zombit_rec->nextGEQ(beg_q + delta_x);
            if (s <= beg_q + block_size - 1) {
                //return (j * block_size) + (s % block_size);
                return (j * block_size) + (s & (block_size-1));
            }
        } else { // normal nextGEQ query on m
            uint64_t next_one_in_m_block = m_rank(beg_q + delta_x) + 1;
            if (next_one_in_m_block <= m_ones) {
                uint64_t s = m_select( next_one_in_m_block );
                if (s <= beg_q + block_size - 1) {
                    //return (j * block_size) + (s % block_size);
                    return (j * block_size) + (s & (block_size-1));
                }
            }
        }
    }
    // x is in fixed block or block of 0s

    // jump func in article
    // next 1 is in next block
    //
    //////////// NORMAL SUCC with rank and select
    uint64_t next_block_with_ones = o_rank(j+1) + 1;
    if (next_block_with_ones > o_ones) {
        return 0; // no bigger or equal value than x
    }
    uint64_t j_p = o_select( next_block_with_ones );

    // next block is uniform full of 1s, return first item of block
    if (u_vector[j_p] == 1) {
        return j_p * block_size;
    }
    uint64_t beg_q1 = u_vector[j] ? q * block_size : (q+1) * block_size;
    uint64_t s_p;
    if (rec_level > 0 ) s_p = zombit_rec->nextGEQ(beg_q1);
    else s_p = m_select( m_rank( beg_q1 ) + 1 );
    uint64_t beg_j_p = j_p * block_size;
    //uint64_t delta_s_p = s_p % block_size;
    uint64_t delta_s_p = s_p & (block_size-1);
    return beg_j_p + delta_s_p;
}




template <
    typename T_u_vec, typename T_u_vec_rank,
    typename T_o_vec, typename T_o_vec_rank, typename T_o_vec_slc,
    typename T_m_vec, typename T_m_vec_rank, typename T_m_vec_slc
    >
const uint64_t Zombit<T_u_vec,T_u_vec_rank,
      T_o_vec,T_o_vec_rank,T_o_vec_slc,
      T_m_vec,T_m_vec_rank,T_m_vec_slc>::size() {
    return orig_bv_size;
}

template <
    typename T_u_vec, typename T_u_vec_rank,
    typename T_o_vec, typename T_o_vec_rank, typename T_o_vec_slc,
    typename T_m_vec, typename T_m_vec_rank, typename T_m_vec_slc
    >
const uint64_t Zombit<T_u_vec,T_u_vec_rank,
      T_o_vec,T_o_vec_rank,T_o_vec_slc,
      T_m_vec,T_m_vec_rank,T_m_vec_slc>::u_vector_size_in_bits() {

    uint64_t this_level_size = (sdsl::size_in_bytes(u_vector) * 8) + \
        (sdsl::size_in_bytes(u_rank) * 8);
    if (rec_level > 0) return this_level_size + zombit_rec->u_vector_size_in_bits();
    return this_level_size;
}

// o_vector size
template <
    typename T_u_vec, typename T_u_vec_rank,
    typename T_o_vec, typename T_o_vec_rank, typename T_o_vec_slc,
    typename T_m_vec, typename T_m_vec_rank, typename T_m_vec_slc
    >
const uint64_t Zombit<T_u_vec,T_u_vec_rank,
          T_o_vec,T_o_vec_rank,T_o_vec_slc,
          T_m_vec,T_m_vec_rank,T_m_vec_slc>::o_vector_size_in_bits() {
    uint64_t this_level_size;
    this_level_size = (sdsl::size_in_bytes(o_vector) * 8) + \
        (sdsl::size_in_bytes(o_rank) * 8) + \
        (sdsl::size_in_bytes(o_select) * 8);


    if (rec_level > 0) return this_level_size + zombit_rec->o_vector_size_in_bits();
    return this_level_size;
}

// m_vector size
template <
    typename T_u_vec, typename T_u_vec_rank,
    typename T_o_vec, typename T_o_vec_rank, typename T_o_vec_slc,
    typename T_m_vec, typename T_m_vec_rank, typename T_m_vec_slc
    >
const uint64_t Zombit<T_u_vec,T_u_vec_rank,
      T_o_vec,T_o_vec_rank,T_o_vec_slc,
      T_m_vec,T_m_vec_rank,T_m_vec_slc>::m_vector_size_in_bits() {
    if (rec_level > 0) return zombit_rec->m_vector_size_in_bits();
    return (sdsl::size_in_bytes(m_vector) * 8) + \
        (sdsl::size_in_bytes(m_rank) * 8) + \
        (sdsl::size_in_bytes(m_select) * 8);
}

// total size
template <
    typename T_u_vec, typename T_u_vec_rank,
    typename T_o_vec, typename T_o_vec_rank, typename T_o_vec_slc,
    typename T_m_vec, typename T_m_vec_rank, typename T_m_vec_slc
    >
const uint64_t Zombit<T_u_vec,T_u_vec_rank,
      T_o_vec,T_o_vec_rank,T_o_vec_slc,
      T_m_vec,T_m_vec_rank,T_m_vec_slc>::size_in_bits() {
    return u_vector_size_in_bits() + o_vector_size_in_bits() + m_vector_size_in_bits();
}





















//////////////////////////////////////////////
//
// Zombit vector with local super block scan nextGEQ O-vector
//
//////////////////////////////////////////////

template <
    typename T_u_vec, typename T_u_vec_rank,
    typename T_o_vec_slc,
    typename T_m_vec, typename T_m_vec_rank, typename T_m_vec_slc
    >
class Zombit<T_u_vec, T_u_vec_rank,
                sdsl::bit_vector, sdsl::rank_support_v5O2<1>, T_o_vec_slc,
                T_m_vec, T_m_vec_rank, T_m_vec_slc> {
    public:

        Zombit();
        ~Zombit();
        void build_zombit(sdsl::bit_vector &X, const uint32_t b, const uint32_t recursio_level, std::string rec_model);
        void build_zombit(sdsl::bit_vector &X, const uint32_t recursio_level, std::string rec_model);
        const uint32_t access(const uint64_t i); // return value of i-th bit
        const uint64_t nextGEQ(const uint64_t x) const; // return next greater or equal than x
        const uint64_t size(); // return lenght of the original bitvector
        const uint64_t size_in_bits();
        const uint64_t u_vector_size_in_bits();
        const uint64_t o_vector_size_in_bits();
        const uint64_t m_vector_size_in_bits();

        uint32_t block_size;
        uint64_t orig_bv_size = 0;
        uint64_t runs_n = 0;
        uint64_t m_blocks = 0;
        uint64_t block_n = 0;
        uint32_t rec_level = 0;
        std::string b_model = "";

        T_u_vec u_vector;
        sdsl::bit_vector o_vector;
        T_m_vec m_vector;
        Zombit<T_u_vec, T_u_vec_rank,\
            sdsl::bit_vector, sdsl::rank_support_v5O2<1>, T_o_vec_slc,\
            T_m_vec, T_m_vec_rank, T_m_vec_slc>* zombit_rec;

        T_u_vec_rank u_rank;
        sdsl::rank_support_v5O2<1> o_rank;
        T_o_vec_slc o_select;
        T_m_vec_rank m_rank;
        T_m_vec_slc m_select;

        size_t o_ones;
        size_t m_ones;

        const uint32_t block_model(const sdsl::bit_vector &M);
};


template <
    typename T_u_vec, typename T_u_vec_rank,
    typename T_o_vec_slc,
    typename T_m_vec, typename T_m_vec_rank, typename T_m_vec_slc
    >
const uint32_t Zombit<T_u_vec,T_u_vec_rank,
        sdsl::bit_vector, sdsl::rank_support_v5O2<1>, T_o_vec_slc,
        T_m_vec,T_m_vec_rank,T_m_vec_slc>::block_model(const sdsl::bit_vector &M) {
    if (b_model == "div2") return block_size/2;
    if (b_model == "div4") return block_size/4;
    if (b_model == "div8") return block_size/8;
    return 1;
}

template <
    typename T_u_vec, typename T_u_vec_rank,
    typename T_o_vec_slc,
    typename T_m_vec, typename T_m_vec_rank, typename T_m_vec_slc
    >
Zombit<T_u_vec,T_u_vec_rank,
        sdsl::bit_vector, sdsl::rank_support_v5O2<1>, T_o_vec_slc,
        T_m_vec,T_m_vec_rank,T_m_vec_slc>::Zombit() {};

template <
    typename T_u_vec, typename T_u_vec_rank,
    typename T_o_vec_slc,
    typename T_m_vec, typename T_m_vec_rank, typename T_m_vec_slc
    >
Zombit<T_u_vec,T_u_vec_rank,
        sdsl::bit_vector, sdsl::rank_support_v5O2<1>, T_o_vec_slc,
        T_m_vec,T_m_vec_rank,T_m_vec_slc>::~Zombit() {
    if (rec_level > 0) delete zombit_rec;
}

/*
 * Note if |X| mod block_size != 0, bitvector X is modified to zero's are added end of X, such that
 * |X| mod block_size = 0 is valid.
 */
template <
    typename T_u_vec, typename T_u_vec_rank,
    typename T_o_vec_slc,
    typename T_m_vec, typename T_m_vec_rank, typename T_m_vec_slc
    >
void Zombit<T_u_vec,T_u_vec_rank,
        sdsl::bit_vector, sdsl::rank_support_v5O2<1>, T_o_vec_slc,
        T_m_vec,T_m_vec_rank,T_m_vec_slc>::build_zombit(
        sdsl::bit_vector &X, const uint32_t b, const uint32_t recursio_level, std::string rec_model ) {
    block_size = b;
    orig_bv_size = X.size();
    rec_level = recursio_level;
    b_model = rec_model;
    if (runs_n == 0) {
        uint32_t prev = 0;
        for (uint64_t i = 0; i < X.size(); i++) {
            if (X[i] && !prev) {
                runs_n++;
                prev = 1;
            } else if (!X[i] && prev) {
                prev = 0;
            }
        }
    }
    // adding 0s end of X, so that all blocks have same length.
    if (X.size() % block_size != 0) {
        X.resize(X.size() + (block_size - (X.size() % block_size)));
    }
    block_n = X.size() / block_size;
    sdsl::bit_vector U_bv = sdsl::bit_vector(block_n);
    sdsl::bit_vector O_bv = sdsl::bit_vector(block_n);
    sdsl::bit_vector M_bv;
    char label = 'O';
    m_blocks = 0;
    std::vector<uint64_t> m_block_idx;

    // building U and O-vector
    for (uint64_t i = 0; i < block_n; i++) {
        //block first val
        if (X[i*block_size] == 0) label = 'Z';
        else label = 'O';
        for (uint64_t j = 1; j < block_size; j++) {

            if (X[i * block_size + j] == 0 && label == 'Z') continue;
            if (X[i * block_size + j] == 1 && label == 'O') continue;
            if (X[i * block_size + j] == 0 && label == 'O') {
                label = 'M';
                continue;
            }
            if (X[ (i * block_size) + j] == 1 && label == 'Z') {
                label = 'M';
                continue;
            }
        }
        // setting value for U-vector at index b_idx
        if (label == 'Z' || label == 'O') U_bv[i] = 1; // all values 0s or 1s
        else U_bv[i] = 0; // mix of 0s and 1s
        // setting value for O-vector at indexx b_idx
        if (label == 'Z') O_bv[i] = 0; // all 0s
        else O_bv[i] = 1; // block contain atleast one 1 bit

        if (label == 'M') {
            m_blocks++;
            m_block_idx.push_back(i);
        }
    }
    // building M-vector
    M_bv = sdsl::bit_vector(block_size * m_blocks);
    uint64_t m_idx = 0;
    for (uint64_t x : m_block_idx) {
        for (uint64_t k = 0; k < block_size; k++) {
            M_bv[m_idx] = X[(x * block_size) + k];
            m_idx++;
        }
    }

    m_block_idx.clear();
    m_block_idx.shrink_to_fit();

    // constructing U and O-vector and it's rank/select support
    u_vector = T_u_vec(U_bv);
    o_vector = sdsl::bit_vector(O_bv);
    u_rank = T_u_vec_rank(&u_vector);
    o_rank = sdsl::rank_support_v5O2<1>(&o_vector);
    o_select = T_o_vec_slc(&o_vector);

    if (rec_level > 0) {
        uint32_t m_b = block_size / 4;
        if (m_b <= 1) {
            rec_level = 0;
            rec_level = false;
        } else {
            zombit_rec = new Zombit<T_u_vec, T_u_vec_rank,\
                         sdsl::bit_vector, sdsl::rank_support_v5O2<1>, T_o_vec_slc,\
                         T_m_vec, T_m_vec_rank, T_m_vec_slc>();
            zombit_rec->build_zombit(M_bv, m_b, rec_level-1, rec_model);
        }
    }
    if (rec_level == 0) {
        m_vector = T_m_vec(M_bv);
        m_rank = T_m_vec_rank(&m_vector);
        m_select = T_m_vec_slc(&m_vector);
        m_ones = T_m_vec_rank(&m_vector)(M_bv.size());
    }

    o_ones = sdsl::rank_support_v5O2<1>(&o_vector)(O_bv.size());
}

// build zombit-vector with optimal block size sqrt(n/runs) according article
template <
    typename T_u_vec, typename T_u_vec_rank,
    typename T_o_vec_slc,
    typename T_m_vec, typename T_m_vec_rank, typename T_m_vec_slc
    >
void Zombit<T_u_vec,T_u_vec_rank,
        sdsl::bit_vector, sdsl::rank_support_v5O2<1>, T_o_vec_slc,
        T_m_vec,T_m_vec_rank,T_m_vec_slc>::build_zombit(
        sdsl::bit_vector &X, const uint32_t recursio_level, std::string rec_model) {
    uint32_t prev = 0;
    rec_level = recursio_level;
    runs_n = 0;
    for (uint64_t i = 0; i < X.size(); i++) {
        if (X[i] && !prev) {
            runs_n++;
            prev = 1;
        } else if (!X[i] && prev) {
            prev = 0;
        }
    }
    if (runs_n == 0) runs_n = 1; // all 0s or 1s
    uint32_t b = (uint32_t) sqrt(X.size()/runs_n);
    //if (b == 1) b++;

    build_zombit(X, b, recursio_level, rec_model);
}


template <
    typename T_u_vec, typename T_u_vec_rank,
    typename T_o_vec_slc,
    typename T_m_vec, typename T_m_vec_rank, typename T_m_vec_slc
    >
const uint32_t Zombit<T_u_vec,T_u_vec_rank,
        sdsl::bit_vector, sdsl::rank_support_v5O2<1>, T_o_vec_slc,
        T_m_vec,T_m_vec_rank,T_m_vec_slc>::access(const uint64_t i) {
    uint64_t j = i / block_size;
    if (u_vector[j] == 1) {
        if (o_vector[j] == 1) return 1;
        else return 0;
    }
    uint64_t q = u_rank(j);
    uint64_t delta_i = i % block_size;
    uint64_t beg_q = q * block_size;
    return rec_level == 0 ? m_vector[beg_q + delta_i] : zombit_rec->access(beg_q + delta_i);
}


// local nextGEQ operation for O-vector is computed with super block scanning.
template <
    typename T_u_vec, typename T_u_vec_rank,
    typename T_o_vec_slc,
    typename T_m_vec, typename T_m_vec_rank, typename T_m_vec_slc
    >
const uint64_t Zombit< T_u_vec, T_u_vec_rank,
        sdsl::bit_vector, sdsl::rank_support_v5O2<1>, T_o_vec_slc,
        T_m_vec, T_m_vec_rank, T_m_vec_slc>::nextGEQ(const uint64_t x) const {
    uint64_t j = x / block_size;
    uint64_t q = u_rank(j);
    // x is in uniform block of 1s
    if (u_vector[j] == 1) {
        if (o_vector[j] == 1) {
            return x;
        }
    } else {
        uint64_t beg_q = q * block_size;
        //uint64_t delta_x = x % block_size;
        uint64_t delta_x = x & (block_size - 1);
        if (rec_level > 0) {
            uint64_t s = zombit_rec->nextGEQ(beg_q + delta_x);
            if (s <= beg_q + block_size - 1) {
                //return (j * block_size) + (s % block_size);
                return (j * block_size) + (s & (block_size-1));
            }
        } else { // normal nextGEQ query on m
            uint64_t next_one_in_m_block = m_rank(beg_q + delta_x) + 1;
            if (next_one_in_m_block <= m_ones) {
                uint64_t s = m_select( next_one_in_m_block );
                if (s <= beg_q + block_size - 1) {
                    //return (j * block_size) + (s % block_size);
                    return (j * block_size) + (s & (block_size-1));
                }
            }
        }
    }
    // x is in fixed block or block of 0s

    // jump func in article
    // next 1 is in next block
    uint64_t j_p = o_rank.nextGEQ(o_vector, j+1);

    // next block is uniform full of 1s, return first item of block
    if (u_vector[j_p] == 1) {
        return j_p * block_size;
    }
    uint64_t beg_q1 = u_vector[j] ? q * block_size : (q+1) * block_size;
    uint64_t s_p;
    if (rec_level > 0 ) s_p = zombit_rec->nextGEQ(beg_q1);
    else s_p = m_select( m_rank( beg_q1 ) + 1 );
    uint64_t beg_j_p = j_p * block_size;
    //uint64_t delta_s_p = s_p % block_size;
    uint64_t delta_s_p = s_p & (block_size-1);
    return beg_j_p + delta_s_p;
}

// length of the original bitvector
template <
    typename T_u_vec, typename T_u_vec_rank,
    typename T_o_vec_slc,
    typename T_m_vec, typename T_m_vec_rank, typename T_m_vec_slc
    >
const uint64_t Zombit<T_u_vec,T_u_vec_rank,
        sdsl::bit_vector, sdsl::rank_support_v5O2<1>, T_o_vec_slc,
      T_m_vec,T_m_vec_rank,T_m_vec_slc>::size() {
    return orig_bv_size;
}

// u_vector size
template <
    typename T_u_vec, typename T_u_vec_rank,
    typename T_o_vec_slc,
    typename T_m_vec, typename T_m_vec_rank, typename T_m_vec_slc
    >
const uint64_t Zombit<T_u_vec,T_u_vec_rank,
        sdsl::bit_vector, sdsl::rank_support_v5O2<1>, T_o_vec_slc,
        T_m_vec,T_m_vec_rank,T_m_vec_slc>::u_vector_size_in_bits() {

    uint64_t this_level_size = (sdsl::size_in_bytes(u_vector) * 8) + \
        (sdsl::size_in_bytes(u_rank) * 8);
    if (rec_level > 0) return this_level_size + zombit_rec->u_vector_size_in_bits();
    return this_level_size;
}

// o_vector size
template <
    typename T_u_vec, typename T_u_vec_rank,
    typename T_o_vec_slc,
    typename T_m_vec, typename T_m_vec_rank, typename T_m_vec_slc
    >
const uint64_t Zombit<T_u_vec,T_u_vec_rank,
        sdsl::bit_vector, sdsl::rank_support_v5O2<1>, T_o_vec_slc,
        T_m_vec,T_m_vec_rank,T_m_vec_slc>::o_vector_size_in_bits() {
    uint64_t this_level_size;
    this_level_size = (sdsl::size_in_bytes(o_vector)*8) + o_rank.nextGEQ_rank_size_in_bits();

    if (rec_level > 0) return this_level_size + zombit_rec->o_vector_size_in_bits();
    return this_level_size;
}

// m_vector size
template <
    typename T_u_vec, typename T_u_vec_rank,
    typename T_o_vec_slc,
    typename T_m_vec, typename T_m_vec_rank, typename T_m_vec_slc
    >
const uint64_t Zombit<T_u_vec,T_u_vec_rank,
        sdsl::bit_vector, sdsl::rank_support_v5O2<1>, T_o_vec_slc,
        T_m_vec,T_m_vec_rank,T_m_vec_slc>::m_vector_size_in_bits() {
    if (rec_level > 0) return zombit_rec->m_vector_size_in_bits();
    return (sdsl::size_in_bytes(m_vector) * 8) + \
        (sdsl::size_in_bytes(m_rank) * 8) + \
        (sdsl::size_in_bytes(m_select) * 8);
}

// total size
template <
    typename T_u_vec, typename T_u_vec_rank,
    typename T_o_vec_slc,
    typename T_m_vec, typename T_m_vec_rank, typename T_m_vec_slc
    >
const uint64_t Zombit<T_u_vec,T_u_vec_rank,
        sdsl::bit_vector, sdsl::rank_support_v5O2<1>, T_o_vec_slc,
        T_m_vec,T_m_vec_rank,T_m_vec_slc>::size_in_bits() {
    return u_vector_size_in_bits() + o_vector_size_in_bits() + m_vector_size_in_bits();
}















//////////////////////////////////////////////
//
// Zombit vector with local super block scan nextGEQ for O,M-vectors
//
//////////////////////////////////////////////

template <
    typename T_u_vec, typename T_u_vec_rank,
    typename T_o_vec_slc,
    typename T_m_vec_slc
    >
class Zombit<T_u_vec, T_u_vec_rank,
            sdsl::bit_vector, sdsl::rank_support_v5O2<1>, T_o_vec_slc,
            sdsl::bit_vector, sdsl::rank_support_v5O2<1>, T_m_vec_slc> {
    public:

        Zombit();
        ~Zombit();
        void build_zombit(sdsl::bit_vector &X, const uint32_t b, const uint32_t recursio_level, std::string rec_model);
        void build_zombit(sdsl::bit_vector &X, const uint32_t recursio_level, std::string rec_model);
        const uint32_t access(const uint64_t i); // return value of i-th bit
        const uint64_t nextGEQ(const uint64_t x) const; // return next greater or equal than x
        const uint64_t size(); // return lenght of the original bitvector
        const uint64_t size_in_bits();
        const uint64_t u_vector_size_in_bits();
        const uint64_t o_vector_size_in_bits();
        const uint64_t m_vector_size_in_bits();

        uint32_t block_size;
        uint64_t orig_bv_size = 0;
        uint64_t runs_n = 0;
        uint64_t m_blocks = 0;
        uint64_t block_n = 0;
        uint32_t rec_level = 0;
        std::string b_model = "";

        T_u_vec u_vector;
        sdsl::bit_vector o_vector;
        sdsl::bit_vector m_vector;
        Zombit<T_u_vec, T_u_vec_rank,\
            sdsl::bit_vector, sdsl::rank_support_v5O2<1>, T_o_vec_slc,\
            sdsl::bit_vector, sdsl::rank_support_v5O2<1>, T_m_vec_slc>* zombit_rec;

        T_u_vec_rank u_rank;
        sdsl::rank_support_v5O2<1> o_rank;
        T_o_vec_slc o_select;
        sdsl::rank_support_v5O2<1> m_rank;
        T_m_vec_slc m_select;

        size_t o_ones;
        size_t m_ones;

        const uint32_t block_model(const sdsl::bit_vector &M);
};


template <
    typename T_u_vec, typename T_u_vec_rank,
    typename T_o_vec_slc,
    typename T_m_vec_slc
    >
const uint32_t Zombit<T_u_vec,T_u_vec_rank,
      sdsl::bit_vector, sdsl::rank_support_v5O2<1>,T_o_vec_slc,
      sdsl::bit_vector, sdsl::rank_support_v5O2<1>,T_m_vec_slc
      >::block_model(const sdsl::bit_vector &M) {
    if (b_model == "div2") return block_size/2;
    if (b_model == "div4") return block_size/4;
    if (b_model == "div8") return block_size/8;
    return 1;
}

template <
    typename T_u_vec, typename T_u_vec_rank,
    typename T_o_vec_slc,
    typename T_m_vec_slc
    >
Zombit<T_u_vec,T_u_vec_rank,
      sdsl::bit_vector, sdsl::rank_support_v5O2<1>,T_o_vec_slc,
      sdsl::bit_vector, sdsl::rank_support_v5O2<1>,T_m_vec_slc
    >::Zombit() {};

template <
    typename T_u_vec, typename T_u_vec_rank,
    typename T_o_vec_slc,
    typename T_m_vec_slc
    >
Zombit<T_u_vec,T_u_vec_rank,
      sdsl::bit_vector, sdsl::rank_support_v5O2<1>,T_o_vec_slc,
      sdsl::bit_vector, sdsl::rank_support_v5O2<1>,T_m_vec_slc
    >::~Zombit() {
    if (rec_level > 0) delete zombit_rec;
}

/*
 * Note if |X| mod block_size != 0, bitvector X is modified to zero's are added end of X, such that
 * |X| mod block_size = 0 is valid.
 */
template <
    typename T_u_vec, typename T_u_vec_rank,
    typename T_o_vec_slc,
    typename T_m_vec_slc
    >
void Zombit<T_u_vec,T_u_vec_rank,
      sdsl::bit_vector, sdsl::rank_support_v5O2<1>,T_o_vec_slc,
      sdsl::bit_vector, sdsl::rank_support_v5O2<1>,T_m_vec_slc
     >::build_zombit(
        sdsl::bit_vector &X, const uint32_t b, const uint32_t recursio_level, std::string rec_model ) {
    block_size = b;
    orig_bv_size = X.size();
    rec_level = recursio_level;
    b_model = rec_model;
    if (runs_n == 0) {
        uint32_t prev = 0;
        for (uint64_t i = 0; i < X.size(); i++) {
            if (X[i] && !prev) {
                runs_n++;
                prev = 1;
            } else if (!X[i] && prev) {
                prev = 0;
            }
        }
    }
    // adding 0s end of X, so that all blocks have same length.
    if (X.size() % block_size != 0) {
        X.resize(X.size() + (block_size - (X.size() % block_size)));
    }
    block_n = X.size() / block_size;
    sdsl::bit_vector U_bv = sdsl::bit_vector(block_n);
    sdsl::bit_vector O_bv = sdsl::bit_vector(block_n);
    sdsl::bit_vector M_bv;
    char label = 'O';
    m_blocks = 0;
    std::vector<uint64_t> m_block_idx;

    // building U and O-vector
    for (uint64_t i = 0; i < block_n; i++) {
        //block first val
        if (X[i*block_size] == 0) label = 'Z';
        else label = 'O';
        for (uint64_t j = 1; j < block_size; j++) {

            if (X[i * block_size + j] == 0 && label == 'Z') continue;
            if (X[i * block_size + j] == 1 && label == 'O') continue;
            if (X[i * block_size + j] == 0 && label == 'O') {
                label = 'M';
                continue;
            }
            if (X[ (i * block_size) + j] == 1 && label == 'Z') {
                label = 'M';
                continue;
            }
        }
        // setting value for U-vector at index b_idx
        if (label == 'Z' || label == 'O') U_bv[i] = 1; // all values 0s or 1s
        else U_bv[i] = 0; // mix of 0s and 1s
        // setting value for O-vector at indexx b_idx
        if (label == 'Z') O_bv[i] = 0; // all 0s
        else O_bv[i] = 1; // block contain atleast one 1 bit

        if (label == 'M') {
            m_blocks++;
            m_block_idx.push_back(i);
        }
    }
    // building M-vector
    M_bv = sdsl::bit_vector(block_size * m_blocks);
    uint64_t m_idx = 0;
    for (uint64_t x : m_block_idx) {
        for (uint64_t k = 0; k < block_size; k++) {
            M_bv[m_idx] = X[(x * block_size) + k];
            m_idx++;
        }
    }

    m_block_idx.clear();
    m_block_idx.shrink_to_fit();

    // constructing U and O-vector and it's rank/select support
    u_vector = T_u_vec(U_bv);
    o_vector = sdsl::bit_vector(O_bv);
    u_rank = T_u_vec_rank(&u_vector);
    o_rank = sdsl::rank_support_v5O2<1>(&o_vector);
    o_select = T_o_vec_slc(&o_vector);

    if (rec_level > 0) {
        uint32_t m_b = block_model(M_bv);
        if (m_b <= 1) {
            rec_level = 0;
            rec_level = false;
        } else {
            zombit_rec = new Zombit<T_u_vec, T_u_vec_rank,\
                         sdsl::bit_vector, sdsl::rank_support_v5O2<1>, T_o_vec_slc,\
                         sdsl::bit_vector, sdsl::rank_support_v5O2<1>, T_m_vec_slc>();
            zombit_rec->build_zombit(M_bv, m_b, rec_level-1, rec_model);
        }
    }
    if (rec_level == 0) {
        m_vector = sdsl::bit_vector(M_bv);
        m_rank = sdsl::rank_support_v5O2<1>(&m_vector);
        m_select = T_m_vec_slc(&m_vector);
        m_ones = sdsl::rank_support_v5O2<1>(&m_vector)(M_bv.size());
    }

    o_ones = sdsl::rank_support_v5O2<1>(&o_vector)(O_bv.size());
}

// build zombit-vector with optimal block size sqrt(n/runs) according article
template <
    typename T_u_vec, typename T_u_vec_rank,
    typename T_o_vec_slc,
    typename T_m_vec_slc
    >
void Zombit<T_u_vec,T_u_vec_rank,
      sdsl::bit_vector, sdsl::rank_support_v5O2<1>,T_o_vec_slc,
      sdsl::bit_vector, sdsl::rank_support_v5O2<1>,T_m_vec_slc
     >::build_zombit(
        sdsl::bit_vector &X, const uint32_t recursio_level, std::string rec_model) {
    uint32_t prev = 0;
    rec_level = recursio_level;
    runs_n = 0;
    for (uint64_t i = 0; i < X.size(); i++) {
        if (X[i] && !prev) {
            runs_n++;
            prev = 1;
        } else if (!X[i] && prev) {
            prev = 0;
        }
    }
    if (runs_n == 0) runs_n = 1; // all 0s or 1s
    uint32_t b = (uint32_t) sqrt(X.size()/runs_n);
    //if (b == 1) b++;

    build_zombit(X, b, recursio_level, rec_model);
}


template <
    typename T_u_vec, typename T_u_vec_rank,
    typename T_o_vec_slc,
    typename T_m_vec_slc
    >
const uint32_t Zombit<T_u_vec,T_u_vec_rank,
      sdsl::bit_vector, sdsl::rank_support_v5O2<1>,T_o_vec_slc,
      sdsl::bit_vector, sdsl::rank_support_v5O2<1>,T_m_vec_slc
      >::access(const uint64_t i) {
    uint64_t j = i / block_size;
    if (u_vector[j] == 1) {
        if (o_vector[j] == 1) return 1;
        else return 0;
    }
    uint64_t q = u_rank(j);
    uint64_t delta_i = i % block_size;
    uint64_t beg_q = q * block_size;
    return rec_level == 0 ? m_vector[beg_q + delta_i] : zombit_rec->access(beg_q + delta_i);
}

// local nextGEQ operation for O and M-vector is computed with super block scanning.
template <
    typename T_u_vec, typename T_u_vec_rank,
    typename T_o_vec_slc,
    typename T_m_vec_slc
    >
const uint64_t Zombit<T_u_vec,T_u_vec_rank,
      sdsl::bit_vector, sdsl::rank_support_v5O2<1>,T_o_vec_slc,
      sdsl::bit_vector, sdsl::rank_support_v5O2<1>,T_m_vec_slc
      >::nextGEQ(const uint64_t x) const {
    uint64_t j = x / block_size;
    uint64_t q = u_rank(j);
    // x is in uniform block of 1s
    if (u_vector[j] == 1) {
        if (o_vector[j] == 1) {
            return x;
        }
    } else {
        uint64_t beg_q = q * block_size;
        //uint64_t delta_x = x % block_size;
        uint64_t delta_x = x & (block_size - 1);
        if (rec_level > 0) {
            uint64_t s = zombit_rec->nextGEQ(beg_q + delta_x);
            if (s <= beg_q + block_size - 1) {
                //return (j * block_size) + (s % block_size);
                return (j * block_size) + (s & (block_size-1));
            }
        } else { // normal nextGEQ query on m
            uint64_t s = m_rank.nextGEQ(m_vector, beg_q + delta_x);
            if (s <= beg_q + block_size - 1) {
                //return (j * block_size) + (s % block_size);
                return (j * block_size) + (s & (block_size-1));
            }
        }
    }
    // x is in fixed block or block of 0s

    // jump func in article
    // next 1 is in next block
    //
    uint64_t j_p = o_rank.nextGEQ(o_vector, j+1);

    // next block is uniform full of 1s, return first item of block
    if (u_vector[j_p] == 1) {
        return j_p * block_size;
    }
    uint64_t beg_q1 = u_vector[j] ? q * block_size : (q+1) * block_size;
    uint64_t s_p;
    if (rec_level > 0 ) s_p = zombit_rec->nextGEQ(beg_q1);
    else s_p = m_rank.nextGEQ(m_vector, beg_q1);
    uint64_t beg_j_p = j_p * block_size;
    //uint64_t delta_s_p = s_p % block_size;
    uint64_t delta_s_p = s_p & (block_size-1);
    return beg_j_p + delta_s_p;
}

// length of the original bitvector
template <
    typename T_u_vec, typename T_u_vec_rank,
    typename T_o_vec_slc,
    typename T_m_vec_slc
    >
const uint64_t Zombit<T_u_vec,T_u_vec_rank,
      sdsl::bit_vector, sdsl::rank_support_v5O2<1>,T_o_vec_slc,
      sdsl::bit_vector, sdsl::rank_support_v5O2<1>,T_m_vec_slc
      >::size() {
    return orig_bv_size;
}

// u_vector size
template <
    typename T_u_vec, typename T_u_vec_rank,
    typename T_o_vec_slc,
    typename T_m_vec_slc
    >
const uint64_t Zombit<T_u_vec,T_u_vec_rank,
      sdsl::bit_vector, sdsl::rank_support_v5O2<1>,T_o_vec_slc,
      sdsl::bit_vector, sdsl::rank_support_v5O2<1>,T_m_vec_slc
      >::u_vector_size_in_bits() {

    uint64_t this_level_size = (sdsl::size_in_bytes(u_vector) * 8) + \
        (sdsl::size_in_bytes(u_rank) * 8);
    if (rec_level > 0) return this_level_size + zombit_rec->u_vector_size_in_bits();
    return this_level_size;
}

// o_vector size
template <
    typename T_u_vec, typename T_u_vec_rank,
    typename T_o_vec_slc,
    typename T_m_vec_slc
    >
const uint64_t Zombit<T_u_vec,T_u_vec_rank,
      sdsl::bit_vector, sdsl::rank_support_v5O2<1>,T_o_vec_slc,
      sdsl::bit_vector, sdsl::rank_support_v5O2<1>,T_m_vec_slc>::o_vector_size_in_bits() {
    uint64_t this_level_size;
    this_level_size = (sdsl::size_in_bytes(o_vector)*8) + o_rank.nextGEQ_rank_size_in_bits();

    if (rec_level > 0) return this_level_size + zombit_rec->o_vector_size_in_bits();
    return this_level_size;
}

// m_vector size
template <
    typename T_u_vec, typename T_u_vec_rank,
    typename T_o_vec_slc,
    typename T_m_vec_slc
    >
const uint64_t Zombit<T_u_vec,T_u_vec_rank,
      sdsl::bit_vector, sdsl::rank_support_v5O2<1>,T_o_vec_slc,
      sdsl::bit_vector, sdsl::rank_support_v5O2<1>,T_m_vec_slc>::m_vector_size_in_bits() {
    if (rec_level > 0) return zombit_rec->m_vector_size_in_bits();
    return (sdsl::size_in_bytes(m_vector)*8) + (m_rank.nextGEQ_rank_size_in_bits());
}

// total size
template <
    typename T_u_vec, typename T_u_vec_rank,
    typename T_o_vec_slc,
    typename T_m_vec_slc
    >
const uint64_t Zombit<T_u_vec,T_u_vec_rank,
      sdsl::bit_vector, sdsl::rank_support_v5O2<1>,T_o_vec_slc,
      sdsl::bit_vector, sdsl::rank_support_v5O2<1>,T_m_vec_slc>::size_in_bits() {
    return u_vector_size_in_bits() + o_vector_size_in_bits() + m_vector_size_in_bits();
}



//// type definitions for Zombit

// bv_bv_bv_O2
typedef Zombit<sdsl::bit_vector, sdsl::rank_support_v5O2<0>,sdsl::bit_vector, sdsl::rank_support_v5O2<1>, sdsl::select_support_mcl<1>,sdsl::bit_vector, sdsl::rank_support_v5O2<1>, sdsl::select_support_mcl<1>> zombit_bv_bv_bv_O2;
// bv_bv_bv
typedef Zombit<sdsl::bit_vector, sdsl::rank_support_v5<0>,sdsl::bit_vector, sdsl::rank_support_v5<1>, sdsl::select_support_mcl<1>,sdsl::bit_vector, sdsl::rank_support_v5<1>, sdsl::select_support_mcl<1>> zombit_bv_bv_bv;
// bv_bv_rrr
typedef Zombit<sdsl::bit_vector, sdsl::rank_support_v5<0>,sdsl::bit_vector, sdsl::rank_support_v5<1>, sdsl::select_support_mcl<1>,sdsl::rrr_vector<>, sdsl::rrr_vector<>::rank_1_type, sdsl::rrr_vector<>::select_1_type> zombit_bv_bv_rrr;
// bv_bv_bvIL
typedef Zombit<sdsl::bit_vector, sdsl::rank_support_v5<0>,sdsl::bit_vector, sdsl::rank_support_v5<1>, sdsl::select_support_mcl<1>,sdsl::bit_vector_il<>, sdsl::rank_support_il<1>, sdsl::select_support_il<1>> zombit_bv_bv_bvIL;
// bv_bv_sd
typedef Zombit<sdsl::bit_vector, sdsl::rank_support_v5<0>,sdsl::bit_vector, sdsl::rank_support_v5<1>, sdsl::select_support_mcl<1>,sdsl::sd_vector<>, sdsl::sd_vector<>::rank_1_type, sdsl::sd_vector<>::select_1_type> zombit_bv_bv_sd;
// rrr_rrr_rrr
typedef Zombit<sdsl::rrr_vector<>, sdsl::rrr_vector<>::rank_0_type,sdsl::rrr_vector<>, sdsl::rrr_vector<>::rank_1_type, sdsl::rrr_vector<>::select_1_type,sdsl::rrr_vector<>, sdsl::rrr_vector<>::rank_1_type, sdsl::rrr_vector<>::select_1_type> zombit_rrr_rrr_rrr;
// hyb_rrr_rrr
typedef Zombit<sdsl::hyb_vector<>, sdsl::hyb_vector<>::rank_0_type,sdsl::rrr_vector<>, sdsl::rrr_vector<>::rank_1_type, sdsl::rrr_vector<>::select_1_type,sdsl::rrr_vector<>, sdsl::rrr_vector<>::rank_1_type, sdsl::rrr_vector<>::select_1_type> zombit_hyb_rrr_rrr;

#endif
