# Zombit-vector

Practical implementation of Zombit-vector from [A. Gómez-Brandón. “Bitvectors with Runs and the Successor/Predecessor Problem”.
In: 2020 Data Compression Conference (DCC). 2020](https://ieeexplore.ieee.org/document/9105872).


## Setup / Install

* Zombit-vector has dependencies on [SDSL](https://github.com/Hiipivahalko/sdsl-lite.git)-libary (see install instructions from link and not it is not the original repo) example to
`<repo-root>/external/cpp_libs`.
* Zombit implementation is template, so to use Zombit-vector, copy content from `./src/include` directory to your include path.

### Example usage

```c++
// main.cpp
#include <iostream>
#include <zombit/zombit_vector.hpp>
#include <sdsl/bit_vectors.hpp>

// bv_bv_bv_O2
typedef Zombit<sdsl::bit_vector, sdsl::rank_support_v5O2<0>,
                sdsl::bit_vector, sdsl::rank_support_v5O2<1>, sdsl::select_support_mcl<1>,
                sdsl::bit_vector, sdsl::rank_support_v5O2<1>, sdsl::select_support_mcl<1>> zombit_bv_bv_bv_O2;

int main() {

    sdsl::bit_vector bv = {1,1,0,1,1,1,1,0,0,0,0,1,0,1,1};
    uint32_t block_size = 32;
    zombit_bv_bv_bv_O2 zombit{};
    zombit.build_zombit(bv, block_size, 0, "");

    std::cout << "nextGEQ(9) = " << zombit.nextGEQ(9) << "\n";
    return 0;
}
```

Then to compile and run the program with `g++`:
```bash
g++ -std=c++11 -O3 -I~/include -L~/lib main.cpp -o main -lsdsl -ldivsufsort -ldivsufsort64
./main
```

### Dev config

For VIM run following command after external lib installation for coc.vim to find sources:

```bash
bear -- make
```
This command creates file `compile_commands.json`.

## TEST libraries

### test setup

```bash
mkdir -p external/cpp_libs
cd external
git clone https://github.com/Hiipivahalko/sdsl-lite.git
cd sdsl-lite
./install ../cpp_libs
```

### running tests

Run scriptfile `run_all_tests.sh` in project root directory:
```bash
./run_all_tests.sh
```

For running a specific testsets:
```bash
./run_all_tests.sh input <prefix_of_test_case>
```
