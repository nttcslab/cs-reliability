#ifndef EQOPT_COMMON_HPP
#define EQOPT_COMMON_HPP

#include <cmath>
#include <cstdint>
#include <cstddef>
#include <cassert>
#include <utility>

using addr_t = long long int;
using Bint = long long int;

// FNV hash constant
constexpr uint64_t FNV_OFFSET_BASIS_64 = 14695981039346656037ULL;
constexpr uint64_t FNV_PRIME_64 = 1099511628211ULL;

// hash function of pair of ints
class HashPI{
public:
  size_t operator()(const std::pair<int, int>& x) const {
    uint64_t h = FNV_OFFSET_BASIS_64;
    h = FNV_PRIME_64* h ^ x.first;
    h = FNV_PRIME_64* h ^ x.second;
    return h;
  }
};

// Bit tricks
constexpr int8_t deBruijnHash32[32] = 
        {0, 1, 28, 2, 29, 14, 24, 3, 30, 22, 20, 15, 25, 17, 4, 8, 31, 27, 13, 23, 21, 19, 16, 7, 26, 12, 18, 6, 11, 5, 10, 9};
constexpr int8_t log2ton(uint32_t val){
  return deBruijnHash32[static_cast<uint32_t>(val * 0x077CB531U) >> 27];
}

constexpr int8_t deBruijnHash64[64] =
        {0, 1, 59, 2, 60, 40, 54, 3, 61, 32, 49, 41, 55, 19, 35, 4, 62, 52, 30, 33, 50, 12, 14, 42, 56, 16, 27, 20, 36, 23, 44, 5, 63, 58, 39, 53, 31, 48, 18, 34, 51, 29, 11, 13, 15, 26, 22, 43, 57, 38, 47, 17, 28, 10, 25, 21, 37, 46, 9, 24, 45, 8, 7, 6};
constexpr int8_t log2ton(uint64_t val){
  return deBruijnHash64[static_cast<uint64_t>(val * 0x03F566ED27179461ULL) >> 58];
}

constexpr bool is2ton(uint32_t val){
  return val && !(val & (val - 1));
}
constexpr bool is2ton(uint64_t val){
  return val && !(val & (val - 1));
}

#endif // EQOPT_COMMON_HPP