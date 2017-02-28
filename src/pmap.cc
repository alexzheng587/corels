#include "pmap.hh"

PrefixPermutationMap::PrefixPermutationMap()
    : pmap(new std::unordered_map<prefix_key, std::pair<double, unsigned char*>, prefix_hash>) {}

CapturedPermutationMap::CapturedPermutationMap()
    : cmap(new std::unordered_map<captured_key, std::pair<std::vector<unsigned short>, double>, captured_hash>) {}
