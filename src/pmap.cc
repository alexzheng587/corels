#include "pmap.hh"

PrefixPermutationMap::PrefixPermutationMap()
    : pmap(new std::unordered_map<prefix_key, prefix_val, prefix_hash, prefix_eq, std::scoped_allocator_adaptor<tracking_allocator<std::pair<const prefix_key, prefix_val> > > >) {}

CapturedPermutationMap::CapturedPermutationMap()
    : cmap(new std::unordered_map<captured_key, std::pair<std::vector<unsigned short>, double>, captured_hash>) {}
