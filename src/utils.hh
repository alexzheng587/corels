#pragma once

#include "alloc.hh"
/* 
 * Prints the final rulelist that CORELS returns.
 * rulelist -- rule ids of optimal rulelist
 * preds -- corresponding predictions of rules (+ default prediction)
 */
void print_final_rulelist(const tracking_vector<unsigned short, DataStruct::Tree> &rulelist,
                          const tracking_vector<bool, DataStruct::Tree> &preds,
                          const bool latex_out,
                          const rule_t rules[],
                          const rule_t labels[],
                          char fname[]);

void print_machine_info();

// From: https://stackoverflow.com/questions/32346861/how-to-write-a-generic-print-for-vector-of-vectors-using-c
template <typename T>
void printVector(const T &t) {
    std::copy(t.cbegin(), t.cend(), std::ostream_iterator<typename T::value_type>(std::cout, ", "));
    std::cout << "\n";
}
