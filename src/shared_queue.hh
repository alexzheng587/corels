#pragma once

#include "alloc.hh"
#include "queue.hh"
#include <queue>
#include <mutex>

typedef std::pair<Node*, tracking_vector<unsigned short, DataStruct::Tree> > internal_root;

class SharedQueue
{
    //
    private:
        /* data */
        std::mutex queue_lk_;
        std::queue<internal_root> q_;
        size_t accs_;

    public:
        SharedQueue(/* args */) {
          accs_ = 0;
        };
        ~SharedQueue();

        inline bool empty();
        inline size_t size();
        inline internal_root pop();
        inline void push(internal_root entry);
        inline void lock();
        inline void unlock();
        inline size_t n_acc() const;
};

inline void SharedQueue::lock() {
    ++accs_;
    queue_lk_.lock();
}

inline void SharedQueue::unlock() {
    queue_lk_.unlock();
}

inline bool SharedQueue::empty() {
    return q_.empty();
}

inline size_t SharedQueue::size() {
    return q_.size();
}

inline void SharedQueue::push(internal_root entry) {
    q_.push(entry);
}

inline internal_root SharedQueue::pop() {
    internal_root entry = q_.front();
    q_.pop();
    return entry;
}

inline size_t SharedQueue::n_acc() const {
    return accs_;
}
