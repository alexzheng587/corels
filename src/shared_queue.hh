#pragma once

#include "alloc.hh"
#include "queue.hh"
#include <queue>
#include <mutex>

class SharedQueue
{
    //
    private:
        /* data */
        std::mutex queue_lk_;
        std::queue<Queue*> q_;
        size_t accs_;

    public:
        SharedQueue(/* args */) {
          accs_ = 0;
        };
        ~SharedQueue();

        inline bool empty();
        inline size_t size();
        inline Queue* pop();
        inline void push(Queue* entry);
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

inline void SharedQueue::push(Queue* entry) {
    q_.push(entry);
}

inline Queue* SharedQueue::pop() {
    Queue* entry = q_.front();
    q_.pop();
    return entry;
}

inline size_t SharedQueue::n_acc() const {
    return accs_;
}
