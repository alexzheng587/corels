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
        inline size_t n_acc() const;
};

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
    if(q_.empty()) {
        // TODO: change to NULL
        return (Queue*)0xDEADBEEF;
    }
    Queue* entry = q_.front();
    q_.pop();
    return entry;
}

inline size_t SharedQueue::n_acc() const {
    return accs_;
}
