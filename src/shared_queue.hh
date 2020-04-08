#pragma once

#include "alloc.hh"
#include "queue.hh"
#include <queue>
#include <mutex>
#include "concurrentqueue.h"

typedef SharedQueue moodycamel::ConcurrentQueue<Queue*>;

class CSharedQueue
{
    //
    private:
        /* data */
        std::queue<Queue*> q_;
        size_t accs_;

    public:
        CSharedQueue(/* args */) {
          accs_ = 0;
        };
        ~CSharedQueue() {};

        inline bool empty();
        inline size_t size();
        inline Queue* pop();
        inline void push(Queue* entry);
        inline size_t n_acc() const;
};

inline bool CSharedQueue::empty() {
    return q_.empty();
}

inline size_t CSharedQueue::size() {
    return q_.size();
}

inline void CSharedQueue::push(Queue* entry) {
    q_.push(entry);
}

inline Queue* CSharedQueue::pop() {
    if(q_.empty()) {
        // TODO: change to NULL
        return (Queue*)0xDEADBEEF;
    }
    Queue* entry = q_.front();
    q_.pop();
    return entry;
}

inline size_t CSharedQueue::n_acc() const {
    return accs_;
}
*