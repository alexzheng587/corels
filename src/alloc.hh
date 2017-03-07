#pragma once
#include "utils.hh"

template <class T>  
    struct track_alloc { 
    typedef T value_type;
    track_alloc() noexcept {}
    template <class U> track_alloc (const track_alloc<U>&) noexcept {}
    T* allocate (size_t n) {
        printf("SZ N: %zu\n", n * sizeof(T));
        return static_cast<T*>(malloc(n*sizeof(T)));
    }   
    void deallocate (T* p, size_t n) {
        printf("DEALLOCATE\n");
        free(p);
    }   
};

template <class T>
struct pmap_alloc : track_alloc<T> {
    typedef T value_type;
    T* allocate (size_t n) {
        //printf("SZ N: %zu\n", n * sizeof(T));
        logger.addToPmapMemory(sizeof(T) * n); 
        return static_cast<T*>(malloc(n*sizeof(T)));
    }   
    void deallocate (T* p, size_t n) {
        logger.removeFromPmapMemory(sizeof(*p) * n); 
        free(p);
    }   
};

template <class T>
struct cache_alloc : track_alloc<T> {
    typedef T value_type;
    cache_alloc() noexcept {}
    template <class U> cache_alloc (const cache_alloc<U>&) noexcept {}
    T* allocate (size_t n) {
        logger.addToTreeMemory(n * sizeof(T));
        return static_cast<T*>(malloc(n*sizeof(T)));
    }
    void deallocate (T* p, size_t n) {
        logger.removeFromTreeMemory(n * sizeof(*p));
        free(p);
    }
};

template <class T>  
struct queue_alloc : track_alloc<T> {
    typedef T value_type;
    T* allocate (size_t n) {
        logger.addToQueueMemory(n * sizeof(T));
        return static_cast<T*>(malloc(n*sizeof(T)));
    }   
    void deallocate (T* p, size_t n) {
        logger.removeFromQueueMemory(n * sizeof(*p));
        free(p);
    }   
};
