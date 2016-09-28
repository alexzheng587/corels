#include <sys/time.h>

inline double timestamp() {
    struct timeval now;
    gettimeofday(&now, 0);
    return now.tv_sec + now.tv_usec * 0.000001;
}
