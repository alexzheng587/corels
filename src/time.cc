#include <string.h>
#include "time.hh"

void clear_time(struct time* t) {
    memset(t, 0, sizeof(*t));
}
