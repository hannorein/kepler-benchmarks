#ifndef _KB_HELPERS_H_
#define _KB_HELPERS_H_

#include <math.h>

#define max(a, b) (((a) > (b)) ? (a) : (b))
#define min(a, b) (((a) < (b)) ? (a) : (b))

inline int sign(double x) { return (x > 0) - (x < 0); }

#endif
