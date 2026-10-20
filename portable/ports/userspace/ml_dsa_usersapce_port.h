// Porting ML-DSA in userspace
#ifndef ML_DSA_USERSPACE_PORT_H
#define ML_DSA_USERSPACE_PORT_H

// @>@@>@@@>@@@@>@@@@@>@@@@@@>@@@@@@@>@@@@@@@@>@@@@@@@@@>@@@@@@@@@@>@@@@@@@@@@@>@@@@@@@@@@@@@>@@@@@@@@@@@@@
// Including libraries from userspace OS based Linux kernel
#include <stdlib.h>   
#include <string.h>   
#include <sys/random.h>
#include <stdint.h>
#include <stddef.h>
#include <errno.h>


// ^>^^>^^^>^^^^>^^^^^>^^^^^^>^^^^^^^>^^^^^^^^>^^^^^^^^^>^^^^^^^^^^>^^^^^^^^^^^>^^^^^^^^^^^^>^^^^^^^^^^^^^
// Prototypes of functions that depend on a this platform

// Fixed-width integer aliases (kernel-style, userspace-safe)
typedef uint8_t  u8;
typedef uint16_t u16;
typedef uint32_t u32;
typedef uint64_t u64;
typedef int64_t s64;
typedef int32_t s32;
typedef int16_t s16;
typedef int8_t s8;

#endif // ML_DSA_USERSPACE_PORT_H
