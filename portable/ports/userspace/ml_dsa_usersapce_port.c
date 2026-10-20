#include "ml_dsa_userspace_port.h"

// Functions for allocations
void *ml_dsa_alloc(size_t size_alloc) { return malloc(size_alloc); }
void ml_dsa_free(void *ptr) { free(ptr); }

// Local function entropy ML-DSA (Using in Linux OS userspace)
int ml_dsa_entropy(void *buf, size_t len)
{
	ssize_t ret = 0;
	size_t offset = 0;
	while(offset < len)
	{
		ret = getrandom((u8 *)buf + offset, len - offset, 0);
		if(ret < 0) 
		{
			if(errno == EINTR) { continue; }
			return -1;
		}
		offset += ret;
	}
	return 0;
}

// Secure memory wipe function.
//
// Uses a volatile pointer to prevent compiler optimizations
// that could remove the memory clearing operation.
//
// NOTE:
//   Effectiveness depends on compiler behavior.
//   For high-assurance environments, consider platform-specific primitives
//   (e.g., explicit_bzero, memset_s, or kernel equivalents).
void ml_dsa_memzero(void *ptr, size_t len)
{
	if(!ptr) { return; }
	
	volatile unsigned char *inside_ptr = (volatile unsigned char *)ptr;

	for(size_t i = 0; i < len; i++)
	{
		inside_ptr[i] = 0;
	}
}

// Error handling
struct ml_dsa_keys *err_create_keys(const u16 err)
{
	(void)err;
	return NULL;
}
