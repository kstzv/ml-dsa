// SPDX-License-Identifier: GPL-2.0 OR MIT
// Copyright (c) 2026 K.S.Zavertailo
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

// For atomic operations, 1 - slot is free, 0 - slot is busy
bool ml_dsa_try_acquire_slot(ml_dsa_atomic_t *slot)
{
    int expected = 1;
    return atomic_compare_exchange_strong_explicit(slot, &expected, 0, memory_order_acquire, memory_order_relaxed);
}

void ml_dsa_release_slot(ml_dsa_atomic_t *slot) { atomic_store_explicit(slot, 1, memory_order_release); }

void ml_dsa_atomic_init_slot(ml_dsa_atomic_t *slot) { atomic_init(slot, 1); }



// Error handling
struct ml_dsa_keys *err_create_keys(const int err)
{
	(void)err;
	return NULL;
}

int err_ml_dsa_sign(const int err) { return err; }

int err_ml_dsa_verify(const int err) { return err; }

struct ml_dsa_verify_ctx *err_ml_dsa_alloc_verify_pool(const int err)
{
	(void)err;
	return NULL;
}

int err_ml_dsa_get_pk(const int err) { return err; }

int err_ml_dsa_get_sk(const int err) { return err; }

struct ml_dsa_keys *err_ml_dsa_create_keys_object_from_imported_sk(const int err)
{
	(void)err;
	return NULL;
}
