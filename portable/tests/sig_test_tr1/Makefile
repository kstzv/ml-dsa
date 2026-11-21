CC ?= gcc

CFLAGS = -O3 -Wall -Wextra -std=c11 -ljansson \
		 -DUSERSPACE \
		 -DTEST \
          -I. \
          -I../../src \
          -I../../sha2_32_core \
          -I../../sha2_64_core \
          -I../../ports/userspace \
          -DML_DSA_FULL_MATRIX_BUFFER=0 \
          -DML_DSA_FULL_MATRIX_BUFFER_SIGN=0 \
          -DML_DSA_SIMD=1 \
          -DSHA2_32=1 \
          -DSHA2_64=1

SRC_CORE = \
	../../ports/userspace/ml_dsa_userspace_port.c \
    ../../src/keccak_1600.c \
    ../../src/ml_dsa.c \
    ../../src/ml_dsa_create_keys.c \
    ../../src/ml_dsa_get_verify_mem.c \
    ../../src/ml_dsa_matrix.c \
    ../../src/ml_dsa_sign_core.c \
    ../../src/ml_dsa_sign_verify_common.c \
    ../../src/ml_dsa_verify_core.c \
    ../../src/shake.c \
    ../../sha2_32_core/sha2_32.c \
    ../../sha2_32_core/sha2_32_core.c \
    ../../sha2_64_core/sha2_64.c \
    ../../sha2_64_core/sha2_64_core.c
    

TEST = test_external.c

TARGET = test_external

all:
	$(CC) $(CFLAGS) $(SRC_CORE) $(TEST) -lm -o $(TARGET)

clean:
	rm -f $(TARGET)
