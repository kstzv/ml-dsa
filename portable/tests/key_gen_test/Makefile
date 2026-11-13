CC ?= gcc

CFLAGS = -O3 -Wall -Wextra -std=c11 \
		 -DUSERSPACE \
          -I. \
          -I../../src \
          -I../../ports/userspace \
          -DML_DSA_FULL_MATRIX_BUFFER=0 \
          -DML_DSA_SIMD=1

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
    ../../src/shake.c
    

TEST = test.c

TARGET = test

all:
	$(CC) $(CFLAGS) $(SRC_CORE) $(TEST) -lm -o $(TARGET)

clean:
	rm -f $(TARGET)
