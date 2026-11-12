# Build from /portable/tests/ml_dsa_stress_test
CC ?= cc
SRC_DIR := ../../src
PORT_DIR := ../../ports/userspace

TARGET := ml_dsa_stress_test
ASAN_TARGET := ml_dsa_stress_test_asan

CPPFLAGS += -DUSERSPACE -DML_DSA_SIMD=1 -DML_DSA_FULL_MATRIX_BUFFER=0 \
            -I$(SRC_DIR) -I$(PORT_DIR)
CFLAGS ?= -std=c11 -O2 -g -Wall -Wextra

CORE_SRC := $(SRC_DIR)/ml_dsa.c \
            $(SRC_DIR)/ml_dsa_create_keys.c \
            $(SRC_DIR)/ml_dsa_get_verify_mem.c \
            $(SRC_DIR)/ml_dsa_matrix.c \
            $(SRC_DIR)/ml_dsa_sign_core.c \
            $(SRC_DIR)/ml_dsa_sign_verify_common.c \
            $(SRC_DIR)/ml_dsa_verify_core.c \
            $(SRC_DIR)/shake.c \
            $(SRC_DIR)/keccak_1600.c
PORT_SRC := $(PORT_DIR)/ml_dsa_usersapce_port.c
TEST_SRC := ml_dsa_stress_test.c
SOURCES := $(TEST_SRC) $(CORE_SRC) $(PORT_SRC)

LEVEL ?= all
ITERATIONS ?= 100
SANFLAGS := -O1 -g -fno-omit-frame-pointer -fsanitize=address,undefined

.PHONY: all run asan run-asan clean
all: $(TARGET)

$(TARGET): $(SOURCES) $(SRC_DIR)/ml_dsa.h $(SRC_DIR)/ml_dsa_core_header.h $(SRC_DIR)/shake.h $(PORT_DIR)/ml_dsa_userspace_port.h
	$(CC) $(CPPFLAGS) $(CFLAGS) $(SOURCES) -o $@ $(LDFLAGS) $(LDLIBS)

run: $(TARGET)
	./$(TARGET) $(LEVEL) $(ITERATIONS)

asan: $(ASAN_TARGET)

$(ASAN_TARGET): $(SOURCES) $(SRC_DIR)/ml_dsa.h $(SRC_DIR)/ml_dsa_core_header.h $(SRC_DIR)/shake.h $(PORT_DIR)/ml_dsa_userspace_port.h
	$(CC) $(CPPFLAGS) -std=c11 -Wall -Wextra $(SANFLAGS) $(SOURCES) -o $@ $(SANFLAGS) $(LDFLAGS) $(LDLIBS)

run-asan: $(ASAN_TARGET)
	./$(ASAN_TARGET) $(LEVEL) $(ITERATIONS)

clean:
	rm -f $(TARGET) $(ASAN_TARGET)
