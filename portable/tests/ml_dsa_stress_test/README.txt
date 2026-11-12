ML-DSA stress test update (based on uploaded sources, 2026-11-12)

Copy portable/tests/ml_dsa_stress_test/{ml_dsa_stress_test.c,Makefile} into your repo.
Apply ml_dsa_fixes.patch from the repository root:
    patch -p1 < ml_dsa_fixes.patch

Patch adjusts 4 files:
  portable/src/ml_dsa.c                   fixed pool uses cached rho when full matrix disabled
  portable/src/ml_dsa.h                   correct independent include guard
  portable/src/ml_dsa_core_header.h       const-correct unpack_h declaration
  portable/src/ml_dsa_verify_core.c       const-correct unpack_h implementation

Build and run:
    cd portable/tests/ml_dsa_stress_test
    make clean && make
    make run LEVEL=all ITERATIONS=1000
    make run-asan LEVEL=all ITERATIONS=100

Note: patch is based on the CURRENT uploaded ml_dsa.c (destructor_verify_pool has one argument).
