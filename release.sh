cmake -S . -B build-static \
  -DCMAKE_BUILD_TYPE=Release \
  -DBEAUTIZE_OPENMP=OFF \
  -DBEAUTIZE_BUILD_TESTS=OFF \
  -DBLA_STATIC=ON \
  -DBLA_VENDOR=Generic \
  -DCMAKE_LIBRARY_PATH="/usr/lib/x86_64-linux-gnu/blas;/usr/lib/x86_64-linux-gnu/lapack" \
  -DCMAKE_C_FLAGS="-ffunction-sections -fdata-sections" \
  -DCMAKE_CXX_FLAGS="-ffunction-sections -fdata-sections" \
  -DCMAKE_Fortran_FLAGS="-ffunction-sections -fdata-sections" \
  -DCMAKE_EXE_LINKER_FLAGS="-static -Wl,--gc-sections -s"

cmake --build build-static -j"$(nproc)"

