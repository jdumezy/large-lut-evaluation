FROM ubuntu:24.04

RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates git cmake ninja-build clang libomp-dev \
    && rm -rf /var/lib/apt/lists/*

ARG BUILD_JOBS=2
ARG NATIVEOPT=ON
ARG BUILD_TYPE=Release
RUN git clone --depth 1 --branch v1.5.1 https://github.com/openfheorg/openfhe-development.git /opt/openfhe \
    && cmake -S /opt/openfhe -B /opt/openfhe/build -G Ninja \
       -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ \
       -DCMAKE_BUILD_TYPE=${BUILD_TYPE} -DBUILD_UNITTESTS=OFF -DBUILD_EXAMPLES=OFF \
       -DBUILD_BENCHMARKS=OFF -DNATIVE_SIZE=64 -DWITH_NATIVEOPT=${NATIVEOPT} \
    && cmake --build /opt/openfhe/build -j ${BUILD_JOBS} \
    && cmake --install /opt/openfhe/build && ldconfig \
    && rm -rf /opt/openfhe

WORKDIR /app
COPY CMakeLists.txt LICENSE ./
COPY include/ include/
COPY src/ src/
COPY tests/ tests/
RUN cmake -S . -B build -G Ninja -DCMAKE_CXX_COMPILER=clang++ \
       -DCMAKE_BUILD_TYPE=${BUILD_TYPE} -DNATIVEOPT=${NATIVEOPT} \
    && cmake --build build -j ${BUILD_JOBS} \
    && ctest --test-dir build --output-on-failure
CMD ["ctest", "--test-dir", "build", "--output-on-failure"]
