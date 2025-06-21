#!/bin/bash

# Build script for Bliss C library
# Works on Linux, macOS, and Windows (MSYS2/MinGW)

set -e  # Exit on any error

echo "Building Bliss C Library"
echo "======================="

# Create build directory
if [ ! -d "build" ]; then
    mkdir build
fi

cd build

# Configure with CMake
echo "Configuring with CMake..."
cmake -DCMAKE_BUILD_TYPE=Release \
      -DBUILD_TESTS=ON \
      -DBUILD_BENCHMARKS=ON \
      ..

# Build
echo "Building..."
cmake --build . --config Release

# Check if we're on Windows (MSYS2/MinGW)
if [[ "$OSTYPE" == "msys" || "$OSTYPE" == "cygwin" ]]; then
    echo "Windows build completed"
    echo "Executables are in build/Release/ or build/"
else
    echo "Unix build completed"
    
    # Run tests if available
    if command -v ctest &> /dev/null; then
        echo "Running tests..."
        ctest --verbose
    fi
fi

echo ""
echo "Build completed successfully!"
echo "To run benchmarks: ./benchmark_main"
echo "To run command-line tool: ./bliss --help"