@echo off
echo Simple Bliss C Library Build
echo ===========================

REM Create build directory
if not exist build mkdir build
cd build

REM Simple build without CMake (for testing)
echo Building library...
gcc -c -O3 -Wall -Wextra -std=c11 ../bliss_core.c -o bliss_core.o
gcc -c -O3 -Wall -Wextra -std=c11 ../bliss_io.c -o bliss_io.o

echo Creating static library...
ar rcs libbliss.a bliss_core.o bliss_io.o

echo Building command-line tool...
gcc -O3 -Wall -Wextra -std=c11 ../tools/bliss_main.c libbliss.a -o bliss.exe -lm

echo Building basic test...
gcc -O3 -Wall -Wextra -std=c11 ../tests/test_basic.c libbliss.a -o test_basic.exe -lm

echo.
echo Build completed!
echo Run: test_basic.exe
echo Run: bliss.exe --help

cd ..