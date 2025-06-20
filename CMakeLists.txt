cmake_minimum_required(VERSION 3.10)
project(bliss_c VERSION 0.77.0 LANGUAGES C)

# Set C standard
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

# Compiler-specific optimizations
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra -Wpedantic")
    set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3 -DNDEBUG -march=native -mtune=native")
    set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -flto -fwhole-program")
    set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -funroll-loops -ffast-math")
    set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -finline-functions -finline-limit=1000")
    
    # IPO/LTO support
    if(CMAKE_BUILD_TYPE STREQUAL "Release")
        include(CheckIPOSupported)
        check_ipo_supported(RESULT ipo_supported OUTPUT error)
        if(ipo_supported)
            message(STATUS "IPO / LTO enabled")
            set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
        else()
            message(STATUS "IPO / LTO not supported: <${error}>")
        endif()
    endif()
    
    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -g -O0 -DDEBUG")
    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -fsanitize=address -fsanitize=undefined")
    set(CMAKE_EXE_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS_DEBUG} -fsanitize=address -fsanitize=undefined")
endif()

if(MSVC)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /W4")
    set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} /O2 /DNDEBUG /GL")
    set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} /LTCG")
endif()

# Build options
option(BUILD_SHARED_LIBS "Build shared libraries" OFF)  # Changed to OFF
option(BUILD_TESTS "Build test programs" ON)
option(BUILD_BENCHMARKS "Build benchmark programs" ON)
option(ENABLE_PROFILING "Enable profiling support" OFF)

if(ENABLE_PROFILING)
    add_definitions(-DBLISS_PROFILING)
    if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pg")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
    endif()
endif()

# Find required packages
find_package(PkgConfig)

# Math library
find_library(MATH_LIBRARY m)

# Source files
# Update BLISS_SOURCES
set(BLISS_SOURCES
    src/bliss_core.c
    src/bliss_io.c
    src/bliss_partition.c      # NEW
    src/bliss_refinement.c     # NEW
    src/bliss_unified.c
    src/orbit.c
)

set(BLISS_HEADERS
    include/bliss.h
)

# Create the library
add_library(bliss ${BLISS_SOURCES} ${BLISS_HEADERS})

# Set library properties
set_target_properties(bliss PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
    PUBLIC_HEADER "${BLISS_HEADERS}"
)

# Link math library if found
if(MATH_LIBRARY)
    target_link_libraries(bliss ${MATH_LIBRARY})
endif()

# Include directories
target_include_directories(bliss PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

# Compiler-specific optimizations for the library
if(CMAKE_BUILD_TYPE STREQUAL "Release")
    if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(bliss PRIVATE
            -fomit-frame-pointer
            -fstrict-aliasing
            -fno-exceptions
            -fvisibility=hidden
        )
        
        # Profile-guided optimization (commented out - requires two-stage build)
        # target_compile_options(bliss PRIVATE -fprofile-generate)
        # target_link_options(bliss PRIVATE -fprofile-generate)
    endif()
endif()

# Tests
# Tests
if(BUILD_TESTS)
    enable_testing()
    
    # ===================================================================
    # NEW COMPREHENSIVE TEST SUITE (Primary Tests)
    # ===================================================================
    
    # Main comprehensive test executable with all new algorithms
    add_executable(test_bliss_complete tests/test_bliss_complete.c)
    target_link_libraries(test_bliss_complete bliss)
    if(WIN32)
        target_link_libraries(test_bliss_complete -static-libgcc -static-libstdc++)
    endif()
    
    # Comprehensive test categories using the new implementation
    add_test(NAME CompleteBlissTests COMMAND test_bliss_complete all)
    add_test(NAME NewBasicTests COMMAND test_bliss_complete basic)
    add_test(NAME NewAutomorphismTests COMMAND test_bliss_complete auto)
    add_test(NAME NewPerformanceTests COMMAND test_bliss_complete perf)
    #add_test(NAME InteractiveTestMode COMMAND test_bliss_complete interactive)
    
    # ===================================================================
    # LEGACY TESTS (Keep for compatibility and regression testing)
    # ===================================================================
    
    # Original basic functionality tests
    add_executable(test_basic tests/test_basic.c)
    target_link_libraries(test_basic bliss)
    if(WIN32)
        target_link_libraries(test_basic -static-libgcc -static-libstdc++)
    endif()
    add_test(NAME LegacyBasicTests COMMAND test_basic)
    
    # Original graph creation tests
    add_executable(test_graph_creation tests/test_graph_creation.c)
    target_link_libraries(test_graph_creation bliss)
    if(WIN32)
        target_link_libraries(test_graph_creation -static-libgcc -static-libstdc++)
    endif()
    add_test(NAME LegacyGraphCreationTests COMMAND test_graph_creation)
    
    # Original automorphism tests (for comparison)
    add_executable(test_automorphisms tests/test_automorphisms.c)
    target_link_libraries(test_automorphisms bliss)
    if(WIN32)
        target_link_libraries(test_automorphisms -static-libgcc -static-libstdc++)
    endif()
    add_test(NAME LegacyAutomorphismTests COMMAND test_automorphisms)
    
    # I/O tests
    add_executable(test_io tests/test_io.c)
    target_link_libraries(test_io bliss)
    if(WIN32)
        target_link_libraries(test_io -static-libgcc -static-libstdc++)
    endif()
    add_test(NAME IOTests COMMAND test_io)
    
    # Canonical labeling tests
    add_executable(test_canonical tests/test_canonical.c)
    target_link_libraries(test_canonical bliss)
    if(WIN32)
        target_link_libraries(test_canonical -static-libgcc -static-libstdc++)
    endif()
    add_test(NAME LegacyCanonicalTests COMMAND test_canonical)
    
    # Performance comparison tests
    add_executable(test_performance tests/test_performance.c)
    target_link_libraries(test_performance bliss)
    if(WIN32)
        target_link_libraries(test_performance -static-libgcc -static-libstdc++)
    endif()
    add_test(NAME LegacyPerformanceTests COMMAND test_performance)
    
    # Memory tests
    add_executable(test_memory tests/test_memory.c)
    target_link_libraries(test_memory bliss)
    if(WIN32)
        target_link_libraries(test_memory -static-libgcc -static-libstdc++)
    endif()
    add_test(NAME MemoryTests COMMAND test_memory)
    
    # Validation tests
    add_executable(test_validation tests/test_validation.c)
    target_link_libraries(test_validation bliss)
    if(WIN32)
        target_link_libraries(test_validation -static-libgcc -static-libstdc++)
    endif()
    add_test(NAME ValidationTests COMMAND test_validation)
    
     # LP foolding tests
    add_executable(test_lpfold tests/test_lpfold.c)
    target_link_libraries(test_lpfold bliss)
    if(WIN32)
        target_link_libraries(test_lpfold -static-libgcc -static-libstdc++)
    endif()
    add_test(NAME LPFoldingTests COMMAND test_lpfold)

    # Classic graph tests (from original bliss distribution)
    add_executable(test_classic_graphs tests/test_classic_graphs.c)
    target_link_libraries(test_classic_graphs bliss)
    if(WIN32)
        target_link_libraries(test_classic_graphs -static-libgcc -static-libstdc++)
    endif()
    add_test(NAME ClassicGraphTests COMMAND test_classic_graphs)
    
    # ===================================================================
    # TEST PROPERTIES AND TIMEOUTS
    # ===================================================================
    
    # Set properties for new comprehensive tests (longer timeout for complex tests)
    set_tests_properties(
        CompleteBlissTests NewBasicTests NewAutomorphismTests NewPerformanceTests
        PROPERTIES
        TIMEOUT 600  # 10 minutes for comprehensive tests
    )
    
    # Set properties for legacy tests (shorter timeout)
    set_tests_properties(
        LegacyBasicTests LegacyGraphCreationTests LegacyAutomorphismTests IOTests 
        LegacyCanonicalTests LegacyPerformanceTests MemoryTests ValidationTests
        ClassicGraphTests
        PROPERTIES
        TIMEOUT 300  # 5 minutes for legacy tests
    )
    
    # ===================================================================
    # CUSTOM TEST TARGETS FOR CONVENIENCE
    # ===================================================================
    
    # Add custom targets for easy testing
    add_custom_target(test_new
        COMMAND ${CMAKE_CTEST_COMMAND} -R "^New" --verbose
        DEPENDS test_bliss_complete
        COMMENT "Running new comprehensive tests only"
    )
    
    add_custom_target(test_legacy  
        COMMAND ${CMAKE_CTEST_COMMAND} -R "^Legacy" --verbose
        DEPENDS test_basic test_graph_creation test_automorphisms test_canonical test_performance
        COMMENT "Running legacy tests only"
    )
    
    add_custom_target(test_quick
        COMMAND ${CMAKE_CTEST_COMMAND} -R "NewBasicTests|LegacyBasicTests" --verbose
        DEPENDS test_bliss_complete test_basic
        COMMENT "Running quick basic tests"
    )
    
    add_custom_target(test_automorphism_comparison
        COMMAND ${CMAKE_CTEST_COMMAND} -R "AutomorphismTests" --verbose
        DEPENDS test_bliss_complete test_automorphisms
        COMMENT "Running automorphism tests (new vs legacy)"
    )
    
    # ===================================================================
    # TEST INFORMATION SUMMARY
    # ===================================================================
    
    message(STATUS "")
    message(STATUS "Test Configuration:")
    message(STATUS "  New comprehensive tests: test_bliss_complete")
    message(STATUS "  Legacy compatibility tests: individual test_* executables")
    message(STATUS "  Available test targets:")
    message(STATUS "    make test          - Run all tests")
    message(STATUS "    make test_new      - Run new comprehensive tests only") 
    message(STATUS "    make test_legacy   - Run legacy tests only")
    message(STATUS "    make test_quick    - Run basic tests only")
    message(STATUS "    make test_automorphism_comparison - Compare new vs legacy automorphism tests")
    message(STATUS "")
    message(STATUS "  Manual test execution:")
    message(STATUS "    ./test_bliss_complete all     - Full comprehensive test suite")
    message(STATUS "    ./test_bliss_complete basic   - Basic functionality tests")
    message(STATUS "    ./test_bliss_complete auto    - Automorphism detection tests")
    message(STATUS "    ./test_bliss_complete perf    - Performance tests")
    message(STATUS "    ./test_bliss_complete interactive - Interactive test mode")
    message(STATUS "")
    
endif()

# Benchmarks
if(BUILD_BENCHMARKS)
    add_executable(benchmark_main benchmarks/benchmark_main.c)
    target_link_libraries(benchmark_main bliss)
    
    add_executable(benchmark_comparison benchmarks/benchmark_comparison.c)
    target_link_libraries(benchmark_comparison bliss)
endif()

# Command-line tool
add_executable(bliss_tool tools/bliss_main.c)
target_link_libraries(bliss_tool bliss)
set_target_properties(bliss_tool PROPERTIES OUTPUT_NAME bliss)

# Installation
include(GNUInstallDirs)

install(TARGETS bliss bliss_tool
    EXPORT blissTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# Install export targets
install(EXPORT blissTargets
    FILE blissTargets.cmake
    NAMESPACE bliss::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/bliss
)

# Create config file
include(CMakePackageConfigHelpers)
configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/blissConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/blissConfig.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/bliss
)

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/blissConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/blissConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/blissConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/bliss
)

# pkg-config file
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/bliss.pc.in")
    configure_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/bliss.pc.in"
        "${CMAKE_CURRENT_BINARY_DIR}/bliss.pc"
        @ONLY
    )

    install(FILES "${CMAKE_CURRENT_BINARY_DIR}/bliss.pc"
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig
    )
else()
    message(STATUS "bliss.pc.in not found - pkg-config file will not be generated")
endif()

# Documentation
#find_package(Doxygen)
if(DOXYGEN_FOUND)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in 
                   ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
    add_custom_target(doc
        ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen" VERBATIM
    )
endif()

# Print configuration summary
message(STATUS "")
message(STATUS "Configuration Summary:")
message(STATUS "  Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "  C Compiler: ${CMAKE_C_COMPILER_ID} ${CMAKE_C_COMPILER_VERSION}")
message(STATUS "  C Flags: ${CMAKE_C_FLAGS}")
if(CMAKE_BUILD_TYPE STREQUAL "Release")
    message(STATUS "  C Release Flags: ${CMAKE_C_FLAGS_RELEASE}")
    message(STATUS "  IPO/LTO: ${CMAKE_INTERPROCEDURAL_OPTIMIZATION}")
endif()
message(STATUS "  Build shared libs: ${BUILD_SHARED_LIBS}")
message(STATUS "  Build tests: ${BUILD_TESTS}")
message(STATUS "  Build benchmarks: ${BUILD_BENCHMARKS}")
message(STATUS "  Install prefix: ${CMAKE_INSTALL_PREFIX}")
message(STATUS "")