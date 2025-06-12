# Bliss C Implementation

A high-performance C translation of the bliss graph automorphism and canonical labeling library, originally written in C++ by Tommi Junttila and Petteri Kaski.

## Overview

This library computes automorphism groups and canonical labelings of graphs using sophisticated partition refinement and backtracking search algorithms. It's designed for applications requiring graph isomorphism testing, automorphism group analysis, and canonical graph representations.

## Features

### Core Functionality
- **Graph Automorphism Detection**: Find generators of the automorphism group
- **Canonical Labeling**: Compute canonical forms for isomorphism testing
- **Multiple Graph Types**: Support for directed and undirected graphs
- **Vertex Coloring**: Handle graphs with colored vertices
- **Performance Optimizations**: Memory-efficient algorithms with multiple heuristics

### Algorithm Implementation
- **Partition Refinement**: 1-dimensional Weisfeiler-Lehman refinement
- **Backtracking Search**: Sophisticated search tree with pruning
- **Orbit-based Pruning**: Skip equivalent search branches
- **Canonical Path Tracking**: Efficient canonical form computation
- **Multiple Splitting Heuristics**: Optimized for different graph types

### Supported Graph Formats
- **DIMACS Format**: Standard graph competition format
- **DOT Format**: Graphviz visualization format
- **Programmatic API**: Direct graph construction in code

## Installation

### Prerequisites
- C11-compatible compiler (GCC 4.9+, Clang 3.5+, MSVC 2015+)
- CMake 3.10 or higher
- Standard math library

### Building from Source

```bash
git clone <repository-url>
cd bliss-c
mkdir build && cd build
cmake ..
make -j$(nproc)
```

### Build Options

```bash
# Release build with optimizations
cmake -DCMAKE_BUILD_TYPE=Release ..

# Debug build with sanitizers
cmake -DCMAKE_BUILD_TYPE=Debug ..

# Build with tests
cmake -DBUILD_TESTS=ON ..

# Build with benchmarks
cmake -DBUILD_BENCHMARKS=ON ..

# Enable profiling support
cmake -DENABLE_PROFILING=ON ..
```

### Installation

```bash
make install
```

This installs:
- Library: `libbliss.a` (static) or `libbliss.so` (shared)
- Headers: `bliss.h`
- Command-line tool: `bliss`
- pkg-config file: `bliss.pc`

## Quick Start

### Basic Usage

```c
#include "bliss.h"
#include <stdio.h>

int main() {
    // Create a 4-vertex cycle graph
    bliss_graph_t *graph = bliss_new(4);
    bliss_add_edge(graph, 0, 1);
    bliss_add_edge(graph, 1, 2);
    bliss_add_edge(graph, 2, 3);
    bliss_add_edge(graph, 3, 0);
    
    // Find automorphisms
    bliss_stats_t *stats = bliss_stats_new();
    bliss_find_automorphisms(graph, stats, NULL, NULL);
    
    printf("Found %lu generators\n", bliss_stats_get_nof_generators(stats));
    printf("Group size approximately: %.0f\n", 
           bliss_stats_get_group_size_approx(stats));
    
    // Cleanup
    bliss_stats_release(stats);
    bliss_release(graph);
    return 0;
}
```

### Automorphism Callback

```c
void automorphism_callback(void *user_param, unsigned int n, const unsigned int *aut) {
    printf("Automorphism: ");
    for (unsigned int i = 0; i < n; i++) {
        printf("%u->%u ", i, aut[i]);
    }
    printf("\n");
}

// Use the callback
bliss_find_automorphisms(graph, stats, automorphism_callback, NULL);
```

### Graph Comparison

```c
// Create two graphs
bliss_graph_t *g1 = bliss_create_cycle_graph(5);
bliss_graph_t *g2 = bliss_create_path_graph(5);

// Compare using canonical forms
bliss_stats_t *stats1 = bliss_stats_new();
bliss_stats_t *stats2 = bliss_stats_new();

const unsigned int *canon1 = bliss_find_canonical_labeling(g1, stats1, NULL, NULL);
const unsigned int *canon2 = bliss_find_canonical_labeling(g2, stats2, NULL, NULL);

int comparison = bliss_cmp(g1, g2);
printf("Graphs are %s\n", comparison == 0 ? "isomorphic" : "different");
```

## API Reference

### Graph Construction

```c
// Create graphs
bliss_graph_t *bliss_new(unsigned int num_vertices);
bliss_graph_t *bliss_new_directed(unsigned int num_vertices);
void bliss_release(bliss_graph_t *graph);

// Add structure
void bliss_add_vertex(bliss_graph_t *graph, unsigned int color);
void bliss_add_edge(bliss_graph_t *graph, unsigned int v1, unsigned int v2);
void bliss_change_color(bliss_graph_t *graph, unsigned int vertex, unsigned int color);

// Bulk operations
void bliss_add_edges_bulk(bliss_graph_t *graph, 
                          const unsigned int *edge_list, 
                          unsigned int num_edges);
```

### Graph Properties

```c
// Basic properties
unsigned int bliss_get_nof_vertices(const bliss_graph_t *graph);
bool bliss_is_directed(const bliss_graph_t *graph);
unsigned int bliss_get_color(const bliss_graph_t *graph, unsigned int vertex);
uint32_t bliss_get_hash(bliss_graph_t *graph);

// Structural queries
unsigned int bliss_get_degree(const bliss_graph_t *graph, unsigned int vertex);
bool bliss_has_edge(const bliss_graph_t *graph, unsigned int v1, unsigned int v2);
const unsigned int *bliss_get_adjacency_list(const bliss_graph_t *graph, 
                                             unsigned int vertex, 
                                             unsigned int *list_size);
```

### Core Algorithms

```c
// Automorphism detection
void bliss_find_automorphisms(bliss_graph_t *graph, 
                              bliss_stats_t *stats,
                              bliss_automorphism_hook_t hook,
                              void *hook_user_param);

// Canonical labeling
const unsigned int *bliss_find_canonical_labeling(bliss_graph_t *graph,
                                                   bliss_stats_t *stats,
                                                   bliss_automorphism_hook_t hook,
                                                   void *hook_user_param);

// Graph comparison
int bliss_cmp(bliss_graph_t *g1, bliss_graph_t *g2);
bool bliss_is_automorphism(const bliss_graph_t *graph, const unsigned int *perm);
```

### Configuration

```c
// Splitting heuristics
typedef enum {
    BLISS_SH_F,      // First non-singleton cell
    BLISS_SH_FL,     // First largest non-singleton cell  
    BLISS_SH_FS,     // First smallest non-singleton cell
    BLISS_SH_FM,     // First maximally connected non-singleton cell
    BLISS_SH_FLM,    // First largest maximally connected non-singleton cell
    BLISS_SH_FSM     // First smallest maximally connected non-singleton cell
} bliss_splitting_heuristic_t;

void bliss_set_splitting_heuristic(bliss_graph_t *graph, bliss_splitting_heuristic_t h);
void bliss_set_component_recursion(bliss_graph_t *graph, bool enabled);
void bliss_set_failure_recording(bliss_graph_t *graph, bool enabled);
```

### Statistics

```c
// Create and manage statistics
bliss_stats_t *bliss_stats_new(void);
void bliss_stats_release(bliss_stats_t *stats);

// Access results
unsigned long bliss_stats_get_nof_generators(const bliss_stats_t *stats);
unsigned long bliss_stats_get_nof_nodes(const bliss_stats_t *stats);
double bliss_stats_get_group_size_approx(const bliss_stats_t *stats);
const char *bliss_stats_get_group_size_string(const bliss_stats_t *stats);
```

### I/O Operations

```c
// File I/O
bliss_graph_t *bliss_read_dimacs(FILE *fp);
void bliss_write_dimacs(const bliss_graph_t *graph, FILE *fp);
void bliss_write_dot(const bliss_graph_t *graph, FILE *fp);

// Convenience functions
bliss_graph_t *bliss_read_graph(const char *filename);
bool bliss_write_graph(const bliss_graph_t *graph, const char *filename);
bool bliss_write_dot_file(const bliss_graph_t *graph, const char *filename);
```

### Graph Generators

```c
// Standard graph families
bliss_graph_t *bliss_create_complete_graph(unsigned int n);
bliss_graph_t *bliss_create_cycle_graph(unsigned int n);
bliss_graph_t *bliss_create_path_graph(unsigned int n);
bliss_graph_t *bliss_create_star_graph(unsigned int n);
bliss_graph_t *bliss_create_petersen_graph(void);
bliss_graph_t *bliss_create_random_graph(unsigned int n, double edge_prob, unsigned int seed);
```

## Command Line Tool

The `bliss` command-line tool provides access to core functionality:

```bash
# Find automorphisms of a graph
bliss input.dimacs

# Generate canonical labeling
bliss -c input.dimacs

# Use specific heuristic
bliss -h fl input.dimacs

# Output in DOT format
bliss -d input.dimacs > output.dot

# Verbose output with statistics
bliss -v input.dimacs
```

### Command Line Options

- `-c`: Compute canonical labeling only
- `-h <heuristic>`: Set splitting heuristic (f, fl, fs, fm, flm, fsm)
- `-d`: Output in DOT format
- `-v`: Verbose output with detailed statistics
- `-t <seconds>`: Set time limit
- `--help`: Show usage information

## Testing

### Running Tests

```bash
# Run all tests
make test

# Run specific test categories
./test_bliss_comprehensive basic
./test_bliss_comprehensive auto
./test_bliss_comprehensive perf

# Interactive testing mode
./test_bliss_comprehensive interactive
```

### Test Categories

- **Basic Tests**: Graph creation, edge addition, property queries
- **Automorphism Tests**: Algorithm correctness on known graphs
- **Performance Tests**: Timing and memory usage analysis
- **I/O Tests**: File format reading and writing
- **Memory Tests**: Leak detection and stress testing

### Benchmark Graphs

The test suite includes standard benchmark graphs:

- Complete graphs (K_n)
- Cycle graphs (C_n)  
- Path graphs (P_n)
- Petersen graph
- Random graphs with various densities
- Highly symmetric graphs for stress testing

## Performance

### Optimization Features

- **Compiler Optimizations**: Aggressive optimization flags for release builds
- **Memory Management**: Efficient allocation strategies
- **Algorithm Selection**: Adaptive heuristics based on graph properties
- **Pruning Techniques**: Multiple levels of search space reduction

### Performance Characteristics

| Graph Type | Vertices | Typical Performance |
|------------|----------|-------------------|
| Sparse graphs | 1000 | < 1 second |
| Dense graphs | 100 | < 1 second |
| Highly symmetric | 50 | Variable (seconds to minutes) |
| Random graphs | 500 | < 10 seconds |

### Memory Usage

Memory usage scales approximately as O(n^2) for dense graphs and O(n + m) for sparse graphs, where n is vertices and m is edges.

## Implementation Details

### Algorithm Overview

1. **Initial Partition**: Vertices grouped by color
2. **Partition Refinement**: Iterative refinement using neighbor signatures  
3. **Target Cell Selection**: Choose cell to split using heuristic
4. **Individualization**: Create singleton cells for branching
5. **Recursive Search**: Depth-first search with backtracking
6. **Pruning**: Orbit-based and canonical path pruning
7. **Generator Collection**: Store automorphism generators

### Data Structures

- **Graphs**: Adjacency list representation with color arrays
- **Partitions**: Cell-based vertex partitioning with fast lookups
- **Search Tree**: Explicit tree nodes with partition snapshots
- **Statistics**: Comprehensive algorithm performance metrics

### Optimization Techniques

- **Invariant Computation**: Multiple graph invariants for efficiency
- **Memory Pools**: Reduced allocation overhead
- **Canonical Ordering**: Efficient canonical form computation
- **Component Recursion**: Handle disconnected graphs efficiently

## Troubleshooting

### Common Issues

**Build Errors**
- Ensure C11 compiler support
- Check CMake version (3.10+)
- Verify math library availability

**Runtime Crashes**
- Check input graph validity with `bliss_is_valid_graph()`
- Ensure proper memory management in callbacks
- Use debug build for detailed error information

**Performance Issues**
- Try different splitting heuristics
- Enable component recursion for disconnected graphs
- Consider graph preprocessing for very large inputs

### Debug Builds

Debug builds include:
- Memory leak detection
- Extensive validation checks
- Detailed algorithm tracing
- Performance profiling hooks

### Memory Debugging

```bash
# Build with address sanitizer
cmake -DCMAKE_BUILD_TYPE=Debug ..
make

# Run with memory checking
./test_bliss_comprehensive
```

## Contributing

### Development Setup

1. Fork the repository
2. Create feature branch
3. Make changes with tests
4. Ensure all tests pass
5. Submit pull request

### Code Style

- Follow C11 standards
- Use consistent naming conventions
- Add comprehensive documentation
- Include unit tests for new features
- Maintain backward compatibility

### Testing Requirements

- All new features must include tests
- Maintain >95% test coverage
- Performance regression tests for algorithm changes
- Memory leak detection in CI

## License

This project is licensed under the GNU Lesser General Public License v3.0 (LGPL-3.0).

The original bliss C++ library is:
- Copyright (c) 2003-2021 Tommi Junttila and Petteri Kaski
- Licensed under LGPL-3.0

This C translation maintains the same license terms.

## Citation

If you use this library in academic work, please cite:

```
Tommi Junttila and Petteri Kaski: Engineering an efficient canonical labeling 
tool for large and sparse graphs. In Proceedings of the Ninth Workshop on 
Algorithm Engineering and Experiments (ALENEX07), pages 135-149, 2007.
```

## Changelog

### Version 0.77c.0
- Initial C translation from C++ bliss
- Complete API compatibility
- Comprehensive test suite
- Performance optimizations
- Memory leak fixes
- CMake build system

## Support

### Documentation
- API reference: See header file `bliss.h`
- Examples: `examples/` directory
- Test cases: `bliss_tests.c`

### Issues
- Report bugs via issue tracker
- Include minimal reproduction case
- Provide system information and build details

### Performance Questions
- Include graph characteristics (size, density, symmetry)
- Try different heuristics before reporting
- Consider preprocessing for special graph types

## Acknowledgments

- Original bliss authors: Tommi Junttila and Petteri Kaski
- Graph automorphism community for algorithm development
- Contributors to the C translation effort: `amexcompliance`, `taboxdev`