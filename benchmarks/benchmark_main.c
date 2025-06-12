/*
 * benchmark_main.c - Performance benchmarks for bliss C library
 * 
 * Clean implementation with working benchmarks for graph automorphisms.
 */

#include "../bliss.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

/* Global variables for automorphism counting */
static unsigned int automorphism_count = 0;

/* Automorphism counting callback */
static void count_automorphisms(void *user_param, unsigned int n, const unsigned int *aut) {
    (void)user_param;
    (void)n;
    (void)aut;
    automorphism_count++;
}

/* Benchmark a single graph */
static double benchmark_graph_simple(bliss_graph_t *graph, const char *description) {
    printf("  %s: ", description);
    fflush(stdout);
    
    automorphism_count = 0;
    
    clock_t start = clock();
    bliss_stats_t *stats = bliss_stats_new();
    bliss_find_automorphisms(graph, stats, count_automorphisms, NULL);
    clock_t end = clock();
    
    double time_taken = ((double)(end - start)) / CLOCKS_PER_SEC;
    unsigned long nodes = bliss_stats_get_nof_nodes(stats);
    
    printf("%.3fs, %lu nodes, %u automorphisms\n", time_taken, nodes, automorphism_count);
    
    bliss_stats_release(stats);
    return time_taken;
}

/* Benchmark complete graphs */
static void benchmark_complete_graphs(void) {
    printf("\nBenchmarking complete graphs:\n");
    
    for (unsigned int n = 4; n <= 8; n++) {
        bliss_graph_t *graph = bliss_create_complete_graph(n);
        if (!graph) continue;
        
        char desc[64];
        snprintf(desc, sizeof(desc), "K_%u", n);
        benchmark_graph_simple(graph, desc);
        
        bliss_release(graph);
    }
}

/* Benchmark cycle graphs */
static void benchmark_cycle_graphs(void) {
    printf("\nBenchmarking cycle graphs:\n");
    
    unsigned int sizes[] = {10, 20, 50, 100};
    int num_sizes = sizeof(sizes) / sizeof(sizes[0]);
    
    for (int i = 0; i < num_sizes; i++) {
        bliss_graph_t *graph = bliss_create_cycle_graph(sizes[i]);
        if (!graph) continue;
        
        char desc[64];
        snprintf(desc, sizeof(desc), "C_%u", sizes[i]);
        benchmark_graph_simple(graph, desc);
        
        bliss_release(graph);
    }
}

/* Benchmark path graphs */
static void benchmark_path_graphs(void) {
    printf("\nBenchmarking path graphs:\n");
    
    for (unsigned int n = 5; n <= 20; n += 5) {
        bliss_graph_t *graph = bliss_create_path_graph(n);
        if (!graph) continue;
        
        char desc[64];
        snprintf(desc, sizeof(desc), "P_%u", n);
        benchmark_graph_simple(graph, desc);
        
        bliss_release(graph);
    }
}

/* Benchmark star graphs */
static void benchmark_star_graphs(void) {
    printf("\nBenchmarking star graphs:\n");
    
    for (unsigned int n = 5; n <= 15; n += 5) {
        bliss_graph_t *graph = bliss_create_star_graph(n);
        if (!graph) continue;
        
        char desc[64];
        snprintf(desc, sizeof(desc), "S_%u", n);
        benchmark_graph_simple(graph, desc);
        
        bliss_release(graph);
    }
}

/* Benchmark Petersen graph */
static void benchmark_petersen(void) {
    printf("\nBenchmarking Petersen graph:\n");
    
    bliss_graph_t *graph = bliss_create_petersen_graph();
    if (!graph) return;
    
    benchmark_graph_simple(graph, "Petersen");
    bliss_release(graph);
}

/* Benchmark canonical labeling */
static void benchmark_canonical_labeling(void) {
    printf("\nBenchmarking canonical labeling:\n");
    
    unsigned int sizes[] = {8, 12, 16, 20};
    int num_sizes = sizeof(sizes) / sizeof(sizes[0]);
    
    for (int i = 0; i < num_sizes; i++) {
        bliss_graph_t *graph = bliss_create_cycle_graph(sizes[i]);
        if (!graph) continue;
        
        clock_t start = clock();
        bliss_stats_t *stats = bliss_stats_new();
        const unsigned int *labeling = bliss_find_canonical_labeling(graph, stats, NULL, NULL);
        clock_t end = clock();
        
        double time_taken = ((double)(end - start)) / CLOCKS_PER_SEC;
        
        printf("  C_%u canonical: %.3fs, labeling %s\n", 
               sizes[i], time_taken, labeling ? "found" : "failed");
        
        bliss_stats_release(stats);
        bliss_release(graph);
    }
}

/* Benchmark splitting heuristics */
static void benchmark_splitting_heuristics(void) {
    printf("\nBenchmarking splitting heuristics on cycle C_20:\n");
    
    bliss_splitting_heuristic_t heuristics[] = {
        BLISS_SH_F, BLISS_SH_FL, BLISS_SH_FS, 
        BLISS_SH_FM, BLISS_SH_FLM, BLISS_SH_FSM
    };
    const char *names[] = {
        "F", "FL", "FS", "FM", "FLM", "FSM"
    };
    
    int num_heuristics = sizeof(heuristics) / sizeof(heuristics[0]);
    
    for (int i = 0; i < num_heuristics; i++) {
        bliss_graph_t *graph = bliss_create_cycle_graph(20);
        if (!graph) continue;
        
        bliss_set_splitting_heuristic(graph, heuristics[i]);
        
        char desc[64];
        snprintf(desc, sizeof(desc), "C_20_%s", names[i]);
        benchmark_graph_simple(graph, desc);
        
        bliss_release(graph);
    }
}

/* Benchmark memory usage */
static void benchmark_memory_usage(void) {
    printf("\nBenchmarking memory usage:\n");
    
    unsigned int sizes[] = {100, 500, 1000, 5000};
    int num_sizes = sizeof(sizes) / sizeof(sizes[0]);
    
    for (int i = 0; i < num_sizes; i++) {
        bliss_graph_t *graph = bliss_new(sizes[i]);
        
        /* Create a cycle graph */
        for (unsigned int j = 0; j < sizes[i]; j++) {
            bliss_add_edge(graph, j, (j + 1) % sizes[i]);
        }
        
        size_t memory = bliss_estimate_memory_usage(graph);
        printf("  %u vertices: %zu bytes (%.2f KB)\n", 
               sizes[i], memory, memory / 1024.0);
        
        bliss_release(graph);
    }
}

/* Benchmark graph operations */
static void benchmark_graph_operations(void) {
    printf("\nBenchmarking graph operations:\n");
    
    const unsigned int n = 1000;
    const unsigned int num_ops = 10000;
    
    /* Edge addition benchmark */
    bliss_graph_t *graph = bliss_new(n);
    
    clock_t start = clock();
    for (unsigned int i = 0; i < num_ops; i++) {
        unsigned int v1 = rand() % n;
        unsigned int v2 = rand() % n;
        bliss_add_edge(graph, v1, v2);
    }
    clock_t end = clock();
    
    double edge_time = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("  Edge addition: %.3fs for %u operations (%.0f ops/sec)\n", 
           edge_time, num_ops, num_ops / edge_time);
    
    /* Hash computation benchmark */
    start = clock();
    for (unsigned int i = 0; i < 1000; i++) {
        bliss_get_hash(graph);
    }
    end = clock();
    
    double hash_time = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("  Hash computation: %.3fs for 1000 operations (%.0f ops/sec)\n", 
           hash_time, 1000.0 / hash_time);
    
    bliss_release(graph);
}

/* Print system information */
static void print_system_info(void) {
    printf("System Information:\n");
    printf("  Compiler: ");
#ifdef __GNUC__
    printf("GCC %d.%d.%d", __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__);
#elif defined(__clang__)
    printf("Clang");
#elif defined(_MSC_VER)
    printf("MSVC %d", _MSC_VER);
#else
    printf("Unknown");
#endif
    printf("\n");
    
    printf("  Build type: ");
#ifdef NDEBUG
    printf("Release");
#else
    printf("Debug");
#endif
    printf("\n");
    
    printf("  Bliss version: %s\n", bliss_version_string());
    printf("\n");
}

/* Main benchmark runner */
int main(int argc, char *argv[]) {
    printf("Bliss C Library Performance Benchmarks\n");
    printf("=======================================\n\n");
    
    /* Initialize random seed */
    srand(12345);
    
    print_system_info();
    
    /* Parse command line arguments */
    bool run_all = (argc == 1);
    bool run_basic = run_all || (argc > 1 && strstr(argv[1], "basic"));
    bool run_heuristics = run_all || (argc > 1 && strstr(argv[1], "heuristics"));
    bool run_memory = run_all || (argc > 1 && strstr(argv[1], "memory"));
    bool run_canonical = run_all || (argc > 1 && strstr(argv[1], "canonical"));
    
    if (argc > 1 && strcmp(argv[1], "help") == 0) {
        printf("Usage: %s [test_type]\n", argv[0]);
        printf("Test types:\n");
        printf("  basic      - Basic graph benchmarks\n");
        printf("  heuristics - Splitting heuristic tests\n");
        printf("  memory     - Memory usage tests\n");
        printf("  canonical  - Canonical labeling tests\n");
        printf("  all        - All benchmarks (default)\n\n");
        return 0;
    }
    
    clock_t total_start = clock();
    
    /* Run selected benchmarks */
    if (run_basic) {
        benchmark_complete_graphs();
        benchmark_cycle_graphs();
        benchmark_path_graphs();
        benchmark_star_graphs();
        benchmark_petersen();
        benchmark_graph_operations();
    }
    
    if (run_heuristics) {
        benchmark_splitting_heuristics();
    }
    
    if (run_memory) {
        benchmark_memory_usage();
    }
    
    if (run_canonical) {
        benchmark_canonical_labeling();
    }
    
    clock_t total_end = clock();
    double total_time = ((double)(total_end - total_start)) / CLOCKS_PER_SEC;
    
    printf("\nBenchmark Summary:\n");
    printf("==================\n");
    printf("  Total execution time: %.3fs\n", total_time);
    printf("  All benchmarks completed successfully!\n");
    
    return 0;
}