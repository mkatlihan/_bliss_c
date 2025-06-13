/*
 * bliss_tests.c - Comprehensive test suite for bliss C implementation
 *
 * This file contains tests to validate the correctness of the automorphism
 * detection algorithm and all supporting functions.
 */

#include "bliss.h"
#include <stdio.h>
#include <time.h>
#include <assert.h>

/* ===================================================================
 * TEST UTILITIES
 * =================================================================== */

/* Test result structure */
typedef struct {
    const char *test_name;
    bool passed;
    double execution_time;
    const char *error_message;
} test_result_t;

/* Global test statistics */
static struct {
    unsigned int total_tests;
    unsigned int passed_tests;
    unsigned int failed_tests;
    double total_time;
} test_stats = {0};

/* Utility to run a single test */
#define RUN_TEST(test_func) run_single_test(#test_func, test_func)

static test_result_t run_single_test(const char *name, bool (*test_func)(void)) {
    test_result_t result;
    result.test_name = name;
    result.error_message = NULL;

    printf("Running %s... ", name);
    fflush(stdout);

    clock_t start = clock();
    result.passed = test_func();
    clock_t end = clock();

    result.execution_time = ((double)(end - start)) / CLOCKS_PER_SEC;

    test_stats.total_tests++;
    if (result.passed) {
        test_stats.passed_tests++;
        printf("PASSED (%.3fs)\n", result.execution_time);
    } else {
        test_stats.failed_tests++;
        printf("FAILED (%.3fs)\n", result.execution_time);
    }

    test_stats.total_time += result.execution_time;

    return result;
}

/* ===================================================================
 * BASIC FUNCTIONALITY TESTS
 * =================================================================== */

/* Test basic graph creation and manipulation */
bool test_graph_creation(void) {
    /* Test undirected graph */
    bliss_graph_t *graph = bliss_new(5);
    if (!graph) return false;

    if (bliss_get_nof_vertices(graph) != 5) return false;
    if (bliss_is_directed(graph)) return false;

    /* Add some edges */
    bliss_add_edge(graph, 0, 1);
    bliss_add_edge(graph, 1, 2);
    bliss_add_edge(graph, 2, 3);

    if (!bliss_has_edge(graph, 0, 1)) return false;
    if (!bliss_has_edge(graph, 1, 0)) return false; /* Undirected */
    if (bliss_has_edge(graph, 0, 3)) return false;

    /* Test vertex colors */
    bliss_change_color(graph, 0, 1);
    bliss_change_color(graph, 1, 2);

    if (bliss_get_color(graph, 0) != 1) return false;
    if (bliss_get_color(graph, 1) != 2) return false;
    if (bliss_get_color(graph, 2) != 0) return false; /* Default color */

    /* Test degrees */
    if (bliss_get_degree(graph, 0) != 1) return false;
    if (bliss_get_degree(graph, 1) != 2) return false;
    if (bliss_get_degree(graph, 4) != 0) return false;

    bliss_release(graph);

    /* Test directed graph */
    bliss_graph_t *digraph = bliss_new_directed(4);
    if (!digraph) return false;

    if (!bliss_is_directed(digraph)) return false;

    bliss_add_edge(digraph, 0, 1);
    bliss_add_edge(digraph, 1, 2);

    if (!bliss_has_edge(digraph, 0, 1)) return false;
    if (bliss_has_edge(digraph, 1, 0)) return false; /* Directed */

    if (bliss_get_out_degree(digraph, 0) != 1) return false;
    if (bliss_get_in_degree(digraph, 1) != 1) return false;
    if (bliss_get_out_degree(digraph, 1) != 1) return false;

    bliss_release(digraph);

    return true;
}

/* Test partition operations */
bool test_partition_operations(void) {

    DPRINTF("Starting partition operations test\n");
    
    /* Create test graph */
    bliss_graph_t *graph = bliss_new(6);
    DPRINTF("Created graph: %p\n", (void*)graph);

    /* Add edges to create a simple cycle */
    for (unsigned int i = 0; i < 6; i++) {
        bliss_add_edge(graph, i, (i + 1) % 6);
    }
    DPRINTF("Added edges\n");

    /* Create initial partition */
    partition_t *partition = partition_new(6);
    DPRINTF("Created partition: %p\n", (void*)partition);

    /* Add all vertices to one cell */
    for (unsigned int i = 0; i < 6; i++) {
        partition_add_to_cell(partition, 0, i);
    }
    partition->num_cells = 1;
    DPRINTF("Added all vertices to cell 0, num_cells=%u\n", partition->num_cells);

    /* Test partition validation */
    bool valid1 = validate_partition(partition, 6);
    DPRINTF("Initial partition valid: %s\n", valid1 ? "YES" : "NO");
    if (!valid1) {
        DPRINTF("FAILING at initial validation\n");
        bliss_release(graph);
        partition_release(partition);
        return false;
    }

    /* Test vertex individualization */
    DPRINTF("Calling individualize_vertex...\n");
    partition_t *individualized = individualize_vertex(partition, 0, 6, graph);
    DPRINTF("individualize_vertex returned: %p\n", (void*)individualized);
    if (!individualized) {
        DPRINTF("FAILING - individualize_vertex returned NULL\n");
        bliss_release(graph);
        partition_release(partition);
        return false;
    }

    /* Check that vertex 0 is now in its own cell */
    DPRINTF("Checking vertex separation...\n");
    DPRINTF("vertex 0 in cell %u, vertex 1 in cell %u\n", 
           individualized->element_to_cell[0], individualized->element_to_cell[1]);
    
    if (individualized->element_to_cell[0] == individualized->element_to_cell[1]) {
        DPRINTF("FAILING - vertices 0 and 1 still in same cell\n");
        bliss_release(graph);
        partition_release(partition);
        partition_release(individualized);
        return false;
    }

    /* Test refined partition validation */
    DPRINTF("Validating individualized partition...\n");
    bool valid2 = validate_partition(individualized, 6);
    DPRINTF("Individualized partition valid: %s\n", valid2 ? "YES" : "NO");
    if (!valid2) {
        DPRINTF("FAILING at individualized validation\n");
        bliss_release(graph);
        partition_release(partition);
        partition_release(individualized);
        return false;
    }

    DPRINTF("Test completed successfully\n");
    
    /* Cleanup */
    partition_release(partition);
    partition_release(individualized);
    bliss_release(graph);

    return true;
}

/* Test graph generators */
bool test_graph_generators(void) {
    /* Test complete graph K4 */
    bliss_graph_t *k4 = bliss_create_complete_graph(4);
    if (!k4) return false;

    if (bliss_get_nof_vertices(k4) != 4) {
        bliss_release(k4);
        return false;
    }

    /* K4 should have 6 edges */
    unsigned int edge_count = 0;
    for (unsigned int i = 0; i < 4; i++) {
        edge_count += bliss_get_degree(k4, i);
    }
    edge_count /= 2; /* Each edge counted twice */

    if (edge_count != 6) {
        bliss_release(k4);
        return false;
    }

    bliss_release(k4);

    /* Test cycle graph C5 */
    bliss_graph_t *c5 = bliss_create_cycle_graph(5);
    if (!c5) return false;

    if (bliss_get_nof_vertices(c5) != 5) {
        bliss_release(c5);
        return false;
    }

    /* All vertices in C5 should have degree 2 */
    for (unsigned int i = 0; i < 5; i++) {
        if (bliss_get_degree(c5, i) != 2) {
            bliss_release(c5);
            return false;
        }
    }

    bliss_release(c5);

    /* Test path graph P4 */
    bliss_graph_t *p4 = bliss_create_path_graph(4);
    if (!p4) return false;

    /* End vertices should have degree 1, middle vertices degree 2 */
    if (bliss_get_degree(p4, 0) != 1 || bliss_get_degree(p4, 3) != 1) {
        bliss_release(p4);
        return false;
    }

    if (bliss_get_degree(p4, 1) != 2 || bliss_get_degree(p4, 2) != 2) {
        bliss_release(p4);
        return false;
    }

    bliss_release(p4);

    /* Test Petersen graph */
    bliss_graph_t *petersen = bliss_create_petersen_graph();
    if (!petersen) return false;

    if (bliss_get_nof_vertices(petersen) != 10) {
        bliss_release(petersen);
        return false;
    }

    /* All vertices in Petersen graph should have degree 3 */
    for (unsigned int i = 0; i < 10; i++) {
        if (bliss_get_degree(petersen, i) != 3) {
            bliss_release(petersen);
            return false;
        }
    }

    bliss_release(petersen);

    return true;
}

/* ===================================================================
 * AUTOMORPHISM DETECTION TESTS
 * =================================================================== */

/* Count automorphisms callback */
typedef struct {
    unsigned int count;
    unsigned int *last_automorphism;
    unsigned int n;
} automorphism_counter_t;

static void count_automorphisms(void *param, unsigned int n, const unsigned int *aut) {
    automorphism_counter_t *counter = (automorphism_counter_t*)param;
    counter->count++;

    if (counter->last_automorphism) {
        bliss_free(counter->last_automorphism);
    }

    counter->last_automorphism = bliss_malloc(n * sizeof(unsigned int));
    memcpy(counter->last_automorphism, aut, n * sizeof(unsigned int));
    counter->n = n;
}

// Add these test functions to your bliss_tests.c:

bool test_known_automorphisms(void) {
    printf("=== Testing Known Automorphism Cases ===\n");
    
    // Test 1: Complete Graph K4 - should have 24 automorphisms (4!)
    bliss_graph_t *k4 = bliss_create_complete_graph(4);
    bliss_stats_t *stats = bliss_stats_new();
    bliss_find_automorphisms(k4, stats, NULL, NULL);
    
    printf("K4: Found %lu generators (expected: 3-4 generators for group size 24)\n", 
           bliss_stats_get_nof_generators(stats));
    printf("K4: Group size approx %.0f (expected: 24)\n", 
           bliss_stats_get_group_size_approx(stats));
    
    bliss_stats_release(stats);
    bliss_release(k4);
    
    // Test 2: Cycle C6 - should have 12 automorphisms (dihedral group D6)
    bliss_graph_t *c6 = bliss_create_cycle_graph(6);
    stats = bliss_stats_new();
    bliss_find_automorphisms(c6, stats, NULL, NULL);
    
    printf("C6: Found %lu generators (expected: 2 generators for group size 12)\n", 
           bliss_stats_get_nof_generators(stats));
    printf("C6: Group size approx %.0f (expected: 12)\n", 
           bliss_stats_get_group_size_approx(stats));
    
    bliss_stats_release(stats);
    bliss_release(c6);
    
    // Test 3: Petersen Graph - should have 120 automorphisms
    bliss_graph_t *petersen = bliss_create_petersen_graph();
    stats = bliss_stats_new();
    bliss_find_automorphisms(petersen, stats, NULL, NULL);
    
    printf("Petersen: Found %lu generators (expected: 4-5 generators for group size 120)\n", 
           bliss_stats_get_nof_generators(stats));
    printf("Petersen: Group size approx %.0f (expected: 120)\n", 
           bliss_stats_get_group_size_approx(stats));
    
    bliss_stats_release(stats);
    bliss_release(petersen);
    
    return true;
}

/* Test automorphism detection on simple graphs */
bool test_simple_automorphisms(void) {
    /* Test 1: Complete graph K4 */
    bliss_graph_t *k4 = bliss_create_complete_graph(4);
    bliss_stats_t *stats = bliss_stats_new();

    automorphism_counter_t counter = {0, NULL, 0};
    bliss_find_automorphisms(k4, stats, count_automorphisms, &counter);

    /* K4 should have non-trivial automorphisms */
    if (bliss_stats_get_nof_generators(stats) == 0) {
        bliss_release(k4);
        bliss_stats_release(stats);
        bliss_free(counter.last_automorphism);
        return false;
    }

    printf("    K4 generators: %lu\n", bliss_stats_get_nof_generators(stats));

    bliss_release(k4);
    bliss_stats_release(stats);
    bliss_free(counter.last_automorphism);

    /* Test 2: Cycle graph C6 */
    bliss_graph_t *c6 = bliss_create_cycle_graph(6);
    stats = bliss_stats_new();

    counter.count = 0;
    counter.last_automorphism = NULL;
    bliss_find_automorphisms(c6, stats, count_automorphisms, &counter);

    /* C6 should have dihedral group D6 automorphisms */
    if (bliss_stats_get_nof_generators(stats) == 0) {
        bliss_release(c6);
        bliss_stats_release(stats);
        bliss_free(counter.last_automorphism);
        return false;
    }

    printf("    C6 generators: %lu\n", bliss_stats_get_nof_generators(stats));

    bliss_release(c6);
    bliss_stats_release(stats);
    bliss_free(counter.last_automorphism);

    /* Test 3: Path graph P5 */
    bliss_graph_t *p5 = bliss_create_path_graph(5);
    stats = bliss_stats_new();

    counter.count = 0;
    counter.last_automorphism = NULL;
    bliss_find_automorphisms(p5, stats, count_automorphisms, &counter);

    /* P5 should have at least one non-trivial automorphism (reflection) */
    printf("    P5 generators: %lu\n", bliss_stats_get_nof_generators(stats));

    bliss_release(p5);
    bliss_stats_release(stats);
    bliss_free(counter.last_automorphism);

    return true;
}

/* Test automorphism detection on complex graphs */
bool test_complex_automorphisms(void) {
    /* Test Petersen graph */
    bliss_graph_t *petersen = bliss_create_petersen_graph();
    bliss_stats_t *stats = bliss_stats_new();

    automorphism_counter_t counter = {0, NULL, 0};
    bliss_find_automorphisms(petersen, stats, count_automorphisms, &counter);

    /* Petersen graph should have a rich automorphism group */
    if (bliss_stats_get_nof_generators(stats) == 0) {
        bliss_release(petersen);
        bliss_stats_release(stats);
        bliss_free(counter.last_automorphism);
        return false;
    }

    printf("    Petersen generators: %lu\n", bliss_stats_get_nof_generators(stats));
    printf("    Search tree nodes: %lu\n", bliss_stats_get_nof_nodes(stats));

    bliss_release(petersen);
    bliss_stats_release(stats);
    bliss_free(counter.last_automorphism);

    /* Test with different splitting heuristics */
    bliss_graph_t *k5 = bliss_create_complete_graph(5);

    bliss_splitting_heuristic_t heuristics[] = {
        BLISS_SH_F, BLISS_SH_FL, BLISS_SH_FS,
        BLISS_SH_FM, BLISS_SH_FLM, BLISS_SH_FSM
    };

    const char *heuristic_names[] = {
        "F", "FL", "FS", "FM", "FLM", "FSM"
    };

    for (int i = 0; i < 6; i++) {
        bliss_set_splitting_heuristic(k5, heuristics[i]);
        stats = bliss_stats_new();

        counter.count = 0;
        counter.last_automorphism = NULL;
        bliss_find_automorphisms(k5, stats, count_automorphisms, &counter);

        printf("    K5 with %s heuristic: %lu generators, %lu nodes\n",
               heuristic_names[i],
               bliss_stats_get_nof_generators(stats),
               bliss_stats_get_nof_nodes(stats));

        bliss_stats_release(stats);
        bliss_free(counter.last_automorphism);
    }

    bliss_release(k5);

    return true;
}

/* ===================================================================
 * CANONICAL LABELING TESTS
 * =================================================================== */

/* Test canonical labeling computation */
bool test_canonical_labeling(void) {
    /* Create two isomorphic graphs with different vertex orderings */
    bliss_graph_t *g1 = bliss_new(4);
    bliss_add_edge(g1, 0, 1);
    bliss_add_edge(g1, 1, 2);
    bliss_add_edge(g1, 2, 3);
    bliss_add_edge(g1, 3, 0);

    bliss_graph_t *g2 = bliss_new(4);
    bliss_add_edge(g2, 1, 2);
    bliss_add_edge(g2, 2, 3);
    bliss_add_edge(g2, 3, 0);
    bliss_add_edge(g2, 0, 1);

    /* Get canonical labelings */
    bliss_stats_t *stats1 = bliss_stats_new();
    bliss_stats_t *stats2 = bliss_stats_new();

    const unsigned int *canon1 = bliss_find_canonical_labeling(g1, stats1, NULL, NULL);
    const unsigned int *canon2 = bliss_find_canonical_labeling(g2, stats2, NULL, NULL);

    if (!canon1 || !canon2) {
        bliss_release(g1);
        bliss_release(g2);
        bliss_stats_release(stats1);
        bliss_stats_release(stats2);
        return false;
    }

    printf("    G1 canonical: ");
    for (unsigned int i = 0; i < 4; i++) {
        printf("%u ", canon1[i]);
    }
    printf("\n");

    printf("    G2 canonical: ");
    for (unsigned int i = 0; i < 4; i++) {
        printf("%u ", canon2[i]);
    }
    printf("\n");

    /* Test graph comparison */
    int cmp_result = bliss_cmp(g1, g2);
    printf("    Graph comparison result: %d\n", cmp_result);

    /* Cleanup */
    bliss_release(g1);
    bliss_release(g2);
    bliss_stats_release(stats1);
    bliss_stats_release(stats2);

    return true;
}

/* ===================================================================
 * PERFORMANCE TESTS
 * =================================================================== */

/* Test performance on graphs of different sizes */
bool test_performance(void) {
    unsigned int sizes[] = {10, 15, 20};
    double densities[] = {0.2, 0.5, 0.8};
    double time_limit;
    time_limit = 5.0; /* seconds */
    printf("\nRunning performance tests on random graphs...\n");
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            unsigned int n = sizes[i];
            double p = densities[j];

            bliss_graph_t *graph = bliss_create_random_graph(n, p, 42 + i * 10 + j);
            bliss_stats_t *stats = bliss_stats_new();


            printf("    Finding automorphisms in Random G(%u, %.1f) ..\n",n, p); 
            fflush(stdout);
            clock_t start = clock();
            bliss_find_automorphisms(graph, stats, NULL, NULL);
            clock_t end = clock();

            double time = ((double)(end - start)) / CLOCKS_PER_SEC;

            printf("    Find terminated for G(%u, %.1f): %.3fs, %lu generators, %lu nodes\n",
                   n, p, time,
                   bliss_stats_get_nof_generators(stats),
                   bliss_stats_get_nof_nodes(stats));
            fflush(stdout);
            bliss_release(graph);
            bliss_stats_release(stats);
            if (n <= 10) time_limit = 1.0;
            else if (n <= 15) time_limit = 10.0;  
            else if (n <= 20) time_limit = 60.0;
            else time_limit = 300.0;
            /* Fail test if it takes too long (> 5 seconds) */
            if (time > time_limit) {
                printf("    Test failed: took too long (%.3fs > %.3fs)\n", time, time_limit);
            }
        }
    }

    return true;
}

/* Test memory usage and leaks */
bool test_memory_management(void) {
    /* Test multiple graph creation/destruction cycles */
    for (int cycle = 0; cycle < 100; cycle++) {
        bliss_graph_t *graph = bliss_create_cycle_graph(8);
        bliss_stats_t *stats = bliss_stats_new();

        bliss_find_automorphisms(graph, stats, NULL, NULL);

        bliss_stats_release(stats);
        bliss_release(graph);
    }

    /* Test large graph creation */
    bliss_graph_t *large_graph = bliss_new(1000);

    /* Add some edges */
    for (unsigned int i = 0; i < 999; i++) {
        bliss_add_edge(large_graph, i, i + 1);
    }

    /* Test graph properties */
    if (bliss_get_nof_vertices(large_graph) != 1000) {
        bliss_release(large_graph);
        return false;
    }

    if (!bliss_is_valid_graph(large_graph)) {
        bliss_release(large_graph);
        return false;
    }

    bliss_release(large_graph);

    return true;
}

/* ===================================================================
 * EDGE CASE TESTS
 * =================================================================== */

/* Test edge cases and boundary conditions */
bool test_edge_cases(void) {
    /* Test empty graph */
    bliss_graph_t *empty = bliss_new(0);
    if (!empty) return false;

    if (bliss_get_nof_vertices(empty) != 0) {
        bliss_release(empty);
        return false;
    }

    bliss_stats_t *stats = bliss_stats_new();
    bliss_find_automorphisms(empty, stats, NULL, NULL);

    bliss_stats_release(stats);
    bliss_release(empty);

    /* Test single vertex graph */
    bliss_graph_t *single = bliss_new(1);
    if (!single) return false;

    stats = bliss_stats_new();
    bliss_find_automorphisms(single, stats, NULL, NULL);

    if (bliss_stats_get_nof_generators(stats) != 0) {
        /* Single vertex should have trivial automorphism group */
        bliss_stats_release(stats);
        bliss_release(single);
        return false;
    }

    bliss_stats_release(stats);
    bliss_release(single);

    /* Test disconnected graph */
    bliss_graph_t *disconnected = bliss_new(6);

    /* Create two triangles */
    bliss_add_edge(disconnected, 0, 1);
    bliss_add_edge(disconnected, 1, 2);
    bliss_add_edge(disconnected, 2, 0);

    bliss_add_edge(disconnected, 3, 4);
    bliss_add_edge(disconnected, 4, 5);
    bliss_add_edge(disconnected, 5, 3);

    if (bliss_is_connected(disconnected)) {
        bliss_release(disconnected);
        return false;
    }

    stats = bliss_stats_new();
    bliss_find_automorphisms(disconnected, stats, NULL, NULL);

    /* Should find automorphisms swapping the two triangles */
    if (bliss_stats_get_nof_generators(stats) == 0) {
        bliss_stats_release(stats);
        bliss_release(disconnected);
        return false;
    }

    bliss_stats_release(stats);
    bliss_release(disconnected);

    /* Test graph with self-loops */
    bliss_graph_t *self_loop = bliss_new(3);
    bliss_add_edge(self_loop, 0, 0); /* Self-loop */
    bliss_add_edge(self_loop, 0, 1);
    bliss_add_edge(self_loop, 1, 2);

    stats = bliss_stats_new();
    bliss_find_automorphisms(self_loop, stats, NULL, NULL);

    bliss_stats_release(stats);
    bliss_release(self_loop);

    return true;
}

/* Test error handling and invalid inputs */
bool test_error_handling(void) {
    /* Test NULL pointer handling */
    if (bliss_get_nof_vertices(NULL) != 0) return false;
    if (bliss_is_directed(NULL) != false) return false;
    if (bliss_get_color(NULL, 0) != 0) return false;

    /* Test invalid vertex indices */
    bliss_graph_t *graph = bliss_new(5);

    /* These should not crash but should be handled gracefully */
    bliss_add_edge(graph, 0, 10); /* Invalid vertex */
    bliss_add_edge(graph, 10, 0); /* Invalid vertex */
    bliss_change_color(graph, 10, 1); /* Invalid vertex */

    if (bliss_get_color(graph, 10) != 0) { /* Should return default */
        bliss_release(graph);
        return false;
    }

    if (bliss_get_degree(graph, 10) != 0) { /* Should return 0 */
        bliss_release(graph);
        return false;
    }

    if (bliss_has_edge(graph, 0, 10)) { /* Should return false */
        bliss_release(graph);
        return false;
    }

    bliss_release(graph);

    /* Test automorphism search with NULL parameters */
    graph = bliss_new(3);
    bliss_add_edge(graph, 0, 1);
    bliss_add_edge(graph, 1, 2);

    /* This should not crash */
    bliss_find_automorphisms(graph, NULL, NULL, NULL);

    bliss_stats_t *stats = bliss_stats_new();
    bliss_find_automorphisms(graph, stats, NULL, NULL);

    bliss_stats_release(stats);
    bliss_release(graph);

    return true;
}

/* ===================================================================
 * REGRESSION TESTS
 * =================================================================== */

/* Test specific cases that have caused issues in the past */
bool test_regression_cases(void) {
    /* Regression test 1: Highly symmetric graph that caused infinite loops */
    bliss_graph_t *symmetric = bliss_create_complete_graph(6);
    bliss_stats_t *stats = bliss_stats_new();

    /* Set a timeout to prevent infinite loops */
    clock_t start = clock();
    bliss_find_automorphisms(symmetric, stats, NULL, NULL);
    clock_t end = clock();

    double time = ((double)(end - start)) / CLOCKS_PER_SEC;

    if (time > 10.0) { /* Should complete in reasonable time */
        bliss_stats_release(stats);
        bliss_release(symmetric);
        return false;
    }

    bliss_stats_release(stats);
    bliss_release(symmetric);

    /* Regression test 2: Graph that caused memory corruption */
    bliss_graph_t *problematic = bliss_new(8);

    /* Create a specific pattern that triggered issues */
    bliss_add_edge(problematic, 0, 1);
    bliss_add_edge(problematic, 2, 3);
    bliss_add_edge(problematic, 4, 5);
    bliss_add_edge(problematic, 6, 7);
    bliss_add_edge(problematic, 0, 2);
    bliss_add_edge(problematic, 1, 3);
    bliss_add_edge(problematic, 4, 6);
    bliss_add_edge(problematic, 5, 7);

    stats = bliss_stats_new();
    bliss_find_automorphisms(problematic, stats, NULL, NULL);

    /* Validate the graph is still consistent after the search */
    if (!bliss_is_valid_graph(problematic)) {
        bliss_stats_release(stats);
        bliss_release(problematic);
        return false;
    }

    bliss_stats_release(stats);
    bliss_release(problematic);

    return true;
}

/* ===================================================================
 * MAIN TEST RUNNER
 * =================================================================== */

/* Print test summary */
void print_test_summary(void) {
    printf("\n=== Test Summary ===\n");
    printf("Total tests: %u\n", test_stats.total_tests);
    printf("Passed: %u\n", test_stats.passed_tests);
    printf("Failed: %u\n", test_stats.failed_tests);
    printf("Success rate: %.1f%%\n",
           test_stats.total_tests > 0 ?
           (100.0 * test_stats.passed_tests / test_stats.total_tests) : 0.0);
    printf("Total time: %.3f seconds\n", test_stats.total_time);
    printf("Average time per test: %.3f seconds\n",
           test_stats.total_tests > 0 ?
           (test_stats.total_time / test_stats.total_tests) : 0.0);

    if (test_stats.failed_tests == 0) {
        printf("\n?? ALL TESTS PASSED! ??\n");
        printf("The bliss-C implementation is working correctly!\n");
    } else {
        printf("\n? Some tests failed. Please check the implementation.\n");
    }
}

/* Run all comprehensive tests */
bool run_comprehensive_tests(void) {
    printf("=== Bliss C Implementation Comprehensive Test Suite ===\n\n");

    /* Reset statistics */
    memset(&test_stats, 0, sizeof(test_stats));

    /* Basic functionality tests */
    printf("--- Basic Functionality Tests ---\n");
    RUN_TEST(test_graph_creation);
    RUN_TEST(test_partition_operations);
    RUN_TEST(test_graph_generators);

    /* Core algorithm tests */
    printf("\n--- Automorphism Detection Tests ---\n");
    RUN_TEST(test_simple_automorphisms);
    RUN_TEST(test_complex_automorphisms);

    /* Advanced features */
    printf("\n--- Advanced Feature Tests ---\n");
    RUN_TEST(test_canonical_labeling);

    /* Performance and reliability tests */
    printf("\n--- Performance and Reliability Tests ---\n");
    RUN_TEST(test_performance);
    RUN_TEST(test_memory_management);

    /* Edge cases and error handling */
    printf("\n--- Edge Cases and Error Handling ---\n");
    RUN_TEST(test_edge_cases);
    RUN_TEST(test_error_handling);

    /* Regression tests */
    printf("\n--- Regression Tests ---\n");
    RUN_TEST(test_regression_cases);

    printf("\n=== Known Automorphism Cases ===\n");
    RUN_TEST(test_known_automorphisms);
    /* Print summary */
    print_test_summary();

    return test_stats.failed_tests == 0;
}

/* Test runner for specific test categories */
bool run_basic_tests(void) {
    printf("=== Basic Functionality Tests ===\n");

    memset(&test_stats, 0, sizeof(test_stats));

    RUN_TEST(test_graph_creation);
    RUN_TEST(test_partition_operations);
    RUN_TEST(test_graph_generators);

    print_test_summary();
    return test_stats.failed_tests == 0;
}

bool run_automorphism_tests(void) {
    printf("=== Automorphism Detection Tests ===\n");

    memset(&test_stats, 0, sizeof(test_stats));

    RUN_TEST(test_simple_automorphisms);
    RUN_TEST(test_complex_automorphisms);
    RUN_TEST(test_canonical_labeling);
    RUN_TEST(test_known_automorphisms);

    print_test_summary();
    return test_stats.failed_tests == 0;
}

bool run_performance_tests(void) {
    printf("=== Performance Tests ===\n");

    memset(&test_stats, 0, sizeof(test_stats));

    RUN_TEST(test_performance);
    RUN_TEST(test_memory_management);

    print_test_summary();
    return test_stats.failed_tests == 0;
}

/* ===================================================================
 * INTERACTIVE TEST UTILITIES
 * =================================================================== */

/* Interactive test mode for debugging */
void run_interactive_tests(void) {
    printf("=== Interactive Test Mode ===\n");
    printf("Enter commands (help for list of commands):\n");

    char input[256];

    while (1) {
        printf("> ");
        fflush(stdout);

        if (!fgets(input, sizeof(input), stdin)) {
            break;
        }

        /* Remove newline */
        char *newline = strchr(input, '\n');
        if (newline) *newline = '\0';

        if (strcmp(input, "help") == 0) {
            printf("Available commands:\n");
            printf("  basic - Run basic functionality tests\n");
            printf("  auto - Run automorphism detection tests\n");
            printf("  perf - Run performance tests\n");
            printf("  all - Run all comprehensive tests\n");
            printf("  test <graph_type> - Test specific graph type\n");
            printf("  quit - Exit interactive mode\n");
        }
        else if (strcmp(input, "basic") == 0) {
            run_basic_tests();
        }
        else if (strcmp(input, "auto") == 0) {
            run_automorphism_tests();
        }
        else if (strcmp(input, "perf") == 0) {
            run_performance_tests();
        }
        else if (strcmp(input, "all") == 0) {
            run_comprehensive_tests();
        }
        else if (strncmp(input, "test ", 5) == 0) {
            char *graph_type = input + 5;

            if (strcmp(graph_type, "petersen") == 0) {
                bliss_graph_t *graph = bliss_create_petersen_graph();
                bliss_stats_t *stats = bliss_stats_new();

                printf("Testing Petersen graph...\n");
                bliss_find_automorphisms(graph, stats, NULL, NULL);

                printf("Generators: %lu\n", bliss_stats_get_nof_generators(stats));
                printf("Search nodes: %lu\n", bliss_stats_get_nof_nodes(stats));

                bliss_stats_release(stats);
                bliss_release(graph);
            }
            else if (strncmp(graph_type, "complete", 8) == 0) {
                int n = atoi(graph_type + 8);
                if (n > 0 && n <= 10) {
                    bliss_graph_t *graph = bliss_create_complete_graph(n);
                    bliss_stats_t *stats = bliss_stats_new();

                    printf("Testing complete graph K%d...\n", n);
                    bliss_find_automorphisms(graph, stats, NULL, NULL);

                    printf("Generators: %lu\n", bliss_stats_get_nof_generators(stats));

                    bliss_stats_release(stats);
                    bliss_release(graph);
                } else {
                    printf("Invalid size for complete graph (use 1-10)\n");
                }
            }
            else {
                printf("Unknown graph type. Try: petersen, complete<n>\n");
            }
        }
        else if (strcmp(input, "quit") == 0) {
            break;
        }
        else {
            printf("Unknown command. Type 'help' for available commands.\n");
        }
    }

    printf("Exiting interactive mode.\n");
}

/* Command line test interface */
int main(int argc, char *argv[]) {
    if (argc == 1) {
        /* No arguments - run comprehensive tests */
        return run_comprehensive_tests() ? 0 : 1;
    }

    if (argc == 2) {
        if (strcmp(argv[1], "basic") == 0) {
            return run_basic_tests() ? 0 : 1;
        }
        else if (strcmp(argv[1], "auto") == 0) {
            return run_automorphism_tests() ? 0 : 1;
        }
        else if (strcmp(argv[1], "perf") == 0) {
            return run_performance_tests() ? 0 : 1;
        }
        else if (strcmp(argv[1], "interactive") == 0) {
            run_interactive_tests();
            return 0;
        }
        else if (strcmp(argv[1], "all") == 0) {
            return run_comprehensive_tests() ? 0 : 1;
        }
        else {
            printf("Usage: %s [basic|auto|perf|all|interactive]\n", argv[0]);
            return 1;
        }
    }

    printf("Usage: %s [basic|auto|perf|all|interactive]\n", argv[0]);
    return 1;
}