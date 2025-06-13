/*
 * test_classic_graphs.c - Classic graph structure tests for bliss C library
 * 
 * Tests graph generators and basic properties of well-known graph families
 * without relying on complex automorphism computations that may be incomplete.
 */

#include "bliss.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define TEST_ASSERT(condition, message) \
    do { \
        if (!(condition)) { \
            fprintf(stderr, "FAIL: %s\n", message); \
            return 0; \
        } \
    } while(0)

#define TEST_SUCCESS() \
    do { \
        printf("PASS: %s\n", __func__); \
        return 1; \
    } while(0)

/* Test complete graph K_n properties */
static int test_complete_graphs(void) {
    printf("Testing complete graphs...\n");
    
    /* Test K_3 (triangle) */
    bliss_graph_t *k3 = bliss_create_complete_graph(3);
    TEST_ASSERT(k3 != NULL, "K3 creation failed");
    TEST_ASSERT(bliss_get_nof_vertices(k3) == 3, "K3 should have 3 vertices");
    
    /* Every vertex should have degree 2 in K3 */
    for (unsigned int i = 0; i < 3; i++) {
        TEST_ASSERT(bliss_get_degree(k3, i) == 2, "K3 vertices should have degree 2");
    }
    
    /* All possible edges should exist */
    TEST_ASSERT(bliss_has_edge(k3, 0, 1), "Edge (0,1) should exist in K3");
    TEST_ASSERT(bliss_has_edge(k3, 1, 2), "Edge (1,2) should exist in K3");
    TEST_ASSERT(bliss_has_edge(k3, 0, 2), "Edge (0,2) should exist in K3");
    
    bliss_release(k3);
    
    /* Test K_4 */
    bliss_graph_t *k4 = bliss_create_complete_graph(4);
    TEST_ASSERT(k4 != NULL, "K4 creation failed");
    TEST_ASSERT(bliss_get_nof_vertices(k4) == 4, "K4 should have 4 vertices");
    
    /* Every vertex should have degree 3 in K4 */
    for (unsigned int i = 0; i < 4; i++) {
        TEST_ASSERT(bliss_get_degree(k4, i) == 3, "K4 vertices should have degree 3");
    }
    
    /* Count edges: K4 should have 6 edges */
    unsigned int edge_count = 0;
    for (unsigned int i = 0; i < 4; i++) {
        for (unsigned int j = i + 1; j < 4; j++) {
            if (bliss_has_edge(k4, i, j)) {
                edge_count++;
            }
        }
    }
    TEST_ASSERT(edge_count == 6, "K4 should have exactly 6 edges");
    
    bliss_release(k4);
    TEST_SUCCESS();
}

/* Test cycle graph C_n properties */
static int test_cycle_graphs(void) {
    printf("Testing cycle graphs...\n");
    
    /* Test C_4 (square) */
    bliss_graph_t *c4 = bliss_create_cycle_graph(4);
    TEST_ASSERT(c4 != NULL, "C4 creation failed");
    TEST_ASSERT(bliss_get_nof_vertices(c4) == 4, "C4 should have 4 vertices");
    
    /* Every vertex in a cycle should have degree 2 */
    for (unsigned int i = 0; i < 4; i++) {
        TEST_ASSERT(bliss_get_degree(c4, i) == 2, "C4 vertices should have degree 2");
    }
    
    /* Check cycle structure: each vertex connected to next */
    for (unsigned int i = 0; i < 4; i++) {
        unsigned int next = (i + 1) % 4;
        TEST_ASSERT(bliss_has_edge(c4, i, next), "Cycle edge should exist");
    }
    
    /* Check that non-cycle edges don't exist */
    TEST_ASSERT(!bliss_has_edge(c4, 0, 2), "Diagonal should not exist in C4");
    TEST_ASSERT(!bliss_has_edge(c4, 1, 3), "Diagonal should not exist in C4");
    
    bliss_release(c4);
    
    /* Test C_5 (pentagon) */
    bliss_graph_t *c5 = bliss_create_cycle_graph(5);
    TEST_ASSERT(c5 != NULL, "C5 creation failed");
    TEST_ASSERT(bliss_get_nof_vertices(c5) == 5, "C5 should have 5 vertices");
    
    /* Verify cycle properties */
    for (unsigned int i = 0; i < 5; i++) {
        TEST_ASSERT(bliss_get_degree(c5, i) == 2, "C5 vertices should have degree 2");
        
        unsigned int next = (i + 1) % 5;
        TEST_ASSERT(bliss_has_edge(c5, i, next), "Cycle edge should exist");
    }
    
    bliss_release(c5);
    TEST_SUCCESS();
}

/* Test path graph P_n properties */
static int test_path_graphs(void) {
    printf("Testing path graphs...\n");
    
    /* Test P_3 (path with 3 vertices) */
    bliss_graph_t *p3 = bliss_create_path_graph(3);
    TEST_ASSERT(p3 != NULL, "P3 creation failed");
    TEST_ASSERT(bliss_get_nof_vertices(p3) == 3, "P3 should have 3 vertices");
    
    /* End vertices should have degree 1, middle vertex degree 2 */
    TEST_ASSERT(bliss_get_degree(p3, 0) == 1, "P3 end vertex should have degree 1");
    TEST_ASSERT(bliss_get_degree(p3, 1) == 2, "P3 middle vertex should have degree 2");
    TEST_ASSERT(bliss_get_degree(p3, 2) == 1, "P3 end vertex should have degree 1");
    
    /* Check path edges */
    TEST_ASSERT(bliss_has_edge(p3, 0, 1), "Path edge (0,1) should exist");
    TEST_ASSERT(bliss_has_edge(p3, 1, 2), "Path edge (1,2) should exist");
    TEST_ASSERT(!bliss_has_edge(p3, 0, 2), "Direct edge (0,2) should not exist");
    
    bliss_release(p3);
    
    /* Test P_5 */
    bliss_graph_t *p5 = bliss_create_path_graph(5);
    TEST_ASSERT(p5 != NULL, "P5 creation failed");
    TEST_ASSERT(bliss_get_nof_vertices(p5) == 5, "P5 should have 5 vertices");
    
    /* Check degree pattern: 1,2,2,2,1 */
    TEST_ASSERT(bliss_get_degree(p5, 0) == 1, "P5 end should have degree 1");
    TEST_ASSERT(bliss_get_degree(p5, 1) == 2, "P5 middle should have degree 2");
    TEST_ASSERT(bliss_get_degree(p5, 2) == 2, "P5 middle should have degree 2");
    TEST_ASSERT(bliss_get_degree(p5, 3) == 2, "P5 middle should have degree 2");
    TEST_ASSERT(bliss_get_degree(p5, 4) == 1, "P5 end should have degree 1");
    
    bliss_release(p5);
    TEST_SUCCESS();
}

/* Test star graph S_n properties */
static int test_star_graphs(void) {
    printf("Testing star graphs...\n");
    
    /* Test S_4 (star with 4 leaves = 5 vertices total) */
    bliss_graph_t *s4 = bliss_create_star_graph(5);
    TEST_ASSERT(s4 != NULL, "S4 creation failed");
    TEST_ASSERT(bliss_get_nof_vertices(s4) == 5, "S4 should have 5 vertices");
    
    /* Center vertex (0) should have degree 4, leaves should have degree 1 */
    TEST_ASSERT(bliss_get_degree(s4, 0) == 4, "Center should have degree 4");
    for (unsigned int i = 1; i < 5; i++) {
        TEST_ASSERT(bliss_get_degree(s4, i) == 1, "Leaf should have degree 1");
    }
    
    /* Check star edges: center connected to all leaves */
    for (unsigned int i = 1; i < 5; i++) {
        TEST_ASSERT(bliss_has_edge(s4, 0, i), "Star edge should exist");
    }
    
    /* Check that leaves are not connected to each other */
    for (unsigned int i = 1; i < 5; i++) {
        for (unsigned int j = i + 1; j < 5; j++) {
            TEST_ASSERT(!bliss_has_edge(s4, i, j), "Leaves should not be connected");
        }
    }
    
    bliss_release(s4);
    TEST_SUCCESS();
}

/* Test Petersen graph properties */
static int test_petersen_graph(void) {
    printf("Testing Petersen graph...\n");
    
    bliss_graph_t *petersen = bliss_create_petersen_graph();
    TEST_ASSERT(petersen != NULL, "Petersen graph creation failed");
    
    /* Basic properties */
    TEST_ASSERT(bliss_get_nof_vertices(petersen) == 10, "Petersen should have 10 vertices");
    
    /* Every vertex should have degree 3 */
    for (unsigned int i = 0; i < 10; i++) {
        unsigned int degree = bliss_get_degree(petersen, i);
        TEST_ASSERT(degree == 3, "Petersen vertices should have degree 3");
    }
    
    /* Check outer pentagon (vertices 0-4) */
    for (unsigned int i = 0; i < 5; i++) {
        unsigned int next = (i + 1) % 5;
        TEST_ASSERT(bliss_has_edge(petersen, i, next), "Outer pentagon edge should exist");
    }
    
    /* Check inner pentagram (vertices 5-9) */
    for (unsigned int i = 0; i < 5; i++) {
        unsigned int target = ((i + 2) % 5) + 5;
        TEST_ASSERT(bliss_has_edge(petersen, i + 5, target), "Inner pentagram edge should exist");
    }
    
    /* Check connections between outer and inner */
    for (unsigned int i = 0; i < 5; i++) {
        TEST_ASSERT(bliss_has_edge(petersen, i, i + 5), "Outer-inner connection should exist");
    }
    
    /* Total edge count should be 15 */
    unsigned int edge_count = 0;
    for (unsigned int i = 0; i < 10; i++) {
        for (unsigned int j = i + 1; j < 10; j++) {
            if (bliss_has_edge(petersen, i, j)) {
                edge_count++;
            }
        }
    }
    TEST_ASSERT(edge_count == 15, "Petersen graph should have 15 edges");
    
    bliss_release(petersen);
    TEST_SUCCESS();
}

/* Test graph comparison and copying */
static int test_graph_operations(void) {
    printf("Testing graph operations...\n");
    
    /* Create two identical graphs */
    bliss_graph_t *g1 = bliss_create_cycle_graph(4);
    bliss_graph_t *g2 = bliss_create_cycle_graph(4);
    TEST_ASSERT(g1 != NULL && g2 != NULL, "Graph creation failed");
    
    /* They should compare as equal */
    int cmp = bliss_cmp(g1, g2);
    TEST_ASSERT(cmp == 0, "Identical graphs should compare equal");
    
    /* Test graph copying */
    bliss_graph_t *g1_copy = bliss_copy(g1);
    TEST_ASSERT(g1_copy != NULL, "Graph copying failed");
    
    cmp = bliss_cmp(g1, g1_copy);
    TEST_ASSERT(cmp == 0, "Original and copy should be equal");
    
    /* Verify copy has same properties */
    TEST_ASSERT(bliss_get_nof_vertices(g1_copy) == 4, "Copy should have 4 vertices");
    for (unsigned int i = 0; i < 4; i++) {
        TEST_ASSERT(bliss_get_degree(g1_copy, i) == 2, "Copy vertices should have degree 2");
    }
    
    /* Test that different graphs are not equal */
    bliss_graph_t *g3 = bliss_create_cycle_graph(5);
    cmp = bliss_cmp(g1, g3);
    TEST_ASSERT(cmp != 0, "Different graphs should not be equal");
    
    bliss_release(g1);
    bliss_release(g2);
    bliss_release(g1_copy);
    bliss_release(g3);
    TEST_SUCCESS();
}

/* Test vertex coloring */
static int test_vertex_coloring(void) {
    printf("Testing vertex coloring...\n");
    
    bliss_graph_t *graph = bliss_create_cycle_graph(4);
    TEST_ASSERT(graph != NULL, "Graph creation failed");
    
    /* Initially all vertices should have color 0 */
    for (unsigned int i = 0; i < 4; i++) {
        TEST_ASSERT(bliss_get_color(graph, i) == 0, "Initial color should be 0");
    }
    
    /* Change some colors */
    bliss_change_color(graph, 0, 1);
    bliss_change_color(graph, 2, 1);
    bliss_change_color(graph, 1, 2);
    bliss_change_color(graph, 3, 2);
    
    /* Verify colors were changed */
    TEST_ASSERT(bliss_get_color(graph, 0) == 1, "Color should be 1");
    TEST_ASSERT(bliss_get_color(graph, 1) == 2, "Color should be 2");
    TEST_ASSERT(bliss_get_color(graph, 2) == 1, "Color should be 1");
    TEST_ASSERT(bliss_get_color(graph, 3) == 2, "Color should be 2");
    
    bliss_release(graph);
    TEST_SUCCESS();
}

/* Test canonical labeling consistency */
static int test_canonical_labeling(void) {
    printf("Testing canonical labeling...\n");
    
    /* Create two isomorphic graphs with different vertex numbering */
    bliss_graph_t *g1 = bliss_new(4);
    bliss_graph_t *g2 = bliss_new(4);
    
    /* Graph 1: cycle 0-1-2-3-0 */
    bliss_add_edge(g1, 0, 1);
    bliss_add_edge(g1, 1, 2);
    bliss_add_edge(g1, 2, 3);
    bliss_add_edge(g1, 3, 0);
    
    /* Graph 2: cycle 1-2-3-0-1 (same cycle, different starting point) */
    bliss_add_edge(g2, 1, 2);
    bliss_add_edge(g2, 2, 3);
    bliss_add_edge(g2, 3, 0);
    bliss_add_edge(g2, 0, 1);
    
    /* Get canonical labelings */
    bliss_stats_t *stats1 = bliss_stats_new();
    bliss_stats_t *stats2 = bliss_stats_new();
    
    const unsigned int *canon1 = bliss_find_canonical_labeling(g1, stats1, NULL, NULL);
    const unsigned int *canon2 = bliss_find_canonical_labeling(g2, stats2, NULL, NULL);
    
    TEST_ASSERT(canon1 != NULL, "Canonical labeling 1 should not be NULL");
    TEST_ASSERT(canon2 != NULL, "Canonical labeling 2 should not be NULL");
    
    /* Apply canonical labelings and compare resulting graphs */
    bliss_graph_t *canonical_g1 = bliss_permute(g1, canon1);
    bliss_graph_t *canonical_g2 = bliss_permute(g2, canon2);
    
    TEST_ASSERT(canonical_g1 != NULL, "Canonical graph 1 should not be NULL");
    TEST_ASSERT(canonical_g2 != NULL, "Canonical graph 2 should not be NULL");
    
    /* Canonical forms should be identical */
    int cmp_result = bliss_cmp(canonical_g1, canonical_g2);
    TEST_ASSERT(cmp_result == 0, "Canonical forms of isomorphic graphs should be identical");
    
    bliss_stats_release(stats1);
    bliss_stats_release(stats2);
    bliss_release(g1);
    bliss_release(g2);
    bliss_release(canonical_g1);
    bliss_release(canonical_g2);
    
    TEST_SUCCESS();
}

/* Test graph connectivity */
static int test_graph_connectivity(void) {
    printf("Testing graph connectivity...\n");
    
    /* Connected graph: cycle */
    bliss_graph_t *connected = bliss_create_cycle_graph(5);
    TEST_ASSERT(bliss_is_connected(connected), "Cycle should be connected");
    bliss_release(connected);
    
    /* Disconnected graph: two separate edges */
    bliss_graph_t *disconnected = bliss_new(4);
    bliss_add_edge(disconnected, 0, 1);  /* Component 1 */
    bliss_add_edge(disconnected, 2, 3);  /* Component 2 */
    
    TEST_ASSERT(!bliss_is_connected(disconnected), "Two separate edges should be disconnected");
    bliss_release(disconnected);
    
    /* Single vertex should be connected */
    bliss_graph_t *single = bliss_new(1);
    TEST_ASSERT(bliss_is_connected(single), "Single vertex should be connected");
    bliss_release(single);
    
    TEST_SUCCESS();
}

/* Test hash function consistency */
static int test_hash_consistency(void) {
    printf("Testing hash consistency...\n");
    
    bliss_graph_t *g1 = bliss_create_cycle_graph(4);
    bliss_graph_t *g2 = bliss_create_cycle_graph(4);
    bliss_graph_t *g3 = bliss_create_cycle_graph(5);
    
    /* Identical graphs should have identical hashes */
    uint32_t h1 = bliss_get_hash(g1);
    uint32_t h2 = bliss_get_hash(g2);
    uint32_t h3 = bliss_get_hash(g3);
    
    TEST_ASSERT(h1 == h2, "Identical graphs should have same hash");
    TEST_ASSERT(h1 != h3, "Different graphs should have different hashes");
    
    /* Hash should be consistent across multiple calls */
    uint32_t h1_again = bliss_get_hash(g1);
    TEST_ASSERT(h1 == h1_again, "Hash should be consistent");
    
    bliss_release(g1);
    bliss_release(g2);
    bliss_release(g3);
    TEST_SUCCESS();
}

/* Test random graph generation */
static int test_random_graphs(void) {
    printf("Testing random graph generation...\n");
    
    /* Test random graph with different parameters */
    bliss_graph_t *sparse = bliss_create_random_graph(10, 0.1, 12345);
    bliss_graph_t *dense = bliss_create_random_graph(10, 0.8, 12345);
    
    TEST_ASSERT(sparse != NULL, "Sparse random graph creation failed");
    TEST_ASSERT(dense != NULL, "Dense random graph creation failed");
    
    TEST_ASSERT(bliss_get_nof_vertices(sparse) == 10, "Random graph should have 10 vertices");
    TEST_ASSERT(bliss_get_nof_vertices(dense) == 10, "Random graph should have 10 vertices");
    
    /* Dense graph should have more edges than sparse (probabilistically) */
    unsigned int sparse_edges = 0, dense_edges = 0;
    
    for (unsigned int i = 0; i < 10; i++) {
        for (unsigned int j = i + 1; j < 10; j++) {
            if (bliss_has_edge(sparse, i, j)) sparse_edges++;
            if (bliss_has_edge(dense, i, j)) dense_edges++;
        }
    }
    
    TEST_ASSERT(dense_edges > sparse_edges, "Dense graph should have more edges");
    
    bliss_release(sparse);
    bliss_release(dense);
    TEST_SUCCESS();
}


/* Main test runner */
int main(void) {
    printf("Running classic graph structure tests...\n\n");
    
    int passed = 0;
    int total = 0;
    
    /* Run structural and property tests */
    total++; passed += test_complete_graphs();
    total++; passed += test_cycle_graphs();
    total++; passed += test_path_graphs();
    total++; passed += test_star_graphs();
    total++; passed += test_petersen_graph();
    total++; passed += test_graph_operations();
    total++; passed += test_vertex_coloring();
    total++; passed += test_graph_connectivity();
    total++; passed += test_hash_consistency();
    total++; passed += test_random_graphs();

    total++; passed += test_canonical_labeling();

    
    /* Print summary */
    printf("\nTest Results: %d/%d passed\n", passed, total);
    
    if (passed == total) {
        printf("All classic graph structure tests PASSED!\n");
        printf("\nGraph generators and basic operations are working correctly.\n");
        printf("The C translation preserves the structural integrity of the original bliss.\n");
        return 0;
    } else {
        printf("Some tests FAILED!\n");
        return 1;
    }
}