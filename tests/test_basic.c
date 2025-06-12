/*
 * test_basic.c - Basic functionality tests for bliss C library
 * 
 * Tests fundamental operations like graph creation, edge addition,
 * and basic properties.
 */

#include "../bliss.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

/* Test framework macros */
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

/* Test graph creation and basic properties */
static int test_graph_creation(void) {
    bliss_graph_t *graph = bliss_new(5);
    TEST_ASSERT(graph != NULL, "Graph creation failed");
    TEST_ASSERT(bliss_get_nof_vertices(graph) == 5, "Vertex count incorrect");
    TEST_ASSERT(!bliss_is_directed(graph), "Should be undirected");
    
    bliss_release(graph);
    
    /* Test directed graph */
    bliss_graph_t *digraph = bliss_new_directed(3);
    TEST_ASSERT(digraph != NULL, "Directed graph creation failed");
    TEST_ASSERT(bliss_get_nof_vertices(digraph) == 3, "Directed vertex count incorrect");
    TEST_ASSERT(bliss_is_directed(digraph), "Should be directed");
    
    bliss_release(digraph);
    TEST_SUCCESS();
}

/* Test edge addition and properties */
static int test_edge_operations(void) {
    bliss_graph_t *graph = bliss_new(4);
    
    /* Add edges to form a cycle */
    bliss_add_edge(graph, 0, 1);
    bliss_add_edge(graph, 1, 2);
    bliss_add_edge(graph, 2, 3);
    bliss_add_edge(graph, 3, 0);
    
    /* Test edge existence */
    TEST_ASSERT(bliss_has_edge(graph, 0, 1), "Edge (0,1) should exist");
    TEST_ASSERT(bliss_has_edge(graph, 1, 0), "Edge (1,0) should exist (undirected)");
    TEST_ASSERT(bliss_has_edge(graph, 2, 3), "Edge (2,3) should exist");
    TEST_ASSERT(!bliss_has_edge(graph, 0, 2), "Edge (0,2) should not exist");
    
    /* Test degrees */
    TEST_ASSERT(bliss_get_degree(graph, 0) == 2, "Vertex 0 should have degree 2");
    TEST_ASSERT(bliss_get_degree(graph, 1) == 2, "Vertex 1 should have degree 2");
    
    bliss_release(graph);
    TEST_SUCCESS();
}

/* Test vertex colors */
static int test_vertex_colors(void) {
    bliss_graph_t *graph = bliss_new(4);
    
    /* Test default colors */
    TEST_ASSERT(bliss_get_color(graph, 0) == 0, "Default color should be 0");
    TEST_ASSERT(bliss_get_color(graph, 3) == 0, "Default color should be 0");
    
    /* Change colors */
    bliss_change_color(graph, 0, 1);
    bliss_change_color(graph, 1, 1);
    bliss_change_color(graph, 2, 2);
    bliss_change_color(graph, 3, 2);
    
    TEST_ASSERT(bliss_get_color(graph, 0) == 1, "Color should be 1");
    TEST_ASSERT(bliss_get_color(graph, 1) == 1, "Color should be 1");
    TEST_ASSERT(bliss_get_color(graph, 2) == 2, "Color should be 2");
    TEST_ASSERT(bliss_get_color(graph, 3) == 2, "Color should be 2");
    
    bliss_release(graph);
    TEST_SUCCESS();
}

/* Test graph copying */
static int test_graph_copying(void) {
    bliss_graph_t *original = bliss_new(3);
    
    /* Build a simple graph */
    bliss_add_edge(original, 0, 1);
    bliss_add_edge(original, 1, 2);
    bliss_change_color(original, 0, 5);
    bliss_change_color(original, 2, 7);
    
    /* Copy the graph */
    bliss_graph_t *copy = bliss_copy(original);
    TEST_ASSERT(copy != NULL, "Graph copying failed");
    
    /* Verify copy properties */
    TEST_ASSERT(bliss_get_nof_vertices(copy) == bliss_get_nof_vertices(original), 
                "Vertex count should match");
    TEST_ASSERT(bliss_is_directed(copy) == bliss_is_directed(original), 
                "Direction property should match");
    
    /* Verify edges */
    TEST_ASSERT(bliss_has_edge(copy, 0, 1), "Copied edge (0,1) should exist");
    TEST_ASSERT(bliss_has_edge(copy, 1, 2), "Copied edge (1,2) should exist");
    TEST_ASSERT(!bliss_has_edge(copy, 0, 2), "Non-edge (0,2) should not exist");
    
    /* Verify colors */
    TEST_ASSERT(bliss_get_color(copy, 0) == 5, "Color should be copied");
    TEST_ASSERT(bliss_get_color(copy, 1) == 0, "Default color should be copied");
    TEST_ASSERT(bliss_get_color(copy, 2) == 7, "Color should be copied");
    
    /* Modify original and verify copy is independent */
    bliss_change_color(original, 0, 99);
    TEST_ASSERT(bliss_get_color(copy, 0) == 5, "Copy should be independent");
    
    bliss_release(original);
    bliss_release(copy);
    TEST_SUCCESS();
}

/* Test graph comparison */
static int test_graph_comparison(void) {
    bliss_graph_t *g1 = bliss_new(3);
    bliss_graph_t *g2 = bliss_new(3);
    bliss_graph_t *g3 = bliss_new(4);
    
    /* Initially empty graphs should be equal */
    TEST_ASSERT(bliss_cmp(g1, g2) == 0, "Empty graphs should be equal");
    
    /* Different sizes should not be equal */
    TEST_ASSERT(bliss_cmp(g1, g3) != 0, "Different sized graphs should not be equal");
    
    /* Add same edges to both */
    bliss_add_edge(g1, 0, 1);
    bliss_add_edge(g1, 1, 2);
    bliss_add_edge(g2, 0, 1);
    bliss_add_edge(g2, 1, 2);
    
    TEST_ASSERT(bliss_cmp(g1, g2) == 0, "Graphs with same edges should be equal");
    
    /* Add different edge to g2 */
    bliss_add_edge(g2, 0, 2);
    TEST_ASSERT(bliss_cmp(g1, g2) != 0, "Graphs with different edges should not be equal");
    
    bliss_release(g1);
    bliss_release(g2);
    bliss_release(g3);
    TEST_SUCCESS();
}

/* Test hash function */
static int test_hash_function(void) {
    bliss_graph_t *g1 = bliss_new(3);
    bliss_graph_t *g2 = bliss_new(3);
    
    /* Empty graphs should have same hash */
    uint32_t h1 = bliss_get_hash(g1);
    uint32_t h2 = bliss_get_hash(g2);
    TEST_ASSERT(h1 == h2, "Empty graphs should have same hash");
    
    /* Add edge and check hash changes */
    bliss_add_edge(g1, 0, 1);
    uint32_t h1_new = bliss_get_hash(g1);
    TEST_ASSERT(h1 != h1_new, "Hash should change when graph changes");
    
    /* Same modification to g2 should give same hash */
    bliss_add_edge(g2, 0, 1);
    h2 = bliss_get_hash(g2);
    TEST_ASSERT(h1_new == h2, "Identical graphs should have same hash");
    
    bliss_release(g1);
    bliss_release(g2);
    TEST_SUCCESS();
}

/* Test statistics structure */
static int test_statistics(void) {
    bliss_stats_t *stats = bliss_stats_new();
    TEST_ASSERT(stats != NULL, "Statistics creation failed");
    
    /* Check initial values */
    TEST_ASSERT(bliss_stats_get_nof_nodes(stats) == 0, "Initial node count should be 0");
    TEST_ASSERT(bliss_stats_get_nof_leaf_nodes(stats) == 0, "Initial leaf count should be 0");
    TEST_ASSERT(bliss_stats_get_nof_generators(stats) == 0, "Initial generator count should be 0");
    TEST_ASSERT(bliss_stats_get_group_size_approx(stats) == 0.0, "Initial group size should be 0");
    
    bliss_stats_release(stats);
    TEST_SUCCESS();
}

/* Test error handling */
static int test_error_handling(void) {
    /* Test NULL pointer handling */
    TEST_ASSERT(bliss_get_nof_vertices(NULL) == 0, "NULL graph should return 0 vertices");
    TEST_ASSERT(bliss_get_color(NULL, 0) == 0, "NULL graph should return color 0");
    TEST_ASSERT(bliss_get_hash(NULL) == 0, "NULL graph should return hash 0");
    TEST_ASSERT(!bliss_is_directed(NULL), "NULL graph should not be directed");
    
    /* Test invalid vertex access */
    bliss_graph_t *graph = bliss_new(3);
    TEST_ASSERT(bliss_get_color(graph, 5) == 0, "Invalid vertex should return color 0");
    TEST_ASSERT(bliss_get_degree(graph, 5) == 0, "Invalid vertex should return degree 0");
    TEST_ASSERT(!bliss_has_edge(graph, 5, 0), "Invalid vertex should not have edges");
    
    /* Test edge addition with invalid vertices */
    bliss_add_edge(graph, 0, 5); /* Should be ignored */
    bliss_add_edge(graph, 5, 0); /* Should be ignored */
    TEST_ASSERT(bliss_get_degree(graph, 0) == 0, "Invalid edges should be ignored");
    
    bliss_release(graph);
    TEST_SUCCESS();
}

/* Test memory management */
static int test_memory_management(void) {
    /* Create and destroy many graphs to test for leaks */
    for (int i = 0; i < 100; i++) {
        bliss_graph_t *graph = bliss_new(10);
        
        /* Add some edges */
        for (unsigned int j = 0; j < 10; j++) {
            bliss_add_edge(graph, j, (j + 1) % 10);
        }
        
        /* Change some colors */
        for (unsigned int j = 0; j < 10; j++) {
            bliss_change_color(graph, j, j % 3);
        }
        
        bliss_release(graph);
    }
    
    TEST_SUCCESS();
}

/* Test dynamic vertex addition */
static int test_dynamic_vertices(void) {
    bliss_graph_t *graph = bliss_new(0);
    TEST_ASSERT(bliss_get_nof_vertices(graph) == 0, "Initial graph should be empty");
    
    /* Add vertices dynamically */
    bliss_add_vertex(graph, 0);
    bliss_add_vertex(graph, 1);
    bliss_add_vertex(graph, 0);
    
    TEST_ASSERT(bliss_get_nof_vertices(graph) == 3, "Should have 3 vertices");
    TEST_ASSERT(bliss_get_color(graph, 0) == 0, "Vertex 0 should have color 0");
    TEST_ASSERT(bliss_get_color(graph, 1) == 1, "Vertex 1 should have color 1");
    TEST_ASSERT(bliss_get_color(graph, 2) == 0, "Vertex 2 should have color 0");
    
    /* Add edges between dynamically added vertices */
    bliss_add_edge(graph, 0, 1);
    bliss_add_edge(graph, 1, 2);
    
    TEST_ASSERT(bliss_has_edge(graph, 0, 1), "Edge should exist");
    TEST_ASSERT(bliss_has_edge(graph, 1, 2), "Edge should exist");
    
    bliss_release(graph);
    TEST_SUCCESS();
}

/* Test adjacency list access */
static int test_adjacency_lists(void) {
    bliss_graph_t *graph = bliss_new(4);
    
    /* Create a star graph */
    bliss_add_edge(graph, 0, 1);
    bliss_add_edge(graph, 0, 2);
    bliss_add_edge(graph, 0, 3);
    
    /* Test adjacency list of center vertex */
    unsigned int list_size;
    const unsigned int *adj_list = bliss_get_adjacency_list(graph, 0, &list_size);
    TEST_ASSERT(adj_list != NULL, "Adjacency list should not be NULL");
    TEST_ASSERT(list_size == 3, "Center vertex should have 3 neighbors");
    
    /* Verify neighbors (order may vary) */
    bool found[3] = {false, false, false};
    for (unsigned int i = 0; i < list_size; i++) {
        if (adj_list[i] == 1) found[0] = true;
        else if (adj_list[i] == 2) found[1] = true;
        else if (adj_list[i] == 3) found[2] = true;
    }
    TEST_ASSERT(found[0] && found[1] && found[2], "All neighbors should be found");
    
    /* Test adjacency list of leaf vertex */
    adj_list = bliss_get_adjacency_list(graph, 1, &list_size);
    TEST_ASSERT(list_size == 1, "Leaf vertex should have 1 neighbor");
    TEST_ASSERT(adj_list[0] == 0, "Leaf should be connected to center");
    
    bliss_release(graph);
    TEST_SUCCESS();
}

/* Main test runner */
int main(void) {
    printf("Running basic functionality tests...\n\n");
    
    int passed = 0;
    int total = 0;
    
    /* Run all tests */
    total++; passed += test_graph_creation();
    total++; passed += test_edge_operations();
    total++; passed += test_vertex_colors();
    total++; passed += test_graph_copying();
    total++; passed += test_graph_comparison();
    total++; passed += test_hash_function();
    total++; passed += test_statistics();
    total++; passed += test_error_handling();
    total++; passed += test_memory_management();
    total++; passed += test_dynamic_vertices();
    total++; passed += test_adjacency_lists();
    
    /* Print summary */
    printf("\nTest Results: %d/%d passed\n", passed, total);
    
    if (passed == total) {
        printf("All basic tests PASSED!\n");
        return 0;
    } else {
        printf("Some tests FAILED!\n");
        return 1;
    }
} 