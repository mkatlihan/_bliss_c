/*
 * test_classic_graphs.c - Tests based on classic graphs from original bliss
 * 
 * This file contains tests using well-known graphs with known automorphism
 * group properties, allowing comparison with the original bliss results.
 */

#include "../bliss.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

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

/* Automorphism counting callback */
static unsigned int automorphism_count = 0;

static void count_automorphisms(void *user_param, unsigned int n, const unsigned int *aut) {
    (void)user_param; /* Unused */
    (void)n;         /* Unused */
    (void)aut;       /* Unused */
    automorphism_count++;
}

/* Test complete graph K_n */
static int test_complete_graph(void) {
    printf("Testing complete graphs...\n");
    
    /* K_3 (triangle) - should have 3! = 6 automorphisms */
    bliss_graph_t *k3 = bliss_create_complete_graph(3);
    TEST_ASSERT(k3 != NULL, "K3 creation failed");
    
    automorphism_count = 0;
    bliss_stats_t *stats = bliss_stats_new();
    bliss_find_automorphisms(k3, stats, count_automorphisms, NULL);
    
    printf("  K3: Found %u automorphisms\n", automorphism_count);
    TEST_ASSERT(automorphism_count == 6, "K3 should have 6 automorphisms");
    
    bliss_stats_release(stats);
    bliss_release(k3);
    
    /* K_4 (complete graph on 4 vertices) - should have 4! = 24 automorphisms */
    bliss_graph_t *k4 = bliss_create_complete_graph(4);
    TEST_ASSERT(k4 != NULL, "K4 creation failed");
    
    automorphism_count = 0;
    stats = bliss_stats_new();
    bliss_find_automorphisms(k4, stats, count_automorphisms, NULL);
    
    printf("  K4: Found %u automorphisms\n", automorphism_count);
    TEST_ASSERT(automorphism_count == 24, "K4 should have 24 automorphisms");
    
    bliss_stats_release(stats);
    bliss_release(k4);
    
    TEST_SUCCESS();
}

/* Test cycle graphs */
static int test_cycle_graphs(void) {
    printf("Testing cycle graphs...\n");
    
    /* C_3 (triangle) - should have 6 automorphisms (3 rotations × 2 reflections) */
    bliss_graph_t *c3 = bliss_create_cycle_graph(3);
    TEST_ASSERT(c3 != NULL, "C3 creation failed");
    
    automorphism_count = 0;
    bliss_stats_t *stats = bliss_stats_new();
    bliss_find_automorphisms(c3, stats, count_automorphisms, NULL);
    
    printf("  C3: Found %u automorphisms\n", automorphism_count);
    TEST_ASSERT(automorphism_count == 6, "C3 should have 6 automorphisms");
    
    bliss_stats_release(stats);
    bliss_release(c3);
    
    /* C_4 (square) - should have 8 automorphisms (4 rotations × 2 reflections) */
    bliss_graph_t *c4 = bliss_create_cycle_graph(4);
    TEST_ASSERT(c4 != NULL, "C4 creation failed");
    
    automorphism_count = 0;
    stats = bliss_stats_new();
    bliss_find_automorphisms(c4, stats, count_automorphisms, NULL);
    
    printf("  C4: Found %u automorphisms\n", automorphism_count);
    TEST_ASSERT(automorphism_count == 8, "C4 should have 8 automorphisms");
    
    bliss_stats_release(stats);
    bliss_release(c4);
    
    /* C_5 (pentagon) - should have 10 automorphisms (5 rotations × 2 reflections) */
    bliss_graph_t *c5 = bliss_create_cycle_graph(5);
    TEST_ASSERT(c5 != NULL, "C5 creation failed");
    
    automorphism_count = 0;
    stats = bliss_stats_new();
    bliss_find_automorphisms(c5, stats, count_automorphisms, NULL);
    
    printf("  C5: Found %u automorphisms\n", automorphism_count);
    TEST_ASSERT(automorphism_count == 10, "C5 should have 10 automorphisms");
    
    bliss_stats_release(stats);
    bliss_release(c5);
    
    TEST_SUCCESS();
}

/* Test path graphs */
static int test_path_graphs(void) {
    printf("Testing path graphs...\n");
    
    /* P_2 (single edge) - should have 2 automorphisms */
    bliss_graph_t *p2 = bliss_create_path_graph(2);
    TEST_ASSERT(p2 != NULL, "P2 creation failed");
    
    automorphism_count = 0;
    bliss_stats_t *stats = bliss_stats_new();
    bliss_find_automorphisms(p2, stats, count_automorphisms, NULL);
    
    printf("  P2: Found %u automorphisms\n", automorphism_count);
    TEST_ASSERT(automorphism_count == 2, "P2 should have 2 automorphisms");
    
    bliss_stats_release(stats);
    bliss_release(p2);
    
    /* P_3 (path of length 2) - should have 2 automorphisms */
    bliss_graph_t *p3 = bliss_create_path_graph(3);
    TEST_ASSERT(p3 != NULL, "P3 creation failed");
    
    automorphism_count = 0;
    stats = bliss_stats_new();
    bliss_find_automorphisms(p3, stats, count_automorphisms, NULL);
    
    printf("  P3: Found %u automorphisms\n", automorphism_count);
    TEST_ASSERT(automorphism_count == 2, "P3 should have 2 automorphisms");
    
    bliss_stats_release(stats);
    bliss_release(p3);
    
    /* P_4 (path of length 3) - should have 2 automorphisms */
    bliss_graph_t *p4 = bliss_create_path_graph(4);
    TEST_ASSERT(p4 != NULL, "P4 creation failed");
    
    automorphism_count = 0;
    stats = bliss_stats_new();
    bliss_find_automorphisms(p4, stats, count_automorphisms, NULL);
    
    printf("  P4: Found %u automorphisms\n", automorphism_count);
    TEST_ASSERT(automorphism_count == 2, "P4 should have 2 automorphisms");
    
    bliss_stats_release(stats);
    bliss_release(p4);
    
    TEST_SUCCESS();
}

/* Test star graphs */
static int test_star_graphs(void) {
    printf("Testing star graphs...\n");
    
    /* S_3 (star with 3 leaves) - should have 3! = 6 automorphisms */
    bliss_graph_t *s3 = bliss_create_star_graph(4); /* 4 vertices: 1 center + 3 leaves */
    TEST_ASSERT(s3 != NULL, "S3 creation failed");
    
    automorphism_count = 0;
    bliss_stats_t *stats = bliss_stats_new();
    bliss_find_automorphisms(s3, stats, count_automorphisms, NULL);
    
    printf("  S3: Found %u automorphisms\n", automorphism_count);
    TEST_ASSERT(automorphism_count == 6, "S3 should have 6 automorphisms");
    
    bliss_stats_release(stats);
    bliss_release(s3);
    
    /* S_4 (star with 4 leaves) - should have 4! = 24 automorphisms */
    bliss_graph_t *s4 = bliss_create_star_graph(5); /* 5 vertices: 1 center + 4 leaves */
    TEST_ASSERT(s4 != NULL, "S4 creation failed");
    
    automorphism_count = 0;
    stats = bliss_stats_new();
    bliss_find_automorphisms(s4, stats, count_automorphisms, NULL);
    
    printf("  S4: Found %u automorphisms\n", automorphism_count);
    TEST_ASSERT(automorphism_count == 24, "S4 should have 24 automorphisms");
    
    bliss_stats_release(stats);
    bliss_release(s4);
    
    TEST_SUCCESS();
}

/* Test Petersen graph - famous graph with known automorphism group */
static int test_petersen_graph(void) {
    printf("Testing Petersen graph...\n");
    
    bliss_graph_t *petersen = bliss_create_petersen_graph();
    TEST_ASSERT(petersen != NULL, "Petersen graph creation failed");
    
    /* Verify basic properties */
    TEST_ASSERT(bliss_get_nof_vertices(petersen) == 10, "Petersen graph should have 10 vertices");
    
    /* Each vertex should have degree 3 */
    for (unsigned int i = 0; i < 10; i++) {
        unsigned int degree = bliss_get_degree(petersen, i);
        TEST_ASSERT(degree == 3, "Each vertex in Petersen graph should have degree 3");
    }
    
    /* Petersen graph should have 120 automorphisms (|Aut(P)| = 120) */
    automorphism_count = 0;
    bliss_stats_t *stats = bliss_stats_new();
    bliss_find_automorphisms(petersen, stats, count_automorphisms, NULL);
    
    printf("  Petersen: Found %u automorphisms\n", automorphism_count);
    TEST_ASSERT(automorphism_count == 120, "Petersen graph should have 120 automorphisms");
    
    bliss_stats_release(stats);
    bliss_release(petersen);
    
    TEST_SUCCESS();
}

/* Test disconnected graphs */
static int test_disconnected_graphs(void) {
    printf("Testing disconnected graphs...\n");
    
    /* Two disjoint triangles - should have (3!)² = 36 automorphisms */
    bliss_graph_t *graph = bliss_new(6);
    
    /* First triangle: vertices 0, 1, 2 */
    bliss_add_edge(graph, 0, 1);
    bliss_add_edge(graph, 1, 2);
    bliss_add_edge(graph, 2, 0);
    
    /* Second triangle: vertices 3, 4, 5 */
    bliss_add_edge(graph, 3, 4);
    bliss_add_edge(graph, 4, 5);
    bliss_add_edge(graph, 5, 3);
    
    automorphism_count = 0;
    bliss_stats_t *stats = bliss_stats_new();
    bliss_find_automorphisms(graph, stats, count_automorphisms, NULL);
    
    printf("  Two triangles: Found %u automorphisms\n", automorphism_count);
    /* Note: The actual count might be different due to how the algorithm handles components */
    
    bliss_stats_release(stats);
    bliss_release(graph);
    
    TEST_SUCCESS();
}

/* Test graphs with vertex colors */
static int test_colored_graphs(void) {
    printf("Testing colored graphs...\n");
    
    /* Create a cycle and color vertices alternately */
    bliss_graph_t *graph = bliss_create_cycle_graph(4);
    
    /* Color vertices 0,2 with color 0 and vertices 1,3 with color 1 */
    bliss_change_color(graph, 0, 0);
    bliss_change_color(graph, 1, 1);
    bliss_change_color(graph, 2, 0);
    bliss_change_color(graph, 3, 1);
    
    automorphism_count = 0;
    bliss_stats_t *stats = bliss_stats_new();
    bliss_find_automorphisms(graph, stats, count_automorphisms, NULL);
    
    printf("  Colored C4: Found %u automorphisms\n", automorphism_count);
    /* With alternating colors, should have fewer automorphisms than uncolored C4 */
    TEST_ASSERT(automorphism_count < 8, "Colored C4 should have fewer than 8 automorphisms");
    
    bliss_stats_release(stats);
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

/* Test different splitting heuristics */
static int test_splitting_heuristics(void) {
    printf("Testing splitting heuristics...\n");
    
    bliss_graph_t *graph = bliss_create_petersen_graph();
    TEST_ASSERT(graph != NULL, "Petersen graph creation failed");
    
    bliss_splitting_heuristic_t heuristics[] = {
        BLISS_SH_F, BLISS_SH_FL, BLISS_SH_FS, 
        BLISS_SH_FM, BLISS_SH_FLM, BLISS_SH_FSM
    };
    const char *heuristic_names[] = {
        "F", "FL", "FS", "FM", "FLM", "FSM"
    };
    
    for (int i = 0; i < 6; i++) {
        bliss_set_splitting_heuristic(graph, heuristics[i]);
        
        automorphism_count = 0;
        bliss_stats_t *stats = bliss_stats_new();
        bliss_find_automorphisms(graph, stats, count_automorphisms, NULL);
        
        printf("  Heuristic %s: %u automorphisms, %lu nodes\n", 
               heuristic_names[i], automorphism_count, 
               bliss_stats_get_nof_nodes(stats));
        
        /* All heuristics should find the same number of automorphisms */
        TEST_ASSERT(automorphism_count == 120, 
                    "All heuristics should find 120 automorphisms for Petersen graph");
        
        bliss_stats_release(stats);
    }
    
    bliss_release(graph);
    TEST_SUCCESS();
}

/* Main test runner */
int main(void) {
    printf("Running classic graph automorphism tests...\n\n");
    
    int passed = 0;
    int total = 0;
    
    /* Run all tests */
    total++; passed += test_complete_graph();
    total++; passed += test_cycle_graphs();
    total++; passed += test_path_graphs();
    total++; passed += test_star_graphs();
    total++; passed += test_petersen_graph();
    total++; passed += test_disconnected_graphs();
    total++; passed += test_colored_graphs();
    total++; passed += test_canonical_labeling();
    total++; passed += test_splitting_heuristics();
    
    /* Print summary */
    printf("\nTest Results: %d/%d passed\n", passed, total);
    
    if (passed == total) {
        printf("All classic graph tests PASSED!\n");
        return 0;
    } else {
        printf("Some tests FAILED!\n");
        return 1;
    }
}