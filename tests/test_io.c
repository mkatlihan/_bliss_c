/*
 * test_io.c - Tests for DIMACS and DOT format I/O
 * 
 * Tests reading and writing graphs in various formats,
 * including edge cases and error handling.
 */

#include "bliss.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>
#include <time.h>

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

/* Create a temporary file with given content */
static FILE* create_temp_file(const char* content) {
    static int temp_counter = 0;
    char temp_filename[256];
    
    /* Create a unique temporary filename */
    snprintf(temp_filename, sizeof(temp_filename), "bliss_temp_%d_%d.tmp", 
             (int)time(NULL), temp_counter++);
    
    FILE* fp = fopen(temp_filename, "w+");
    if (fp && content) {
        fwrite(content, 1, strlen(content), fp);
        rewind(fp);
    }
    
    /* Note: We don't delete the file here since we need to read from it.
     * The file will be cleaned up by the OS or can be manually deleted. */
    
    return fp;
}

/* Helper to clean up temporary file */
static void cleanup_temp_file(FILE* fp, const char* content) {
    (void)content; /* Unused parameter */
    if (fp) {
        fclose(fp);
    }
}

/* Test basic DIMACS reading */
static int test_dimacs_reading(void) {
    printf("Testing DIMACS format reading...\n");
    
    /* Simple undirected graph */
    const char* dimacs_content = 
        "c This is a comment\n"
        "c Another comment\n"
        "p edge 4 5\n"
        "n 1 0\n"
        "n 2 1\n"
        "n 3 0\n"
        "n 4 2\n"
        "e 1 2\n"
        "e 2 3\n"
        "e 3 4\n"
        "e 4 1\n"
        "e 1 3\n";
    
    FILE* fp = create_temp_file(dimacs_content);
    TEST_ASSERT(fp != NULL, "Failed to create temporary file");

    bliss_graph_t* graph = bliss_read_dimacs(fp);
    fclose(fp);
    
    TEST_ASSERT(graph != NULL, "Failed to read DIMACS graph");
    TEST_ASSERT(bliss_get_nof_vertices(graph) == 4, "Should have 4 vertices");
    TEST_ASSERT(!bliss_is_directed(graph), "Should be undirected");
    
    /* Check vertex colors (DIMACS uses 1-based indexing) */
    TEST_ASSERT(bliss_get_color(graph, 0) == 0, "Vertex 0 should have color 0");
    TEST_ASSERT(bliss_get_color(graph, 1) == 1, "Vertex 1 should have color 1");
    TEST_ASSERT(bliss_get_color(graph, 2) == 0, "Vertex 2 should have color 0");
    TEST_ASSERT(bliss_get_color(graph, 3) == 2, "Vertex 3 should have color 2");
    
    /* Check edges */
    TEST_ASSERT(bliss_has_edge(graph, 0, 1), "Edge (0,1) should exist");
    TEST_ASSERT(bliss_has_edge(graph, 1, 2), "Edge (1,2) should exist");
    TEST_ASSERT(bliss_has_edge(graph, 2, 3), "Edge (2,3) should exist");
    TEST_ASSERT(bliss_has_edge(graph, 3, 0), "Edge (3,0) should exist");
    TEST_ASSERT(bliss_has_edge(graph, 0, 2), "Edge (0,2) should exist");
    
        
    bliss_release(graph);
    TEST_SUCCESS();
}

/* Test directed graph DIMACS reading */
static int test_dimacs_directed(void) {
    printf("Testing directed DIMACS format...\n");
    
    const char* dimacs_content = 
        "c Directed graph test\n"
        "p arc 3 3\n"
        "a 1 2\n"
        "a 2 3\n"
        "a 3 1\n";
    
    FILE* fp = create_temp_file(dimacs_content);
    TEST_ASSERT(fp != NULL, "Failed to create temporary file");
    
    bliss_graph_t* graph = bliss_read_dimacs(fp);
    cleanup_temp_file(fp, dimacs_content);
    
    TEST_ASSERT(graph != NULL, "Failed to read directed DIMACS graph");
    TEST_ASSERT(bliss_get_nof_vertices(graph) == 3, "Should have 3 vertices");
    TEST_ASSERT(bliss_is_directed(graph), "Should be directed");
    
    /* Check directed edges */
    TEST_ASSERT(bliss_has_edge(graph, 0, 1), "Edge (0,1) should exist");
    TEST_ASSERT(bliss_has_edge(graph, 1, 2), "Edge (1,2) should exist");
    TEST_ASSERT(bliss_has_edge(graph, 2, 0), "Edge (2,0) should exist");
    
    /* Check that reverse edges don't exist (directed) */
    TEST_ASSERT(!bliss_has_edge(graph, 1, 0), "Reverse edge (1,0) should not exist");
    TEST_ASSERT(!bliss_has_edge(graph, 2, 1), "Reverse edge (2,1) should not exist");
    TEST_ASSERT(!bliss_has_edge(graph, 0, 2), "Reverse edge (0,2) should not exist");
    
    bliss_release(graph);
    TEST_SUCCESS();
}

/* Test DIMACS writing */
static int test_dimacs_writing(void) {
    printf("Testing DIMACS format writing...\n");
    
    /* Create a simple graph */
    bliss_graph_t* graph = bliss_new(3);
    bliss_add_edge(graph, 0, 1);
    bliss_add_edge(graph, 1, 2);
    bliss_change_color(graph, 0, 5);
    bliss_change_color(graph, 2, 3);
    
    /* Write to temporary file using a simpler approach */
    char temp_filename[256];
    snprintf(temp_filename, sizeof(temp_filename), "test_write_%d.dimacs", (int)time(NULL));
    
    FILE* fp = fopen(temp_filename, "w");
    TEST_ASSERT(fp != NULL, "Failed to create write test file");
    
    bliss_write_dimacs(graph, fp);
    fclose(fp);
    
    /* Read back and verify */
    fp = fopen(temp_filename, "r");
    TEST_ASSERT(fp != NULL, "Failed to reopen written file");
    
    bliss_graph_t* read_graph = bliss_read_dimacs(fp);
    fclose(fp);
    
    /* Clean up temp file */
    unlink(temp_filename);
    
    TEST_ASSERT(read_graph != NULL, "Failed to read back written graph");
    TEST_ASSERT(bliss_get_nof_vertices(read_graph) == 3, "Vertex count mismatch");
    TEST_ASSERT(bliss_has_edge(read_graph, 0, 1), "Edge (0,1) missing");
    TEST_ASSERT(bliss_has_edge(read_graph, 1, 2), "Edge (1,2) missing");
    TEST_ASSERT(bliss_get_color(read_graph, 0) == 5, "Color mismatch for vertex 0");
    TEST_ASSERT(bliss_get_color(read_graph, 1) == 0, "Color mismatch for vertex 1");
    TEST_ASSERT(bliss_get_color(read_graph, 2) == 3, "Color mismatch for vertex 2");
    
    bliss_release(graph);
    bliss_release(read_graph);
    TEST_SUCCESS();
}

/* Test DOT format writing */
static int test_dot_writing(void) {
    printf("Testing DOT format writing...\n");
    
    /* Create a small graph */
    bliss_graph_t* graph = bliss_new(4);
    bliss_add_edge(graph, 0, 1);
    bliss_add_edge(graph, 1, 2);
    bliss_add_edge(graph, 2, 3);
    bliss_add_edge(graph, 3, 0);
    
    bliss_change_color(graph, 0, 1);
    bliss_change_color(graph, 2, 1);
    
    /* Write to temporary file */
    char temp_filename[256];
    snprintf(temp_filename, sizeof(temp_filename), "test_dot_%d.dot", (int)time(NULL));
    
    FILE* fp = fopen(temp_filename, "w");
    TEST_ASSERT(fp != NULL, "Failed to create DOT test file");
    
    bliss_write_dot(graph, fp);
    fclose(fp);
    
    /* Read content and check basic structure */
    fp = fopen(temp_filename, "r");
    TEST_ASSERT(fp != NULL, "Failed to reopen DOT file");
    
    char buffer[1024];
    size_t bytes_read = fread(buffer, 1, sizeof(buffer) - 1, fp);
    buffer[bytes_read] = '\0';
    fclose(fp);
    
    /* Clean up */
    unlink(temp_filename);
    
    /* Check for expected DOT structure */
    TEST_ASSERT(strstr(buffer, "graph G") != NULL, "Should contain 'graph G'");
    TEST_ASSERT(strstr(buffer, "rankdir=LR") != NULL, "Should contain rankdir");
    TEST_ASSERT(strstr(buffer, "0 --") != NULL, "Should contain edge from vertex 0");
    TEST_ASSERT(strstr(buffer, "fillcolor") != NULL, "Should contain color information");
    
    bliss_release(graph);
    TEST_SUCCESS();
}

/* Test directed DOT format */
static int test_dot_directed(void) {
    printf("Testing directed DOT format...\n");
    
    bliss_graph_t* graph = bliss_new_directed(3);
    bliss_add_edge(graph, 0, 1);
    bliss_add_edge(graph, 1, 2);
    bliss_add_edge(graph, 2, 0);
    
    char temp_filename[256];
    snprintf(temp_filename, sizeof(temp_filename), "test_directed_%d.dot", (int)time(NULL));
    
    FILE* fp = fopen(temp_filename, "w");
    TEST_ASSERT(fp != NULL, "Failed to create directed DOT file");
    
    bliss_write_dot(graph, fp);
    fclose(fp);
    
    fp = fopen(temp_filename, "r");
    TEST_ASSERT(fp != NULL, "Failed to reopen directed DOT file");
    
    char buffer[1024];
    size_t bytes_read = fread(buffer, 1, sizeof(buffer) - 1, fp);
    buffer[bytes_read] = '\0';
    fclose(fp);
    
    unlink(temp_filename);
    
    /* Check for directed graph structure */
    TEST_ASSERT(strstr(buffer, "digraph G") != NULL, "Should contain 'digraph G'");
    TEST_ASSERT(strstr(buffer, " -> ") != NULL, "Should contain directed edge operator");
    
    bliss_release(graph);
    TEST_SUCCESS();
}

/* Test error handling in DIMACS reading */
static int test_dimacs_error_handling(void) {
    printf("Testing DIMACS error handling...\n");
    
    /* Test invalid problem line */
    const char* invalid_content1 = "p invalid 3 2\n";
    FILE* fp1 = create_temp_file(invalid_content1);
    bliss_graph_t* graph1 = bliss_read_dimacs(fp1);
    fclose(fp1);
    TEST_ASSERT(graph1 == NULL, "Should reject invalid format");
    
    /* Test invalid vertex numbers */
    const char* invalid_content2 = 
        "p edge 3 2\n"
        "e 1 5\n"  /* Vertex 5 doesn't exist */
        "e 2 3\n";
    FILE* fp2 = create_temp_file(invalid_content2);
    bliss_graph_t* graph2 = bliss_read_dimacs(fp2);
    fclose(fp2);
    TEST_ASSERT(graph2 == NULL, "Should reject invalid vertex numbers");
    
    /* Test malformed lines */
    const char* invalid_content3 = 
        "p edge 2 1\n"
        "e 1\n";  /* Missing second vertex */
    FILE* fp3 = create_temp_file(invalid_content3);
    bliss_graph_t* graph3 = bliss_read_dimacs(fp3);
    fclose(fp3);
    TEST_ASSERT(graph3 == NULL, "Should reject malformed edge lines");
    
    TEST_SUCCESS();
}

/* Test file I/O functions */
static int test_file_io_functions(void) {
    printf("Testing file I/O convenience functions...\n");
    
    /* Create a test graph */
    bliss_graph_t* original = bliss_create_cycle_graph(5);
    bliss_change_color(original, 0, 1);
    bliss_change_color(original, 2, 2);
    
    /* Write to file */
    const char* filename = "test_graph.dimacs";
    bool write_success = bliss_write_graph(original, filename);
    TEST_ASSERT(write_success, "Writing graph to file should succeed");
    
    /* Read back from file */
    bliss_graph_t* loaded = bliss_read_graph(filename);
    TEST_ASSERT(loaded != NULL, "Reading graph from file should succeed");
    
    /* Compare graphs */
    TEST_ASSERT(bliss_get_nof_vertices(loaded) == bliss_get_nof_vertices(original), 
                "Vertex count should match");
    TEST_ASSERT(bliss_get_color(loaded, 0) == bliss_get_color(original, 0), 
                "Colors should match");
    TEST_ASSERT(bliss_get_color(loaded, 2) == bliss_get_color(original, 2), 
                "Colors should match");
    
    /* Test DOT file writing */
    const char* dot_filename = "test_graph.dot";
    bool dot_success = bliss_write_dot_file(original, dot_filename);
    TEST_ASSERT(dot_success, "Writing DOT file should succeed");
    
    /* Clean up */
    unlink(filename);
    unlink(dot_filename);
    bliss_release(original);
    bliss_release(loaded);
    
    TEST_SUCCESS();
}

/* Test bulk operations */
static int test_bulk_operations(void) {
    printf("Testing bulk operations...\n");
    
    bliss_graph_t* graph = bliss_new(6);
    
    /* Test bulk edge addition */
    unsigned int edges[] = {0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 0};
    bliss_add_edges_bulk(graph, edges, 6);
    
    /* Verify edges were added */
    TEST_ASSERT(bliss_has_edge(graph, 0, 1), "Bulk edge (0,1) should exist");
    TEST_ASSERT(bliss_has_edge(graph, 2, 3), "Bulk edge (2,3) should exist");
    TEST_ASSERT(bliss_has_edge(graph, 5, 0), "Bulk edge (5,0) should exist");
    
    /* Test bulk color setting */
    unsigned int colors[] = {1, 2, 1, 2, 1, 2};
    bliss_set_colors_bulk(graph, colors, 6);
    
    /* Verify colors */
    TEST_ASSERT(bliss_get_color(graph, 0) == 1, "Bulk color 0 should be 1");
    TEST_ASSERT(bliss_get_color(graph, 1) == 2, "Bulk color 1 should be 2");
    TEST_ASSERT(bliss_get_color(graph, 5) == 2, "Bulk color 5 should be 2");
    
    bliss_release(graph);
    TEST_SUCCESS();
}

/* Test large graph I/O */
static int test_large_graph_io(void) {
    printf("Testing large graph I/O...\n");
    
    /* Create a moderately large graph */
    const unsigned int n = 100;
    bliss_graph_t* graph = bliss_new(n);
    
    /* Create a random-ish graph */
    for (unsigned int i = 0; i < n; i++) {
        bliss_add_edge(graph, i, (i + 1) % n);  /* Cycle */
        if (i < n/2) {
            bliss_add_edge(graph, i, i + n/2);   /* Additional edges */
        }
        bliss_change_color(graph, i, i % 5);    /* 5 different colors */
    }
    
    /* Write and read back using file-based approach */
    char temp_filename[256];
    snprintf(temp_filename, sizeof(temp_filename), "test_large_%d.dimacs", (int)time(NULL));
    
    FILE* fp = fopen(temp_filename, "w");
    TEST_ASSERT(fp != NULL, "Failed to create large graph test file");
    
    bliss_write_dimacs(graph, fp);
    fclose(fp);
    
    fp = fopen(temp_filename, "r");
    TEST_ASSERT(fp != NULL, "Failed to reopen large graph file");
    
    bliss_graph_t* loaded = bliss_read_dimacs(fp);
    fclose(fp);
    
    unlink(temp_filename);
    
    TEST_ASSERT(loaded != NULL, "Failed to load large graph");
    TEST_ASSERT(bliss_get_nof_vertices(loaded) == n, "Vertex count mismatch");
    
    /* Spot check some properties */
    TEST_ASSERT(bliss_has_edge(loaded, 0, 1), "Cycle edge should exist");
    TEST_ASSERT(bliss_has_edge(loaded, n-1, 0), "Cycle edge should exist");
    TEST_ASSERT(bliss_get_color(loaded, 10) == 10 % 5, "Color should match");
    
    bliss_release(graph);
    bliss_release(loaded);
    TEST_SUCCESS();
}

/* Test validation and repair */
static int test_validation_and_repair(void) {
    printf("Testing graph validation and repair...\n");
    
    bliss_graph_t* graph = bliss_new(3);
    bliss_add_edge(graph, 0, 1);
    bliss_add_edge(graph, 1, 2);
    
    /* Graph should be valid initially */
    TEST_ASSERT(bliss_is_valid_graph(graph), "Graph should be valid");
    
    bliss_error_t result = bliss_validate_and_repair(graph);
    TEST_ASSERT(result == BLISS_OK, "Validation should succeed");
    
    bliss_release(graph);
    TEST_SUCCESS();
}

/* Test memory usage estimation */
static int test_memory_estimation(void) {
    printf("Testing memory usage estimation...\n");
    
    bliss_graph_t* graph = bliss_new(10);
    
    /* Add some edges */
    for (unsigned int i = 0; i < 10; i++) {
        bliss_add_edge(graph, i, (i + 1) % 10);
    }
    
    size_t memory = bliss_estimate_memory_usage(graph);
    TEST_ASSERT(memory > 0, "Memory usage should be positive");
    
    printf("  Estimated memory usage: %zu bytes\n", memory);
    
    bliss_release(graph);
    TEST_SUCCESS();
}

/* Main test runner */
int main(void) {
    printf("Running I/O format tests...\n\n");
    
    int passed = 0;
    int total = 0;
    
    /* Run all tests */
    total++; passed += test_dimacs_reading();
    total++; passed += test_dimacs_directed();
    total++; passed += test_dimacs_writing();
    total++; passed += test_dot_writing();
    total++; passed += test_dot_directed();
    total++; passed += test_dimacs_error_handling();
    total++; passed += test_file_io_functions();
    total++; passed += test_bulk_operations();
    total++; passed += test_large_graph_io();
    total++; passed += test_validation_and_repair();
    total++; passed += test_memory_estimation();
    
    /* Print summary */
    printf("\nTest Results: %d/%d passed\n", passed, total);
    
    if (passed == total) {
        printf("All I/O tests PASSED!\n");
        return 0;
    } else {
        printf("Some tests FAILED!\n");
        return 1;
    }
}