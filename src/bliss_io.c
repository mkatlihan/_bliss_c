/*
 * bliss_io.c - Input/Output functions for DIMACS format and DOT format
 * 
 * This file handles file I/O operations and graph format conversions.
 * Separated for better code organization and compiler optimization.
 */

#include "bliss.h"
#include <ctype.h>
#include <time.h>

 /* ===================================================================
  * VERBOSE OUTPUT UTILITIES
  * =================================================================== */
void DPRINTF(const char* format, ...) {
  int verbose_level = bliss_get_verbose_level();
  FILE* verbose_file = bliss_get_verbose_file();
  if (verbose_level > 0) {
    va_list args;
    va_start(args, format);
    vfprintf(verbose_file ? verbose_file : stdout, format, args);
    va_end(args);
  }
}
void DPRINTF_IF(int level, const char* format, ...) {
  if (bliss_get_verbose_level() >= level) {
    va_list args;
    va_start(args, format);
    vfprintf(bliss_get_verbose_file() ? bliss_get_verbose_file() : stdout, format, args);
    va_end(args);
  }
}

/* ===================================================================
 * ADVANCED I/O FUNCTIONS AND FORMAT EXTENSIONS
 * =================================================================== */

/* Read graph from file with automatic format detection */
bliss_graph_t *bliss_read_graph(const char *filename) {
    if (BLISS_UNLIKELY(!filename)) {
        return NULL;
    }
    
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        return NULL;
    }
    
    bliss_graph_t *graph = bliss_read_dimacs(fp);
    fclose(fp);
    
    return graph;
}

/* Write graph to file in DIMACS format */
bool bliss_write_graph(const bliss_graph_t *graph, const char *filename) {
    if (BLISS_UNLIKELY(!graph || !filename)) {
        return false;
    }
    
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        return false;
    }
    
    bliss_write_dimacs(graph, fp);
    fclose(fp);
    
    return true;
}

/* Write graph to file in DOT format */
bool bliss_write_dot_file(const bliss_graph_t *graph, const char *filename) {
    if (BLISS_UNLIKELY(!graph || !filename)) {
        return false;
    }
    
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        return false;
    }
    
    bliss_write_dot(graph, fp);
    fclose(fp);
    
    return true;
}

/* ===================================================================
 * GRAPH VALIDATION AND CONSISTENCY CHECKING
 * =================================================================== */

BLISS_HOT
bool bliss_is_valid_graph(const bliss_graph_t *graph) {
    if (!graph) {
        return false;
    }
    
    /* Check basic consistency */
    if (graph->num_vertices > graph->vertex_capacity) {
        return false;
    }
    
    /* Verify adjacency list consistency */
    unsigned int edge_count = 0;
    
    for (unsigned int i = 0; i < graph->num_vertices; i++) {
        /* Check adjacency list bounds */
        if (graph->adj_list_sizes[i] > graph->adj_list_caps[i]) {
            return false;
        }
        
        /* Check that all neighbors are valid vertices */
        for (unsigned int j = 0; j < graph->adj_list_sizes[i]; j++) {
            unsigned int neighbor = graph->adj_lists[i][j];
            if (neighbor >= graph->num_vertices) {
                return false;
            }
        }
        
        edge_count += graph->adj_list_sizes[i];
    }
    
    /* For undirected graphs, each edge is counted twice */
    if (!graph->is_directed) {
        if (edge_count % 2 != 0) {
            return false; /* Odd number of directed edges in undirected graph */
        }
        edge_count /= 2;
    }
    
    /* Check edge count consistency */
    if (edge_count != graph->num_edges) {
        return false;
    }
    
    /* For directed graphs, verify incoming adjacency lists */
    if (graph->is_directed) {
        for (unsigned int i = 0; i < graph->num_vertices; i++) {
            if (graph->in_adj_list_sizes[i] > graph->in_adj_list_caps[i]) {
                return false;
            }
            
            for (unsigned int j = 0; j < graph->in_adj_list_sizes[i]; j++) {
                unsigned int source = graph->in_adj_lists[i][j];
                if (source >= graph->num_vertices) {
                    return false;
                }
                
                /* Verify that corresponding outgoing edge exists */
                bool found = false;
                for (unsigned int k = 0; k < graph->adj_list_sizes[source]; k++) {
                    if (graph->adj_lists[source][k] == i) {
                        found = true;
                        break;
                    }
                }
                
                if (!found) {
                    return false;
                }
            }
        }
    }
    
    return true;
}

/* ===================================================================
 * GRAPH STATISTICS AND PROPERTIES
 * =================================================================== */

/* Get degree of a vertex */
unsigned int bliss_get_degree(const bliss_graph_t *graph, unsigned int vertex) {
    if (BLISS_UNLIKELY(!graph || vertex >= graph->num_vertices)) {
        return 0;
    }
    
    return graph->adj_list_sizes[vertex];
}

/* Get in-degree of a vertex (for directed graphs) */
unsigned int bliss_get_in_degree(const bliss_graph_t *graph, unsigned int vertex) {
    if (BLISS_UNLIKELY(!graph || vertex >= graph->num_vertices || !graph->is_directed)) {
        return 0;
    }
    
    return graph->in_adj_list_sizes[vertex];
}

/* Get out-degree of a vertex (for directed graphs) */
unsigned int bliss_get_out_degree(const bliss_graph_t *graph, unsigned int vertex) {
    if (BLISS_UNLIKELY(!graph || vertex >= graph->num_vertices)) {
        return 0;
    }
    
    return graph->adj_list_sizes[vertex];
}

/* Check if there's an edge between two vertices */
bool bliss_has_edge(const bliss_graph_t *graph, unsigned int v1, unsigned int v2) {
    if (BLISS_UNLIKELY(!graph || v1 >= graph->num_vertices || v2 >= graph->num_vertices)) {
        return false;
    }
    
    /* Search in adjacency list of v1 */
    for (unsigned int i = 0; i < graph->adj_list_sizes[v1]; i++) {
        if (graph->adj_lists[v1][i] == v2) {
            return true;
        }
    }
    
    return false;
}

/* Get adjacency list of a vertex */
const unsigned int *bliss_get_adjacency_list(const bliss_graph_t *graph, 
                                             unsigned int vertex, 
                                             unsigned int *list_size) {
    if (BLISS_UNLIKELY(!graph || vertex >= graph->num_vertices)) {
        if (list_size) *list_size = 0;
        return NULL;
    }
    
    if (list_size) {
        *list_size = graph->adj_list_sizes[vertex];
    }
    
    return graph->adj_lists[vertex];
}

/* ===================================================================
 * MEMORY-EFFICIENT BULK OPERATIONS
 * =================================================================== */

/* Add multiple edges efficiently */
void bliss_add_edges_bulk(bliss_graph_t *graph, 
                          const unsigned int *edge_list, 
                          unsigned int num_edges) {
    if (BLISS_UNLIKELY(!graph || !edge_list)) {
        return;
    }
    
    /* Pre-allocate space for efficiency */
    for (unsigned int i = 0; i < num_edges; i++) {
        unsigned int v1 = edge_list[2 * i];
        unsigned int v2 = edge_list[2 * i + 1];
        
        if (v1 < graph->num_vertices && v2 < graph->num_vertices) {
            /* Ensure capacity for v1 */
            if (graph->adj_list_sizes[v1] >= graph->adj_list_caps[v1]) {
                graph->adj_list_caps[v1] = max_uint(INITIAL_ADJ_CAPACITY, 
                                                    graph->adj_list_caps[v1] * 2);
                graph->adj_lists[v1] = bliss_realloc(graph->adj_lists[v1], 
                                                     graph->adj_list_caps[v1] * sizeof(unsigned int));
            }
            
            /* For undirected graphs, ensure capacity for v2 as well */
            if (!graph->is_directed && v1 != v2) {
                if (graph->adj_list_sizes[v2] >= graph->adj_list_caps[v2]) {
                    graph->adj_list_caps[v2] = max_uint(INITIAL_ADJ_CAPACITY, 
                                                        graph->adj_list_caps[v2] * 2);
                    graph->adj_lists[v2] = bliss_realloc(graph->adj_lists[v2], 
                                                         graph->adj_list_caps[v2] * sizeof(unsigned int));
                }
            }
        }
    }
    
    /* Add all edges */
    for (unsigned int i = 0; i < num_edges; i++) {
        unsigned int v1 = edge_list[2 * i];
        unsigned int v2 = edge_list[2 * i + 1];
        bliss_add_edge(graph, v1, v2);
    }
}

/* Set vertex colors efficiently */
void bliss_set_colors_bulk(bliss_graph_t *graph, 
                           const unsigned int *colors, 
                           unsigned int num_vertices) {
    if (BLISS_UNLIKELY(!graph || !colors)) {
        return;
    }
    
    unsigned int count = min_uint(num_vertices, graph->num_vertices);
    memcpy(graph->vertex_colors, colors, count * sizeof(unsigned int));
    
    graph->hash_valid = false;
    graph->canonical_labeling_valid = false;
}

/* ===================================================================
 * ERROR HANDLING AND DIAGNOSTICS
 * =================================================================== */

/* Get human-readable error message */
const char *bliss_error_string(bliss_error_t error) {
    switch (error) {
        case BLISS_OK:
            return "Success";
        case BLISS_ERROR_INVALID_ARGUMENT:
            return "Invalid argument";
        case BLISS_ERROR_OUT_OF_MEMORY:
            return "Out of memory";
        case BLISS_ERROR_IO:
            return "I/O error";
        case BLISS_ERROR_INVALID_GRAPH:
            return "Invalid graph";
        default:
            return "Unknown error";
    }
}

/* Print graph statistics to file */
void bliss_print_graph_info(const bliss_graph_t *graph, FILE *fp) {
    if (BLISS_UNLIKELY(!graph || !fp)) {
        return;
    }
    
    fprintf(fp, "Graph Information:\n");
    fprintf(fp, "  Type: %s\n", graph->is_directed ? "Directed" : "Undirected");
    fprintf(fp, "  Vertices: %u\n", graph->num_vertices);
    fprintf(fp, "  Edges: %u\n", graph->num_edges);
    fprintf(fp, "  Hash: 0x%08x (%s)\n", 
            graph->graph_hash, 
            graph->hash_valid ? "valid" : "invalid");
    
    /* Color distribution */
    unsigned int max_color = 0;
    for (unsigned int i = 0; i < graph->num_vertices; i++) {
        if (graph->vertex_colors[i] > max_color) {
            max_color = graph->vertex_colors[i];
        }
    }
    
    fprintf(fp, "  Color classes: %u\n", max_color + 1);
    
    /* Degree statistics */
    if (graph->num_vertices > 0) {
        unsigned int min_degree = graph->adj_list_sizes[0];
        unsigned int max_degree = graph->adj_list_sizes[0];
        unsigned long total_degree = 0;
        
        for (unsigned int i = 0; i < graph->num_vertices; i++) {
            unsigned int degree = graph->adj_list_sizes[i];
            if (degree < min_degree) min_degree = degree;
            if (degree > max_degree) max_degree = degree;
            total_degree += degree;
        }
        
        fprintf(fp, "  Degree range: [%u, %u]\n", min_degree, max_degree);
        fprintf(fp, "  Average degree: %.2f\n", 
                (double)total_degree / graph->num_vertices);
    }
    
    fprintf(fp, "  Memory usage: ~%zu bytes\n", 
            sizeof(bliss_graph_t) + 
            graph->vertex_capacity * sizeof(unsigned int) + /* colors */
            graph->vertex_capacity * 3 * sizeof(unsigned int) + /* adj metadata */
            graph->num_edges * sizeof(unsigned int) * (graph->is_directed ? 2 : 2)); /* adj lists */
}

/* Validate adjacency list consistency and repair if possible */
bliss_error_t bliss_validate_and_repair(bliss_graph_t *graph) {
    if (!graph) {
        return BLISS_ERROR_INVALID_ARGUMENT;
    }
    
    bool repaired = false;
    
    /* Remove invalid neighbors */
    for (unsigned int i = 0; i < graph->num_vertices; i++) {
        unsigned int write_pos = 0;
        
        for (unsigned int read_pos = 0; read_pos < graph->adj_list_sizes[i]; read_pos++) {
            unsigned int neighbor = graph->adj_lists[i][read_pos];
            
            if (neighbor < graph->num_vertices) {
                graph->adj_lists[i][write_pos++] = neighbor;
            } else {
                repaired = true;
            }
        }
        
        graph->adj_list_sizes[i] = write_pos;
    }
    
    /* Recount edges */
    unsigned int new_edge_count = 0;
    for (unsigned int i = 0; i < graph->num_vertices; i++) {
        new_edge_count += graph->adj_list_sizes[i];
    }
    
    if (!graph->is_directed) {
        new_edge_count /= 2;
    }
    
    if (new_edge_count != graph->num_edges) {
        graph->num_edges = new_edge_count;
        repaired = true;
    }
    
    /* Invalidate cached data if repairs were made */
    if (repaired) {
        graph->hash_valid = false;
        graph->canonical_labeling_valid = false;
    }
    
    return bliss_is_valid_graph(graph) ? BLISS_OK : BLISS_ERROR_INVALID_GRAPH;
}

/* ===================================================================
 * SPECIALIZED GRAPH GENERATORS FOR TESTING
 * =================================================================== */

/* Create a complete graph (clique) */
bliss_graph_t *bliss_create_complete_graph(unsigned int n) {
    bliss_graph_t *graph = bliss_new(n);
    if (!graph) return NULL;
    
    for (unsigned int i = 0; i < n; i++) {
        for (unsigned int j = i + 1; j < n; j++) {
            bliss_add_edge(graph, i, j);
        }
    }
    
    return graph;
}

/* Create a cycle graph */
bliss_graph_t *bliss_create_cycle_graph(unsigned int n) {
    if (n < 3) return NULL;
    
    bliss_graph_t *graph = bliss_new(n);
    if (!graph) return NULL;
    
    for (unsigned int i = 0; i < n; i++) {
        bliss_add_edge(graph, i, (i + 1) % n);
    }
    
    return graph;
}

/* Create a path graph */
bliss_graph_t *bliss_create_path_graph(unsigned int n) {
    if (n < 2) return NULL;
    
    bliss_graph_t *graph = bliss_new(n);
    if (!graph) return NULL;
    
    for (unsigned int i = 0; i < n - 1; i++) {
        bliss_add_edge(graph, i, i + 1);
    }
    
    return graph;
}

/* Create a star graph */
bliss_graph_t *bliss_create_star_graph(unsigned int n) {
    if (n < 2) return NULL;
    
    bliss_graph_t *graph = bliss_new(n);
    if (!graph) return NULL;
    
    for (unsigned int i = 1; i < n; i++) {
        bliss_add_edge(graph, 0, i);
    }
    
    return graph;
}

/* Create Petersen graph (classic test case) */
bliss_graph_t *bliss_create_petersen_graph(void) {
    bliss_graph_t *graph = bliss_new(10);
    if (!graph) return NULL;
    
    /* Outer pentagon */
    for (unsigned int i = 0; i < 5; i++) {
        bliss_add_edge(graph, i, (i + 1) % 5);
    }
    
    /* Inner pentagram */
    for (unsigned int i = 0; i < 5; i++) {
        bliss_add_edge(graph, i + 5, ((i + 2) % 5) + 5);
    }
    
    /* Connections between outer and inner */
    for (unsigned int i = 0; i < 5; i++) {
        bliss_add_edge(graph, i, i + 5);
    }
    
    return graph;
}

/* Create a random graph (for testing) */
bliss_graph_t *bliss_create_random_graph(unsigned int n, double edge_prob, unsigned int seed) {
    bliss_graph_t *graph = bliss_new(n);
    if (!graph) return NULL;
    
    /* Simple linear congruential generator for reproducible results */
    unsigned int rng_state = seed;
    
    for (unsigned int i = 0; i < n; i++) {
        for (unsigned int j = i + 1; j < n; j++) {
            rng_state = rng_state * 1103515245U + 12345U;
            double rand_val = (double)(rng_state % 32768) / 32768.0;
            
            if (rand_val < edge_prob) {
                bliss_add_edge(graph, i, j);
            }
        }
    }
    
    return graph;
}

/* ===================================================================
 * GRAPH ANALYSIS UTILITIES
 * =================================================================== */

/* Check if graph is connected (for undirected graphs) */
bool bliss_is_connected(const bliss_graph_t *graph) {
    if (!graph || graph->num_vertices == 0) {
        return true;
    }
    
    if (graph->is_directed) {
        return false; /* Not implemented for directed graphs */
    }
    
    bool *visited = bliss_malloc(graph->num_vertices * sizeof(bool));
    memset(visited, 0, graph->num_vertices * sizeof(bool));
    
    /* DFS from vertex 0 */
    unsigned int *stack = bliss_malloc(graph->num_vertices * sizeof(unsigned int));
    unsigned int stack_size = 0;
    
    stack[stack_size++] = 0;
    visited[0] = true;
    unsigned int visited_count = 1;
    
    while (stack_size > 0) {
        unsigned int current = stack[--stack_size];
        
        for (unsigned int i = 0; i < graph->adj_list_sizes[current]; i++) {
            unsigned int neighbor = graph->adj_lists[current][i];
            if (!visited[neighbor]) {
                visited[neighbor] = true;
                visited_count++;
                stack[stack_size++] = neighbor;
            }
        }
    }
    
    bliss_free(visited);
    bliss_free(stack);
    
    return visited_count == graph->num_vertices;
}

/* Count number of triangles in graph */
unsigned int bliss_count_triangles(const bliss_graph_t *graph) {
    if (!graph || graph->is_directed) {
        return 0;
    }
    
    unsigned int triangle_count = 0;
    
    for (unsigned int i = 0; i < graph->num_vertices; i++) {
        for (unsigned int j_idx = 0; j_idx < graph->adj_list_sizes[i]; j_idx++) {
            unsigned int j = graph->adj_lists[i][j_idx];
            if (j <= i) continue; /* Avoid double counting */
            
            for (unsigned int k_idx = 0; k_idx < graph->adj_list_sizes[j]; k_idx++) {
                unsigned int k = graph->adj_lists[j][k_idx];
                if (k <= j) continue;
                
                /* Check if edge (i,k) exists */
                if (bliss_has_edge(graph, i, k)) {
                    triangle_count++;
                }
            }
        }
    }
    
    return triangle_count;
}

/* Compute clustering coefficient */
double bliss_clustering_coefficient(const bliss_graph_t *graph) {
    if (!graph || graph->is_directed || graph->num_vertices < 3) {
        return 0.0;
    }
    
    double total_coefficient = 0.0;
    unsigned int vertices_with_degree_ge_2 = 0;
    
    for (unsigned int i = 0; i < graph->num_vertices; i++) {
        unsigned int degree = graph->adj_list_sizes[i];
        if (degree < 2) continue;
        
        vertices_with_degree_ge_2++;
        
        /* Count triangles involving vertex i */
        unsigned int triangles = 0;
        for (unsigned int j_idx = 0; j_idx < degree; j_idx++) {
            unsigned int j = graph->adj_lists[i][j_idx];
            
            for (unsigned int k_idx = j_idx + 1; k_idx < degree; k_idx++) {
                unsigned int k = graph->adj_lists[i][k_idx];
                
                if (bliss_has_edge(graph, j, k)) {
                    triangles++;
                }
            }
        }
        
        /* Local clustering coefficient */
        double possible_edges = (double)(degree * (degree - 1)) / 2.0;
        total_coefficient += triangles / possible_edges;
    }
    
    return vertices_with_degree_ge_2 > 0 ? 
           total_coefficient / vertices_with_degree_ge_2 : 0.0;
}

/* ===================================================================
 * PERFORMANCE BENCHMARKING UTILITIES
 * =================================================================== */

/* Measure time for automorphism computation */
double bliss_benchmark_automorphisms(bliss_graph_t *graph, unsigned int iterations) {
    if (!graph || iterations == 0) {
        return 0.0;
    }
    
    clock_t start = clock();
    
    for (unsigned int i = 0; i < iterations; i++) {
        bliss_stats_t *stats = bliss_stats_new();
        bliss_find_automorphisms(graph, stats, NULL, NULL);
        bliss_stats_release(stats);
    }
    
    clock_t end = clock();
    
    return ((double)(end - start)) / CLOCKS_PER_SEC / iterations;
}

/* Memory usage estimation */
size_t bliss_estimate_memory_usage(const bliss_graph_t *graph) {
    if (!graph) return 0;
    
    size_t base_size = sizeof(bliss_graph_t);
    size_t vertex_arrays = graph->vertex_capacity * sizeof(unsigned int) * 4; /* colors, sizes, caps, lists */
    
    size_t adjacency_memory = 0;
    for (unsigned int i = 0; i < graph->num_vertices; i++) {
        adjacency_memory += graph->adj_list_caps[i] * sizeof(unsigned int);
        if (graph->is_directed) {
            adjacency_memory += graph->in_adj_list_caps[i] * sizeof(unsigned int);
        }
    }
    
    size_t canonical_memory = graph->canonical_labeling ? 
                              graph->num_vertices * sizeof(unsigned int) : 0;
    
    return base_size + vertex_arrays + adjacency_memory + canonical_memory;
} 

/* DIMACS FORMAT I/O FUNCTIONS
 * =================================================================== */

BLISS_HOT
static void skip_whitespace_and_comments(FILE *fp) {
    int c;
    while ((c = fgetc(fp)) != EOF) {
        if (c == 'c') {
            /* Skip comment line */
            while ((c = fgetc(fp)) != EOF && c != '\n') {
                /* consume comment */
            }
        } else if (isspace(c)) {
            /* Skip whitespace */
            continue;
        } else {
            /* Put back non-whitespace, non-comment character */
            ungetc(c, fp);
            break;
        }
    }
}

BLISS_HOT
static bool read_dimacs_problem_line(FILE *fp, bool *is_directed, 
                                     unsigned int *num_vertices, 
                                     unsigned int *num_edges) {
    skip_whitespace_and_comments(fp);
    
    char line_type;
    char format[32];
    
    if (fscanf(fp, "%c %31s %u %u", &line_type, format, num_vertices, num_edges) != 4) {
        return false;
    }
    
    if (line_type != 'p') {
        return false;
    }
    
    /* Check format - bliss supports both 'edge' and 'col' variants */
    if (strcmp(format, "edge") == 0) {
        *is_directed = false;
    } else if (strcmp(format, "col") == 0) {
        *is_directed = false;
    } else if (strcmp(format, "arc") == 0) {
        *is_directed = true;
    } else {
        /* Unknown format */
        return false;
    }
    
    return true;
}

BLISS_HOT
bliss_graph_t *bliss_read_dimacs(FILE *fp) {
    if (BLISS_UNLIKELY(!fp)) {
        return NULL;
    }
    
    bool is_directed;
    unsigned int num_vertices, num_edges;
    
    if (!read_dimacs_problem_line(fp, &is_directed, &num_vertices, &num_edges)) {
        return NULL;
    }
    
    bliss_graph_t *graph = is_directed ? 
        bliss_new_directed(num_vertices) : 
        bliss_new(num_vertices);
    
    if (!graph) {
        return NULL;
    }
    
    char line_type;
    unsigned int vertex, color;
    unsigned int v1, v2;
    
    /* Read vertex colors and edges */
    while (fscanf(fp, " %c", &line_type) == 1) {
        switch (line_type) {
            case 'c':
                /* Comment line - skip to end of line */
                while (fgetc(fp) != '\n' && !feof(fp)) {
                    /* consume comment */
                }
                break;
                
            case 'n':
            case 'v':
                /* Vertex color definition: n vertex color */
                if (fscanf(fp, "%u %u", &vertex, &color) != 2) {
                    bliss_release(graph);
                    return NULL;
                }
                
                /* Convert from 1-based to 0-based indexing */
                if (vertex == 0 || vertex > num_vertices) {
                    bliss_release(graph);
                    return NULL;
                }
                vertex--;
                
                bliss_change_color(graph, vertex, color);
                break;
                
            case 'e':
            case 'a':
                /* Edge definition: e v1 v2 */
                if (fscanf(fp, "%u %u", &v1, &v2) != 2) {
                    bliss_release(graph);
                    return NULL;
                }
                
                /* Convert from 1-based to 0-based indexing */
                if (v1 == 0 || v1 > num_vertices || v2 == 0 || v2 > num_vertices) {
                    bliss_release(graph);
                    return NULL;
                }
                v1--;
                v2--;
                
                bliss_add_edge(graph, v1, v2);
                break;
                
            default:
                /* Unknown line type - ignore */
                while (fgetc(fp) != '\n' && !feof(fp)) {
                    /* consume line */
                }
                break;
        }
    }
    
    return graph;
}

BLISS_HOT
void bliss_write_dimacs(const bliss_graph_t *graph, FILE *fp) {
    if (BLISS_UNLIKELY(!graph || !fp)) {
        return;
    }
    
    /* Write header comment */
    fprintf(fp, "c Graph with %u vertices and %u edges\n", 
            graph->num_vertices, graph->num_edges);
    fprintf(fp, "c Generated by bliss C implementation\n");
    
    /* Write problem line */
    fprintf(fp, "p %s %u %u\n", 
            graph->is_directed ? "arc" : "edge",
            graph->num_vertices, 
            graph->num_edges);
    
    /* Write vertex colors (only non-zero colors) */
    for (unsigned int i = 0; i < graph->num_vertices; i++) {
        if (graph->vertex_colors[i] != 0) {
            fprintf(fp, "n %u %u\n", i + 1, graph->vertex_colors[i]);
        }
    }
    
    /* Write edges */
    for (unsigned int i = 0; i < graph->num_vertices; i++) {
        for (unsigned int j = 0; j < graph->adj_list_sizes[i]; j++) {
            unsigned int neighbor = graph->adj_lists[i][j];
            
            /* For undirected graphs, only write each edge once */
            if (!graph->is_directed && i > neighbor) {
                continue;
            }
            
            fprintf(fp, "%c %u %u\n", 
                    graph->is_directed ? 'a' : 'e',
                    i + 1, neighbor + 1);
        }
    }
}

/* ===================================================================
 * DOT FORMAT OUTPUT (Graphviz)
 * =================================================================== */

BLISS_HOT
void bliss_write_dot(const bliss_graph_t *graph, FILE *fp) {
    if (BLISS_UNLIKELY(!graph || !fp)) {
        return;
    }
    
    /* Write DOT header */
    fprintf(fp, "%s G {\n", graph->is_directed ? "digraph" : "graph");
    fprintf(fp, "  rankdir=LR;\n");
    fprintf(fp, "  node [shape=circle];\n");
    
    /* Define vertex colors */
    const char *colors[] = {
        "white", "lightblue", "lightgreen", "lightcoral", 
        "lightyellow", "lightpink", "lightgray", "lightcyan"
    };
    const int num_colors = sizeof(colors) / sizeof(colors[0]);
    
    /* Write vertices with colors */
    for (unsigned int i = 0; i < graph->num_vertices; i++) {
        const char *color = colors[graph->vertex_colors[i] % num_colors];
        fprintf(fp, "  %u [label=\"%u\" fillcolor=\"%s\" style=\"filled\"];\n", 
                i, i, color);
    }
    
    /* Write edges */
    const char *edge_op = graph->is_directed ? " -> " : " -- ";
    
    for (unsigned int i = 0; i < graph->num_vertices; i++) {
        for (unsigned int j = 0; j < graph->adj_list_sizes[i]; j++) {
            unsigned int neighbor = graph->adj_lists[i][j];
            
            /* For undirected graphs, only write each edge once */
            if (!graph->is_directed && i > neighbor) {
                continue;
            }
            
            fprintf(fp, "  %u%s%u;\n", i, edge_op, neighbor);
        }
    }
    
    fprintf(fp, "}\n");
}

/* ===================================================================
 * GRAPH CONSTRUCTION HELPERS
 * =================================================================== */

void bliss_add_vertex(bliss_graph_t *graph, unsigned int color) {
    if (BLISS_UNLIKELY(!graph)) {
        return;
    }
    
    /* Expand capacity if needed */
    if (graph->num_vertices >= graph->vertex_capacity) {
        unsigned int new_capacity = graph->vertex_capacity * 2;
        
        graph->vertex_colors = bliss_realloc(graph->vertex_colors, 
                                             new_capacity * sizeof(unsigned int));
        graph->adj_lists = bliss_realloc(graph->adj_lists, 
                                         new_capacity * sizeof(unsigned int*));
        graph->adj_list_sizes = bliss_realloc(graph->adj_list_sizes, 
                                              new_capacity * sizeof(unsigned int));
        graph->adj_list_caps = bliss_realloc(graph->adj_list_caps, 
                                             new_capacity * sizeof(unsigned int));
        
        if (graph->is_directed) {
            graph->in_adj_lists = bliss_realloc(graph->in_adj_lists, 
                                                new_capacity * sizeof(unsigned int*));
            graph->in_adj_list_sizes = bliss_realloc(graph->in_adj_list_sizes, 
                                                     new_capacity * sizeof(unsigned int));
            graph->in_adj_list_caps = bliss_realloc(graph->in_adj_list_caps, 
                                                    new_capacity * sizeof(unsigned int));
        }
        
        /* Initialize new entries */
        for (unsigned int i = graph->vertex_capacity; i < new_capacity; i++) {
            graph->vertex_colors[i] = 0;
            graph->adj_lists[i] = NULL;
            graph->adj_list_sizes[i] = 0;
            graph->adj_list_caps[i] = 0;
            
            if (graph->is_directed) {
                graph->in_adj_lists[i] = NULL;
                graph->in_adj_list_sizes[i] = 0;
                graph->in_adj_list_caps[i] = 0;
            }
        }
        
        graph->vertex_capacity = new_capacity;
    }
    
    /* Set color for new vertex */
    graph->vertex_colors[graph->num_vertices] = color;
    graph->num_vertices++;
    
    /* Invalidate cached data */
    graph->hash_valid = false;
    graph->canonical_labeling_valid = false;
}

/* ===================================================================
 * GRAPH COPYING
 * =================================================================== */

bliss_graph_t *bliss_copy(const bliss_graph_t *graph) {
    if (BLISS_UNLIKELY(!graph)) {
        return NULL;
    }
    
    bliss_graph_t *copy = graph->is_directed ? 
        bliss_new_directed(graph->num_vertices) : 
        bliss_new(graph->num_vertices);
    
    if (!copy) {
        return NULL;
    }
    
    /* Copy basic properties */
    copy->num_edges = graph->num_edges;
    
    /* Copy vertex colors */
    memcpy(copy->vertex_colors, graph->vertex_colors, 
           graph->num_vertices * sizeof(unsigned int));
    
    /* Copy adjacency lists */
    for (unsigned int i = 0; i < graph->num_vertices; i++) {
        if (graph->adj_list_sizes[i] > 0) {
            copy->adj_list_sizes[i] = graph->adj_list_sizes[i];
            copy->adj_list_caps[i] = graph->adj_list_caps[i];
            copy->adj_lists[i] = bliss_malloc(copy->adj_list_caps[i] * sizeof(unsigned int));
            memcpy(copy->adj_lists[i], graph->adj_lists[i], 
                   graph->adj_list_sizes[i] * sizeof(unsigned int));
        }
        
        if (graph->is_directed && graph->in_adj_list_sizes[i] > 0) {
            copy->in_adj_list_sizes[i] = graph->in_adj_list_sizes[i];
            copy->in_adj_list_caps[i] = graph->in_adj_list_caps[i];
            copy->in_adj_lists[i] = bliss_malloc(copy->in_adj_list_caps[i] * sizeof(unsigned int));
            memcpy(copy->in_adj_lists[i], graph->in_adj_lists[i], 
                   graph->in_adj_list_sizes[i] * sizeof(unsigned int));
        }
    }
    
    /* Copy configuration */
    copy->splitting_heuristic = graph->splitting_heuristic;
    copy->use_component_recursion = graph->use_component_recursion;
    copy->use_failure_recording = graph->use_failure_recording;
    copy->use_long_prune = graph->use_long_prune;
    
    /* Copy hash if valid */
    if (graph->hash_valid) {
        copy->graph_hash = graph->graph_hash;
        copy->hash_valid = true;
    }
    
    return copy;
}

/* ===================================================================
 */