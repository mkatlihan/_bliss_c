/*
 * bliss.h - Graph automorphism and canonical labeling library
 * 
 * C translation of the bliss C++ library
 * Copyright (c) 2025 - Plain C Translation
 * Based on original bliss by Tommi Junttila and Petteri Kaski
 * Licensed under GNU LGPL v3
 */

#ifndef BLISS_H
#define BLISS_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include <stdarg.h>

/* ===================================================================
 * CONSTANTS AND MACROS
 * =================================================================== */

/* Maximum number of refinement rounds to prevent infinite loops */
#ifndef MAX_REFINEMENT_ROUNDS
#define MAX_REFINEMENT_ROUNDS 1000
#endif

/* Initial capacity for adjacency lists */
#ifndef INITIAL_ADJ_CAPACITY  
#define INITIAL_ADJ_CAPACITY 8
#endif

/* Initial capacity for partition cells */
#ifndef INITIAL_CELL_CAPACITY
#define INITIAL_CELL_CAPACITY 64
#endif

/* Hash multiplier for computing graph hashes */
#ifndef HASH_MULTIPLIER
#define HASH_MULTIPLIER 2654435761U
#endif

#ifndef INITIAL_VERTEX_CAPACITY
#define INITIAL_VERTEX_CAPACITY 8
#endif

/* Forward declarations for opaque types */
typedef struct bliss_graph bliss_graph_t;
typedef struct bliss_stats bliss_stats_t;

/* Error codes */
typedef enum {
    BLISS_OK = 0,
    BLISS_ERROR_INVALID_ARGUMENT = -1,
    BLISS_ERROR_OUT_OF_MEMORY = -2,
    BLISS_ERROR_IO = -3,
    BLISS_ERROR_INVALID_GRAPH = -4
} bliss_error_t;

/* Splitting heuristics enumeration */
typedef enum {
    BLISS_SH_F,      /* First non-singleton cell */
    BLISS_SH_FL,     /* First largest non-singleton cell */
    BLISS_SH_FS,     /* First smallest non-singleton cell */
    BLISS_SH_FM,     /* First maximally non-trivially connected non-singleton cell */
    BLISS_SH_FLM,    /* First largest maximally non-trivially connected non-singleton cell */
    BLISS_SH_FSM     /* First smallest maximally non-trivially connected non-singleton cell */
} bliss_splitting_heuristic_t;

/* Type definitions for callbacks */
typedef void (*bliss_automorphism_hook_t)(void *user_param, unsigned int n, const unsigned int *aut);
typedef bool (*bliss_termination_hook_t)(void *user_param);

/* Add these structure definitions to bliss.h */

/* Statistics structure */
struct bliss_stats {
    char *group_size_str;        /* String representation of group size */
    double group_size_approx;    /* Approximation of group size */
    unsigned long nof_nodes;     /* Number of search tree nodes */
    unsigned long nof_leaf_nodes;/* Number of leaf nodes */
    unsigned long nof_bad_nodes; /* Number of bad nodes */
    unsigned long nof_canupdates;/* Number of canonical updates */
    unsigned long nof_generators;/* Number of generators found */
    double max_level;            /* Maximum level in search tree */
};

/* Partition cell structure */
typedef struct partition_cell {
    unsigned int *elements;      
    unsigned int size;           
    unsigned int first;          
    bool is_unit;               
} partition_cell_t;

/* Partition structure */
typedef struct partition {
    partition_cell_t *cells;            
    unsigned int num_cells;             
    unsigned int capacity;              
    unsigned int *element_to_cell;      
    unsigned int *position_in_cell;     
    unsigned int num_discrete_cells;    
} partition_t;

/* Search tree node structure */
typedef struct search_node {
    partition_t partition;              
    unsigned int splitting_vertex;     
    unsigned int target_cell;          
    unsigned int level;                
    struct search_node *parent;        
    struct search_node *child;         
    struct search_node *sibling;       
} search_node_t;

/* Complete search state structure */
typedef struct search_state {
    search_node_t *root;               
    search_node_t *current;            
    unsigned int *best_path;           
    unsigned int path_length;          
    unsigned int max_level;            
    unsigned int **generators;         
    unsigned int num_generators;       
    unsigned int generator_capacity;   
    unsigned int *orbit;               
    unsigned int *orbit_size;          
    bliss_stats_t *stats;              
    bliss_termination_hook_t terminate_func; 
    void *terminate_param;             
    unsigned int *workspace;           
    unsigned int workspace_size;       
} search_state_t;


/* Core graph structure - designed for compiler optimization */
struct bliss_graph {
    /* Basic graph properties */
    unsigned int num_vertices;
    unsigned int num_edges;
    bool is_directed;
    
    /* Vertex colors - dense array for cache efficiency */
    unsigned int *vertex_colors;
    
    /* Adjacency representation - optimized for both sparse and dense graphs */
    unsigned int **adj_lists;     /* Array of adjacency lists */
    unsigned int *adj_list_sizes; /* Sizes of adjacency lists */
    unsigned int *adj_list_caps;  /* Capacities of adjacency lists */
    
    /* For directed graphs - separate in/out adjacency */
    unsigned int **in_adj_lists;     /* Incoming edges for directed graphs */
    unsigned int *in_adj_list_sizes;
    unsigned int *in_adj_list_caps;
    
    /* Search state and optimization data */
    unsigned int *canonical_labeling; /* Current canonical labeling */
    bool canonical_labeling_valid;
    
    /* Configuration options */
    bliss_splitting_heuristic_t splitting_heuristic;
    bool use_component_recursion;
    bool use_failure_recording;
    bool use_long_prune;
    
    /* Internal search structures - kept opaque for optimization */
    void *search_state;
    
    /* Hash value for quick comparisons */
    uint32_t graph_hash;
    bool hash_valid;
    
    /* Memory management */
    unsigned int vertex_capacity;
	
	/* Add this field: */
    search_state_t *internal_search_state;   /* Internal search state */
};

/* Performance monitoring */
typedef struct {
    unsigned long partition_refinements;
    unsigned long cell_splits;
    unsigned long individualization_operations;
    unsigned long orbit_pruning_skips;
    double total_refinement_time;
    double total_individualization_time;
} bliss_performance_stats_t;

/* ===================================================================
 * PARTITION MANIPULATION FUNCTIONS - ADD THESE TO bliss.h
 * =================================================================== */

 /* Forward declarations for functions from bliss_core.c */
void* bliss_malloc(size_t size);
void* bliss_realloc(void* ptr, size_t size);
void bliss_free(void* ptr);
unsigned int min_uint(unsigned int a, unsigned int b);
unsigned int max_uint(unsigned int a, unsigned int b);

/* Partition operations */
partition_t *partition_new(unsigned int num_vertices);
void partition_release(partition_t *partition);
partition_t *partition_copy(const partition_t *original, unsigned int num_vertices);
void partition_add_to_cell(partition_t *partition, unsigned int cell_idx, unsigned int vertex);
bool partition_is_discrete(const partition_t *partition, unsigned int num_vertices);
unsigned int partition_get_cell_of_vertex(const partition_t *partition, unsigned int vertex);
unsigned int partition_get_position_in_cell(const partition_t *partition, unsigned int vertex);

/* Individualization and refinement */
partition_t *individualize_vertex(const partition_t *original_partition,
                                 unsigned int vertex_to_individualize,
                                 unsigned int num_vertices,
                                 const bliss_graph_t *graph);
bool refine_partition_complete(bliss_graph_t *graph, partition_t *partition);
unsigned int select_target_cell(const bliss_graph_t *graph,
                                const partition_t *partition,
                                bliss_splitting_heuristic_t heuristic);

/* Canonical ordering and utilities */
void canonicalize_cell_ordering(partition_t *partition, 
                                unsigned int cell_idx,
                                const bliss_graph_t *graph);
void extract_labeling_from_partition(const partition_t *partition,
                                     unsigned int *labeling,
                                     unsigned int n);

/* Search and pruning */
bool vertices_in_same_orbit(unsigned int v1, unsigned int v2,
                           const search_state_t *state,
                           const bliss_graph_t *graph);
unsigned int get_canonical_representative(const partition_cell_t *cell,
                                         const bliss_graph_t *graph);

/* Validation and debugging */
bool validate_partition(const partition_t *partition, unsigned int num_vertices);
void print_partition(const partition_t *partition, FILE *fp);

/* Main algorithm */
void bliss_find_automorphisms_unified(bliss_graph_t *graph, 
                                     bliss_stats_t *stats,
                                     bliss_automorphism_hook_t hook,
                                     void *hook_user_param);

/* Memory pool management */
void bliss_init_memory_pool(size_t initial_size);
void bliss_reset_memory_pool(void);
void bliss_cleanup_memory_pool(void);
void *bliss_pool_alloc(size_t size);

/* Enhanced automorphism search */
void bliss_find_automorphisms_complete(bliss_graph_t *graph, 
                                      bliss_stats_t *stats,
                                      bliss_automorphism_hook_t hook,
                                      void *hook_user_param);

/* Performance monitoring */
void bliss_reset_performance_stats(void);
const bliss_performance_stats_t *bliss_get_performance_stats(void);
void bliss_print_performance_stats(FILE *fp);

/* Test suite */
bool run_comprehensive_tests(void);

/* ===================================================================
 * CORE API FUNCTIONS
 * =================================================================== */

/* Graph creation and destruction */
bliss_graph_t *bliss_new(unsigned int num_vertices);
bliss_graph_t *bliss_new_directed(unsigned int num_vertices);
void bliss_release(bliss_graph_t *graph);
bliss_graph_t *bliss_copy(const bliss_graph_t *graph);

/* Graph modification */
void bliss_add_vertex(bliss_graph_t *graph, unsigned int color);
void bliss_add_edge(bliss_graph_t *graph, unsigned int v1, unsigned int v2);
void bliss_change_color(bliss_graph_t *graph, unsigned int vertex, unsigned int color);

/* Graph properties and analysis */
unsigned int bliss_get_nof_vertices(const bliss_graph_t *graph);
unsigned int bliss_get_color(const bliss_graph_t *graph, unsigned int vertex);
uint32_t bliss_get_hash(bliss_graph_t *graph);
bool bliss_is_directed(const bliss_graph_t *graph);
unsigned int bliss_get_degree(const bliss_graph_t *graph, unsigned int vertex);
unsigned int bliss_get_in_degree(const bliss_graph_t *graph, unsigned int vertex);
unsigned int bliss_get_out_degree(const bliss_graph_t *graph, unsigned int vertex);
bool bliss_has_edge(const bliss_graph_t *graph, unsigned int v1, unsigned int v2);
const unsigned int *bliss_get_adjacency_list(const bliss_graph_t *graph, 
                                             unsigned int vertex, 
                                             unsigned int *list_size);

/* Graph comparison */
int bliss_cmp(bliss_graph_t *g1, bliss_graph_t *g2);

/* Core algorithms */
void bliss_find_automorphisms(bliss_graph_t *graph, 
                              bliss_stats_t *stats,
                              bliss_automorphism_hook_t hook,
                              void *hook_user_param);

const unsigned int *bliss_find_canonical_labeling(bliss_graph_t *graph,
                                                   bliss_stats_t *stats,
                                                   bliss_automorphism_hook_t hook,
                                                   void *hook_user_param);
bool update_canonical_labeling(search_state_t *state,
                                      const partition_t *current_partition,
                                      unsigned int n);
bool is_canonical_path(const search_state_t *state,
                              const partition_t *current_partition,
                              unsigned int level);
/* Graph transformation */
bliss_graph_t *bliss_permute(const bliss_graph_t *graph, const unsigned int *perm);
bool bliss_is_automorphism(const bliss_graph_t *graph, const unsigned int *perm);

/* Configuration */
void bliss_set_splitting_heuristic(bliss_graph_t *graph, bliss_splitting_heuristic_t heuristic);
void bliss_set_component_recursion(bliss_graph_t *graph, bool enabled);
void bliss_set_failure_recording(bliss_graph_t *graph, bool enabled);
void bliss_set_long_prune(bliss_graph_t *graph, bool enabled);

/* I/O functions */
bliss_graph_t *bliss_read_dimacs(FILE *fp);
void bliss_write_dimacs(const bliss_graph_t *graph, FILE *fp);
void bliss_write_dot(const bliss_graph_t *graph, FILE *fp);
bliss_graph_t *bliss_create_cycle_graph(unsigned int n);
bliss_graph_t *bliss_create_complete_graph(unsigned int n);
bliss_graph_t *bliss_create_path_graph(unsigned int n);
bliss_graph_t *bliss_create_petersen_graph(void);
bliss_graph_t *bliss_create_star_graph(unsigned int n);
bliss_graph_t *bliss_create_random_graph(unsigned int n, double edge_prob, unsigned int seed);
void bliss_print_graph_info(const bliss_graph_t *graph, FILE *fp);
/* File I/O convenience functions */
bliss_graph_t *bliss_read_graph(const char *filename);
bool bliss_write_graph(const bliss_graph_t *graph, const char *filename);
bool bliss_write_dot_file(const bliss_graph_t *graph, const char *filename);

/* Bulk operations */
void bliss_add_edges_bulk(bliss_graph_t *graph, 
                          const unsigned int *edge_list, 
                          unsigned int num_edges);
void bliss_set_colors_bulk(bliss_graph_t *graph, 
                           const unsigned int *colors, 
                           unsigned int num_vertices);

size_t bliss_estimate_memory_usage(const bliss_graph_t *graph);

/* Graph validation */
bool bliss_is_valid_graph(const bliss_graph_t *graph);
bool bliss_is_connected(const bliss_graph_t *graph);
/* Statistics */
bliss_stats_t *bliss_stats_new(void);
void bliss_stats_release(bliss_stats_t *stats);
const char *bliss_stats_get_group_size_string(const bliss_stats_t *stats);
double bliss_stats_get_group_size_approx(const bliss_stats_t *stats);
unsigned long bliss_stats_get_nof_nodes(const bliss_stats_t *stats);
unsigned long bliss_stats_get_nof_leaf_nodes(const bliss_stats_t *stats);
unsigned long bliss_stats_get_nof_generators(const bliss_stats_t *stats);


bliss_error_t bliss_validate_and_repair(bliss_graph_t *graph);
void search_automorphisms_with_individualization(bliss_graph_t* graph,
  search_state_t* state,
  bliss_automorphism_hook_t hook,
  void* hook_param);

search_node_t* create_child_node_unified(const search_node_t* parent,
  unsigned int vertex_to_individualize,
  unsigned int target_cell_idx,
  const bliss_graph_t* graph);

/* Memory and performance optimization hints for compiler */
#ifdef __GNUC__
#define BLISS_LIKELY(x)     __builtin_expect(!!(x), 1)
#define BLISS_UNLIKELY(x)   __builtin_expect(!!(x), 0)
#define BLISS_INLINE        static inline __attribute__((always_inline))
#define BLISS_HOT           __attribute__((hot))
#define BLISS_COLD          __attribute__((cold))
#define BLISS_PURE          __attribute__((pure))
#define BLISS_CONST         __attribute__((const))
#else
#define BLISS_LIKELY(x)     (x)
#define BLISS_UNLIKELY(x)   (x)
#define BLISS_INLINE        static inline
#define BLISS_HOT
#define BLISS_COLD
#define BLISS_PURE
#define BLISS_CONST
#endif


/* Version information */
#define BLISS_VERSION_MAJOR 0
#define BLISS_VERSION_MINOR 77
#define BLISS_VERSION_PATCH 0

const char *bliss_version_string(void);
int bliss_set_verbose_level(int level);
int bliss_set_verbose_file(FILE* fp);
FILE* bliss_get_verbose_file(void);
int bliss_get_verbose_level(void);
void DPRINTF(const char* format, ...);
void DPRINTF_IF(int level, const char* format, ...);

#endif /* BLISS_H */