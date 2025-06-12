/*
 * bliss_core.c - Core implementation of graph automorphism algorithms
 * 
 * This file contains the main algorithmic implementations optimized for
 * interprocedural compiler optimizations and performance.
 */

#include "bliss.h"
#include <time.h>
#include <math.h>



/* Global state for verbose output */
static int verbose_level = 0;
static FILE *verbose_file = NULL;


/* ===================================================================
 * UTILITY FUNCTIONS (INLINED FOR OPTIMIZATION)
 * =================================================================== */

BLISS_INLINE BLISS_HOT 
uint32_t hash_combine(uint32_t h1, uint32_t h2) {
    return h1 ^ (h2 + 0x9e3779b9 + (h1 << 6) + (h1 >> 2));
}

BLISS_INLINE BLISS_HOT
uint32_t hash_array(const unsigned int *arr, unsigned int size) {
    uint32_t hash = 2166136261U;
    for (unsigned int i = 0; i < size; i++) {
        hash = (hash ^ arr[i]) * 16777619U;
    }
    return hash;
}

BLISS_INLINE BLISS_HOT
bool arrays_equal(const unsigned int *a1, const unsigned int *a2, unsigned int size) {
    return memcmp(a1, a2, size * sizeof(unsigned int)) == 0;
}

BLISS_PURE
unsigned int min_uint(unsigned int a, unsigned int b) {
    return a < b ? a : b;
}

BLISS_PURE
unsigned int max_uint(unsigned int a, unsigned int b) {
    return a > b ? a : b;
}

/* ===================================================================
 * MEMORY MANAGEMENT FUNCTIONS
 * =================================================================== */

void *bliss_malloc(size_t size) {
    void *ptr = malloc(size);
    if (BLISS_UNLIKELY(!ptr && size > 0)) {
        fprintf(stderr, "bliss: Out of memory\n");
        abort();
    }
    return ptr;
}

void *bliss_realloc(void *ptr, size_t size) {
    void *new_ptr = realloc(ptr, size);
    if (BLISS_UNLIKELY(!new_ptr && size > 0)) {
        fprintf(stderr, "bliss: Out of memory during realloc\n");
        abort();
    }
    return new_ptr;
}

void bliss_free(void *ptr) {
    if (ptr) free(ptr);
}


/* ===================================================================
 * GRAPH HASH COMPUTATION
 * =================================================================== */

BLISS_HOT
static uint32_t compute_graph_hash(const bliss_graph_t *graph) {
    uint32_t hash = 2166136261U;
    
    /* Hash basic properties */
    hash = hash_combine(hash, graph->num_vertices);
    hash = hash_combine(hash, graph->num_edges);
    hash = hash_combine(hash, graph->is_directed ? 1 : 0);
    
    /* Hash vertex colors */
    for (unsigned int i = 0; i < graph->num_vertices; i++) {
        hash = hash_combine(hash, graph->vertex_colors[i]);
    }
    
    /* Hash adjacency structure */
    for (unsigned int i = 0; i < graph->num_vertices; i++) {
        hash = hash_combine(hash, graph->adj_list_sizes[i]);
        if (graph->adj_list_sizes[i] > 0) {
            hash = hash_combine(hash, hash_array(graph->adj_lists[i], graph->adj_list_sizes[i]));
        }
    }
    
    return hash;
}

/* ===================================================================
 * GRAPH CREATION AND BASIC OPERATIONS
 * =================================================================== */

bliss_graph_t *bliss_new(unsigned int num_vertices) {
    bliss_graph_t *graph = bliss_malloc(sizeof(bliss_graph_t));
    
    graph->num_vertices = num_vertices;
    graph->num_edges = 0;
    graph->is_directed = false;
    graph->vertex_capacity = max_uint(num_vertices, INITIAL_VERTEX_CAPACITY);
    
    /* Allocate vertex arrays */
    graph->vertex_colors = bliss_malloc(graph->vertex_capacity * sizeof(unsigned int));
    graph->adj_lists = bliss_malloc(graph->vertex_capacity * sizeof(unsigned int*));
    graph->adj_list_sizes = bliss_malloc(graph->vertex_capacity * sizeof(unsigned int));
    graph->adj_list_caps = bliss_malloc(graph->vertex_capacity * sizeof(unsigned int));
    
    /* Initialize arrays */
    memset(graph->vertex_colors, 0, graph->vertex_capacity * sizeof(unsigned int));
    memset(graph->adj_lists, 0, graph->vertex_capacity * sizeof(unsigned int*));
    memset(graph->adj_list_sizes, 0, graph->vertex_capacity * sizeof(unsigned int));
    memset(graph->adj_list_caps, 0, graph->vertex_capacity * sizeof(unsigned int));
    
    /* Directed graph arrays */
    graph->in_adj_lists = NULL;
    graph->in_adj_list_sizes = NULL;
    graph->in_adj_list_caps = NULL;
    
    /* Search state */
    graph->canonical_labeling = NULL;
    graph->canonical_labeling_valid = false;
    graph->search_state = NULL;
    
    /* Configuration */
    graph->splitting_heuristic = BLISS_SH_F;
    graph->use_component_recursion = true;
    graph->use_failure_recording = true;
    graph->use_long_prune = true;
    
    /* Hash */
    graph->graph_hash = 0;
    graph->hash_valid = false;
    
    return graph;
}

bliss_graph_t *bliss_new_directed(unsigned int num_vertices) {
    bliss_graph_t *graph = bliss_new(num_vertices);
    graph->is_directed = true;
    
    /* Allocate incoming adjacency arrays for directed graphs */
    graph->in_adj_lists = bliss_malloc(graph->vertex_capacity * sizeof(unsigned int*));
    graph->in_adj_list_sizes = bliss_malloc(graph->vertex_capacity * sizeof(unsigned int));
    graph->in_adj_list_caps = bliss_malloc(graph->vertex_capacity * sizeof(unsigned int));
    
    memset(graph->in_adj_lists, 0, graph->vertex_capacity * sizeof(unsigned int*));
    memset(graph->in_adj_list_sizes, 0, graph->vertex_capacity * sizeof(unsigned int));
    memset(graph->in_adj_list_caps, 0, graph->vertex_capacity * sizeof(unsigned int));
    
    return graph;
}

void bliss_release(bliss_graph_t *graph) {
    if (!graph) return;
    
    /* Free adjacency lists */
    for (unsigned int i = 0; i < graph->vertex_capacity; i++) {
        bliss_free(graph->adj_lists[i]);
        if (graph->is_directed) {
            bliss_free(graph->in_adj_lists[i]);
        }
    }
    
    bliss_free(graph->vertex_colors);
    bliss_free(graph->adj_lists);
    bliss_free(graph->adj_list_sizes);
    bliss_free(graph->adj_list_caps);
    
    if (graph->is_directed) {
        bliss_free(graph->in_adj_lists);
        bliss_free(graph->in_adj_list_sizes);
        bliss_free(graph->in_adj_list_caps);
    }
    
    bliss_free(graph->canonical_labeling);
    bliss_free(graph->search_state);
    bliss_free(graph);
}

/* ===================================================================
 * GRAPH MODIFICATION FUNCTIONS
 * =================================================================== */

BLISS_HOT
void bliss_add_edge(bliss_graph_t *graph, unsigned int v1, unsigned int v2) {
    if (BLISS_UNLIKELY(v1 >= graph->num_vertices || v2 >= graph->num_vertices)) {
        return; /* Invalid vertices */
    }
    
    /* Invalidate cached data */
    graph->hash_valid = false;
    graph->canonical_labeling_valid = false;
    
    /* Add to adjacency list of v1 */
    if (graph->adj_list_sizes[v1] >= graph->adj_list_caps[v1]) {
        graph->adj_list_caps[v1] = max_uint(INITIAL_ADJ_CAPACITY, graph->adj_list_caps[v1] * 2);
        graph->adj_lists[v1] = bliss_realloc(graph->adj_lists[v1], 
                                             graph->adj_list_caps[v1] * sizeof(unsigned int));
    }
    
    graph->adj_lists[v1][graph->adj_list_sizes[v1]++] = v2;
    
    if (graph->is_directed) {
        /* For directed graphs, add to incoming adjacency list of v2 */
        if (graph->in_adj_list_sizes[v2] >= graph->in_adj_list_caps[v2]) {
            graph->in_adj_list_caps[v2] = max_uint(INITIAL_ADJ_CAPACITY, graph->in_adj_list_caps[v2] * 2);
            graph->in_adj_lists[v2] = bliss_realloc(graph->in_adj_lists[v2], 
                                                    graph->in_adj_list_caps[v2] * sizeof(unsigned int));
        }
        graph->in_adj_lists[v2][graph->in_adj_list_sizes[v2]++] = v1;
    } else {
        /* For undirected graphs, add symmetric edge */
        if (v1 != v2) { /* Avoid duplicate self-loops */
            if (graph->adj_list_sizes[v2] >= graph->adj_list_caps[v2]) {
                graph->adj_list_caps[v2] = max_uint(INITIAL_ADJ_CAPACITY, graph->adj_list_caps[v2] * 2);
                graph->adj_lists[v2] = bliss_realloc(graph->adj_lists[v2], 
                                                     graph->adj_list_caps[v2] * sizeof(unsigned int));
            }
            graph->adj_lists[v2][graph->adj_list_sizes[v2]++] = v1;
        }
    }
    
    graph->num_edges++;
}

void bliss_change_color(bliss_graph_t *graph, unsigned int vertex, unsigned int color) {
    if (BLISS_UNLIKELY(vertex >= graph->num_vertices)) {
        return;
    }
    
    graph->vertex_colors[vertex] = color;
    graph->hash_valid = false;
    graph->canonical_labeling_valid = false;
}

/* ===================================================================
 * GRAPH PROPERTY FUNCTIONS
 * =================================================================== */

unsigned int bliss_get_nof_vertices(const bliss_graph_t *graph) {
    return graph ? graph->num_vertices : 0;
}

unsigned int bliss_get_color(const bliss_graph_t *graph, unsigned int vertex) {
    if (BLISS_UNLIKELY(!graph || vertex >= graph->num_vertices)) {
        return 0;
    }
    return graph->vertex_colors[vertex];
}

uint32_t bliss_get_hash(bliss_graph_t *graph) {
    if (BLISS_UNLIKELY(!graph)) {
        return 0;
    }
    
    if (!graph->hash_valid) {
        graph->graph_hash = compute_graph_hash(graph);
        graph->hash_valid = true;
    }
    
    return graph->graph_hash;
}

bool bliss_is_directed(const bliss_graph_t *graph) {
    return graph ? graph->is_directed : false;
}

/* ===================================================================
 * REFINEMENT ALGORITHM - CORE OF BLISS
 * =================================================================== */

BLISS_HOT
static bool refine_partition(bliss_graph_t *graph, partition_t *partition) {
    bool refined = false;
    unsigned int round = 0;
    
    while (round < MAX_REFINEMENT_ROUNDS) {
        bool changed_this_round = false;
        
        /* Process each cell */
        for (unsigned int cell_idx = 0; cell_idx < partition->num_cells; cell_idx++) {
            partition_cell_t *cell = &partition->cells[cell_idx];
            
            if (cell->is_unit) continue; /* Skip unit cells */
            
            /* Count degrees for vertices in this cell */
            unsigned int *degrees = bliss_malloc(cell->size * sizeof(unsigned int));
            memset(degrees, 0, cell->size * sizeof(unsigned int));
            
            /* Compute degrees relative to already processed cells */
            for (unsigned int i = 0; i < cell->size; i++) {
                unsigned int vertex = cell->elements[i];
                
                for (unsigned int j = 0; j < graph->adj_list_sizes[vertex]; j++) {
                    unsigned int neighbor = graph->adj_lists[vertex][j];
                    unsigned int neighbor_cell = partition->element_to_cell[neighbor];
                    
                    /* Only count edges to vertices in smaller-indexed cells */
                    if (neighbor_cell < cell_idx) {
                        degrees[i]++;
                    }
                }
            }
            
            /* Sort vertices by degree (stable sort to maintain canonical order) */
            /* Simple insertion sort for small cells, efficient for this use case */
            for (unsigned int i = 1; i < cell->size; i++) {
                unsigned int key_vertex = cell->elements[i];
                unsigned int key_degree = degrees[i];
                int j = i - 1;
                
                while (j >= 0 && degrees[j] > key_degree) {
                    cell->elements[j + 1] = cell->elements[j];
                    degrees[j + 1] = degrees[j];
                    j--;
                }
                cell->elements[j + 1] = key_vertex;
                degrees[j + 1] = key_degree;
            }
            
            /* Split cell if degrees differ */
            unsigned int split_point = 0;
            for (unsigned int i = 1; i < cell->size; i++) {
                if (degrees[i] != degrees[split_point]) {
                    /* Found a split point */
                    changed_this_round = true;
                    refined = true;
                    /* TODO: Implement actual cell splitting */
                    break;
                }
            }
            
            bliss_free(degrees);
        }
        
        if (!changed_this_round) break;
        round++;
    }
    
    return refined;
}

/* ===================================================================
 * MAIN AUTOMORPHISM SEARCH ALGORITHM
 * =================================================================== */

BLISS_HOT
static void search_automorphisms_recursive(bliss_graph_t *graph, 
                                           search_state_t *state,
                                           bliss_automorphism_hook_t hook,
                                           void *hook_param) {
    /* Check termination condition */
    if (state->terminate_func && state->terminate_func(state->terminate_param)) {
        return;
    }
    
    state->stats->nof_nodes++;
    
    /* Refine current partition */
    bool refined = refine_partition(graph, &state->current->partition);
    (void)refined; /* Mark as used to avoid warning */
    
    /* Check if partition is discrete (canonical labeling found) */
    if (state->current->partition.num_discrete_cells == graph->num_vertices) {
        state->stats->nof_leaf_nodes++;
        
        /* Extract automorphism */
        unsigned int *automorphism = bliss_malloc(graph->num_vertices * sizeof(unsigned int));
        for (unsigned int i = 0; i < graph->num_vertices; i++) {
            automorphism[i] = state->current->partition.cells[i].elements[0];
        }
        
        if (hook) {
            hook(hook_param, graph->num_vertices, automorphism);
        }
        
        /* Store generator if this is a non-trivial automorphism */
        bool is_identity = true;
        for (unsigned int i = 0; i < graph->num_vertices; i++) {
            if (automorphism[i] != i) {
                is_identity = false;
                break;
            }
        }
        
        if (!is_identity) {
            /* Add to generator set */
            if (state->num_generators >= state->generator_capacity) {
                state->generator_capacity *= 2;
                state->generators = bliss_realloc(state->generators, 
                    state->generator_capacity * sizeof(unsigned int*));
            }
            
            state->generators[state->num_generators] = automorphism;
            state->num_generators++;
            state->stats->nof_generators++;
        } else {
            bliss_free(automorphism);
        }
        
        return;
    }
    
    /* Find best cell to split according to splitting heuristic */
    unsigned int best_cell = UINT_MAX;
    unsigned int best_score = 0;
    
    for (unsigned int i = 0; i < state->current->partition.num_cells; i++) {
        partition_cell_t *cell = &state->current->partition.cells[i];
        
        if (cell->is_unit) continue;
        
        unsigned int score = 0;
        switch (graph->splitting_heuristic) {
            case BLISS_SH_F:   /* First non-singleton */
                if (best_cell == UINT_MAX) {
                    best_cell = i;
                }
                break;
                
            case BLISS_SH_FL:  /* First largest */
                score = cell->size;
                if (score > best_score) {
                    best_score = score;
                    best_cell = i;
                }
                break;
                
            case BLISS_SH_FS:  /* First smallest */
                score = UINT_MAX - cell->size;
                if (score > best_score) {
                    best_score = score;
                    best_cell = i;
                }
                break;
                
            case BLISS_SH_FM:  /* First maximally connected */
            case BLISS_SH_FLM: /* First largest maximally connected */
            case BLISS_SH_FSM: /* First smallest maximally connected */
                /* Compute connectivity score */
                for (unsigned int j = 0; j < cell->size; j++) {
                    unsigned int vertex = cell->elements[j];
                    score += graph->adj_list_sizes[vertex];
                }
                
                if (graph->splitting_heuristic == BLISS_SH_FLM) {
                    score = score * cell->size;
                } else if (graph->splitting_heuristic == BLISS_SH_FSM) {
                    score = score * (UINT_MAX - cell->size);
                }
                
                if (score > best_score) {
                    best_score = score;
                    best_cell = i;
                }
                break;
        }
    }
    
    if (best_cell == UINT_MAX) {
        /* All cells are unit - should not happen here */
        return;
    }
    
    /* Branch on each vertex in the selected cell */
    partition_cell_t *split_cell = &state->current->partition.cells[best_cell];
    
    for (unsigned int i = 0; i < split_cell->size; i++) {
        unsigned int split_vertex = split_cell->elements[i];
        
        /* Create child node */
        search_node_t *child = bliss_malloc(sizeof(search_node_t));
        child->parent = state->current;
        child->child = NULL;
        child->sibling = NULL;
        child->splitting_vertex = split_vertex;
        child->target_cell = best_cell;
        child->level = state->current->level + 1;
        #if BLISS_DEBUG&2
        debug_partition_state(&parent->partition, "BEFORE child creation");  
        #endif
        /* Initialize child partition */
        /* Copy partition and make split_vertex first in its cell */
        child->partition = *partition_new(graph->num_vertices);
        #if BLISS_DEBUG&2
        debug_partition_state(&child->partition, "AFTER child creation");
        #endif
        /* TODO: Implement partition copying and vertex movement */
        
        /* Link child to parent */
        if (!state->current->child) {
            state->current->child = child;
        } else {
            search_node_t *sibling = state->current->child;
            while (sibling->sibling) {
                sibling = sibling->sibling;
            }
            sibling->sibling = child;
        }
        
        /* Recurse on child */
        search_node_t *prev_current = state->current;
        state->current = child;
        search_automorphisms_recursive(graph, state, hook, hook_param);
        state->current = prev_current;
    }
}

/* ===================================================================
 * PUBLIC API IMPLEMENTATION
 * =================================================================== */

void bliss_find_automorphisms_incomplete(bliss_graph_t *graph, 
                              bliss_stats_t *stats,
                              bliss_automorphism_hook_t hook,
                              void *hook_user_param) {
    if (BLISS_UNLIKELY(!graph || !stats)) {
        return;
    }
    
    /* Initialize search state */
    search_state_t *state = bliss_malloc(sizeof(search_state_t));
    state->stats = stats;
    state->terminate_func = NULL;
    state->terminate_param = NULL;
    
    /* Initialize generator storage */
    state->generator_capacity = 16;
    state->generators = bliss_malloc(state->generator_capacity * sizeof(unsigned int*));
    state->num_generators = 0;
    
    /* Create initial partition based on vertex colors */
    state->root = bliss_malloc(sizeof(search_node_t));
    state->root->parent = NULL;
    state->root->child = NULL;
    state->root->sibling = NULL;
    state->root->level = 0;
    state->root->partition = *partition_new(graph->num_vertices);
    
    /* Group vertices by color */
    unsigned int max_color = 0;
    for (unsigned int i = 0; i < graph->num_vertices; i++) {
        if (graph->vertex_colors[i] > max_color) {
            max_color = graph->vertex_colors[i];
        }
    }
    
    /* Create color-based partition */
    for (unsigned int color = 0; color <= max_color; color++) {
        for (unsigned int vertex = 0; vertex < graph->num_vertices; vertex++) {
            if (graph->vertex_colors[vertex] == color) {
                /* Find or create cell for this color */
                unsigned int cell_idx = color;
                if (cell_idx >= state->root->partition.num_cells) {
                    state->root->partition.num_cells = cell_idx + 1;
                }
                partition_add_to_cell(&state->root->partition, cell_idx, vertex);
            }
        }
    }
    
    /* Initialize statistics */
    stats->nof_nodes = 0;
    stats->nof_leaf_nodes = 0;
    stats->nof_bad_nodes = 0;
    stats->nof_canupdates = 0;
    stats->nof_generators = 0;
    stats->max_level = 0;
    stats->group_size_str = NULL;
    stats->group_size_approx = 1.0;
    
    /* Start search */
    state->current = state->root;
    search_automorphisms_recursive(graph, state, hook, hook_user_param);
    
    /* Compute group size approximation */
    stats->group_size_approx = 1.0;
    for (unsigned int i = 0; i < state->num_generators; i++) {
        stats->group_size_approx *= 2.0; /* Rough approximation */
    }
    
    /* Clean up generators */
    for (unsigned int i = 0; i < state->num_generators; i++) {
        bliss_free(state->generators[i]);
    }
    bliss_free(state->generators);
    
    partition_release(&state->root->partition);
    bliss_free(state->root);
    bliss_free(state);
}

/* Enhanced search with proper individualization and orbit pruning */
void search_automorphisms_with_individualization(bliss_graph_t *graph,
                                                       search_state_t *state,
                                                       bliss_automorphism_hook_t hook,
                                                       void *hook_param) {
    /* Check termination condition */
    if (state->terminate_func && state->terminate_func(state->terminate_param)) {
        return;
    }

    state->stats->nof_nodes++;

    /* Refine current partition to equitable form */
    refine_partition_complete(graph, &state->current->partition);

    /* Check if we've reached a discrete partition */
    if (state->current->partition.num_discrete_cells == graph->num_vertices) {
        state->stats->nof_leaf_nodes++;

        /* Extract automorphism and process it */
        unsigned int *automorphism = bliss_malloc(graph->num_vertices * sizeof(unsigned int));
        extract_labeling_from_partition(&state->current->partition, automorphism,
                                       graph->num_vertices);

        if (hook) {
            hook(hook_param, graph->num_vertices, automorphism);
        }

        /* Store generator if not identity */
        bool is_identity = true;
        for (unsigned int i = 0; i < graph->num_vertices; i++) {
            if (automorphism[i] != i) {
                is_identity = false;
                break;
            }
        }

        if (!is_identity) {
            if (state->num_generators >= state->generator_capacity) {
                state->generator_capacity *= 2;
                state->generators = bliss_realloc(state->generators,
                    state->generator_capacity * sizeof(unsigned int*));
            }

            state->generators[state->num_generators] = automorphism;
            state->num_generators++;
            state->stats->nof_generators++;
        } else {
            bliss_free(automorphism);
        }

        return;
    }

    /* Select target cell for branching */
    unsigned int target_cell = select_target_cell(graph, &state->current->partition,
                                                 graph->splitting_heuristic);

    if (target_cell == UINT_MAX) {
        return; /* No suitable cell found */
    }

    /* Get canonical representative for orbit pruning */
    partition_cell_t *split_cell = &state->current->partition.cells[target_cell];
    unsigned int canonical_rep = get_canonical_representative(split_cell, graph);
    if (canonical_rep == UINT_MAX) {
      return; /* No canonical representative found */
    }

    /* Branch on each vertex in the selected cell */
    for (unsigned int i = 0; i < split_cell->size; i++) {
        unsigned int split_vertex = split_cell->elements[i];

        /* Orbit pruning: skip vertices in the same orbit as already processed ones */
        bool skip_vertex = false;
        for (unsigned int j = 0; j < i; j++) {
            if (vertices_in_same_orbit(split_vertex, split_cell->elements[j], state, graph)) {
                skip_vertex = true;
                break;
            }
        }

        if (skip_vertex) {
            continue;
        }

        /* Create child node with individualized partition */
        search_node_t *child = create_child_node_unified(state->current, split_vertex,
                                                 target_cell, graph);

        if (!child) {
            continue; /* Skip if individualization failed */
        }

        /* Update maximum level reached */
        if (child->level > state->stats->max_level) {
            state->stats->max_level = child->level;
        }

        /* Recursively search the child */
        search_node_t *prev_current = state->current;
        state->current = child;
        search_automorphisms_with_individualization(graph, state, hook, hook_param);
        state->current = prev_current;

        /* Clean up child node */
        partition_release(&child->partition);
        bliss_free(child);
    }
}

/* Updated main function that uses proper individualization */
void bliss_find_automorphisms_complete(bliss_graph_t *graph,
                                      bliss_stats_t *stats,
                                      bliss_automorphism_hook_t hook,
                                      void *hook_user_param) {
    if (BLISS_UNLIKELY(!graph || !stats)) {
        return;
    }

    /* Initialize search state */
    search_state_t *state = bliss_malloc(sizeof(search_state_t));
    state->stats = stats;
    state->terminate_func = NULL;
    state->terminate_param = NULL;
    state->best_path = NULL;
    state->path_length = 0;

    /* Initialize generator storage */
    state->generator_capacity = 16;
    state->generators = bliss_malloc(state->generator_capacity * sizeof(unsigned int*));
    state->num_generators = 0;

    /* Create initial partition based on vertex colors */
    state->root = bliss_malloc(sizeof(search_node_t));
    state->root->parent = NULL;
    state->root->child = NULL;
    state->root->sibling = NULL;
    state->root->level = 0;
    state->root->partition = *partition_new(graph->num_vertices);

    /* Build initial color-based partition */
    unsigned int max_color = 0;
    for (unsigned int i = 0; i < graph->num_vertices; i++) {
        if (graph->vertex_colors[i] > max_color) {
            max_color = graph->vertex_colors[i];
        }
    }

    /* Create partition cells for each color */
    for (unsigned int color = 0; color <= max_color; color++) {
        for (unsigned int vertex = 0; vertex < graph->num_vertices; vertex++) {
            if (graph->vertex_colors[vertex] == color) {
                if (color >= state->root->partition.num_cells) {
                    state->root->partition.num_cells = color + 1;
                }
                partition_add_to_cell(&state->root->partition, color, vertex);
            }
        }
    }

    /* Canonicalize initial partition */
    for (unsigned int i = 0; i < state->root->partition.num_cells; i++) {
        canonicalize_cell_ordering(&state->root->partition, i, graph);
    }

    /* Initialize statistics */
    stats->nof_nodes = 0;
    stats->nof_leaf_nodes = 0;
    stats->nof_bad_nodes = 0;
    stats->nof_canupdates = 0;
    stats->nof_generators = 0;
    stats->max_level = 0;
    stats->group_size_str = NULL;
    stats->group_size_approx = 1.0;

    /* Start search with proper individualization */
    state->current = state->root;
    search_automorphisms_with_individualization(graph, state, hook, hook_user_param);

    /* Compute group size approximation */
    stats->group_size_approx = 1.0;
    for (unsigned int i = 0; i < state->num_generators; i++) {
        stats->group_size_approx *= 2.0;
    }

    /* Clean up */
    for (unsigned int i = 0; i < state->num_generators; i++) {
        bliss_free(state->generators[i]);
    }
    bliss_free(state->generators);
    bliss_free(state->best_path);

    partition_release(&state->root->partition);
    bliss_free(state->root);
    bliss_free(state);
}

/* Improved automorphism search with proper pruning */
BLISS_HOT
static void search_automorphisms_improved(bliss_graph_t* graph,
  search_state_t* state,
  bliss_automorphism_hook_t hook,
  void* hook_param) {
  /* Check termination condition */
  if (state->terminate_func && state->terminate_func(state->terminate_param)) {
    return;
  }

  state->stats->nof_nodes++;

  /* Refine current partition to equitable form */
  bool refined = refine_partition_complete(graph, &state->current->partition);
  (void)refined; /* Avoid unused variable warning */

  /* Check if we've reached a discrete partition (all cells are units) */
  if (state->current->partition.num_discrete_cells == graph->num_vertices) {
    state->stats->nof_leaf_nodes++;

    /* Update canonical labeling if this is the best so far */
    bool is_canonical = update_canonical_labeling(state, &state->current->partition,
      graph->num_vertices);

    if (is_canonical) {
      /* This is the new canonical form */
      state->stats->nof_canupdates++;
    }
    else {
      /* Found an automorphism */
      unsigned int* automorphism = bliss_malloc(graph->num_vertices * sizeof(unsigned int));
      extract_labeling_from_partition(&state->current->partition, automorphism,
        graph->num_vertices);

      if (hook) {
        hook(hook_param, graph->num_vertices, automorphism);
      }

      /* Store generator */
      if (state->num_generators >= state->generator_capacity) {
        state->generator_capacity *= 2;
        state->generators = bliss_realloc(state->generators,
          state->generator_capacity * sizeof(unsigned int*));
      }

      state->generators[state->num_generators] = automorphism;
      state->num_generators++;
      state->stats->nof_generators++;
    }

    return;
  }

  /* Canonical path pruning */
  if (!is_canonical_path(state, &state->current->partition, state->current->level)) {
    state->stats->nof_bad_nodes++;
    return;
  }

  /* Select target cell for branching */
  unsigned int target_cell = select_target_cell(graph, &state->current->partition,
    graph->splitting_heuristic);

  if (target_cell == UINT_MAX) {
    /* All cells are units - should not happen here */
    return;
  }

  /* Branch on each vertex in the selected cell */
  partition_cell_t* split_cell = &state->current->partition.cells[target_cell];

  for (unsigned int i = 0; i < split_cell->size; i++) {
    unsigned int split_vertex = split_cell->elements[i];

    /* Create child node - simplified for this example */
    search_node_t* child = bliss_malloc(sizeof(search_node_t));
    child->parent = state->current;
    child->child = NULL;
    child->sibling = NULL;
    child->splitting_vertex = split_vertex;
    child->target_cell = target_cell;
    child->level = state->current->level + 1;

    /* Copy partition and individualize the selected vertex */
    child->partition = *partition_new(graph->num_vertices);

    /* TODO: Implement proper partition copying and vertex individualization */
    /* This is a critical missing piece that requires careful implementation */

    /* Update maximum level reached */
    if (child->level > state->stats->max_level) {
      state->stats->max_level = child->level;
    }

    /* Recursively search the child */
    search_node_t* prev_current = state->current;
    state->current = child;
    search_automorphisms_improved(graph, state, hook, hook_param);
    state->current = prev_current;

    /* Clean up child node */
    partition_release(&child->partition);
    bliss_free(child);
  }
}

/* ===================================================================
 * REPLACEMENT FOR MAIN AUTOMORPHISM FUNCTION
 * =================================================================== */

 /* Replace the existing bliss_find_automorphisms with this improved version */
void bliss_find_automorphisms_improved(bliss_graph_t* graph,
  bliss_stats_t* stats,
  bliss_automorphism_hook_t hook,
  void* hook_user_param) {
  if (BLISS_UNLIKELY(!graph || !stats)) {
    return;
  }

  /* Initialize search state */
  search_state_t* state = bliss_malloc(sizeof(search_state_t));
  state->stats = stats;
  state->terminate_func = NULL;
  state->terminate_param = NULL;
  state->best_path = NULL;
  state->path_length = 0;

  /* Initialize generator storage */
  state->generator_capacity = 16;
  state->generators = bliss_malloc(state->generator_capacity * sizeof(unsigned int*));
  state->num_generators = 0;

  /* Create initial partition based on vertex colors */
  state->root = bliss_malloc(sizeof(search_node_t));
  state->root->parent = NULL;
  state->root->child = NULL;
  state->root->sibling = NULL;
  state->root->level = 0;
  state->root->partition = *partition_new(graph->num_vertices);

  /* Build initial color-based partition */
  unsigned int max_color = 0;
  for (unsigned int i = 0; i < graph->num_vertices; i++) {
    if (graph->vertex_colors[i] > max_color) {
      max_color = graph->vertex_colors[i];
    }
  }

  /* Create partition cells for each color */
  for (unsigned int color = 0; color <= max_color; color++) {
    bool color_found = false;
    for (unsigned int vertex = 0; vertex < graph->num_vertices; vertex++) {
      if (graph->vertex_colors[vertex] == color) {
        if (!color_found) {
          /* Create new cell for this color */
          if (state->root->partition.num_cells <= color) {
            state->root->partition.num_cells = color + 1;
          }
          color_found = true;
        }
        partition_add_to_cell(&state->root->partition, color, vertex);
      }
    }
  }

  /* Initialize statistics */
  stats->nof_nodes = 0;
  stats->nof_leaf_nodes = 0;
  stats->nof_bad_nodes = 0;
  stats->nof_canupdates = 0;
  stats->nof_generators = 0;
  stats->max_level = 0;
  stats->group_size_str = NULL;
  stats->group_size_approx = 1.0;

  /* Start improved search */
  state->current = state->root;
  search_automorphisms_improved(graph, state, hook, hook_user_param);

  /* Compute better group size approximation */
  /* This is still simplified - proper implementation would use Schreier-Sims */
  stats->group_size_approx = 1.0;
  for (unsigned int i = 0; i < state->num_generators; i++) {
    stats->group_size_approx *= 2.0; /* Very rough approximation */
  }

  /* Clean up */
  for (unsigned int i = 0; i < state->num_generators; i++) {
    bliss_free(state->generators[i]);
  }
  bliss_free(state->generators);
  bliss_free(state->best_path);

  partition_release(&state->root->partition);
  bliss_free(state->root);
  bliss_free(state);
}

/* ===================================================================
 * GRAPH COMPARISON AND PERMUTATION
 * =================================================================== */

int bliss_cmp(bliss_graph_t *g1, bliss_graph_t *g2) {
    if (!g1 && !g2) return 0;
    if (!g1) return -1;
    if (!g2) return 1;
    
    /* Compare basic properties */
    if (g1->num_vertices != g2->num_vertices) {
        return g1->num_vertices < g2->num_vertices ? -1 : 1;
    }
    
    if (g1->num_edges != g2->num_edges) {
        return g1->num_edges < g2->num_edges ? -1 : 1;
    }
    
    if (g1->is_directed != g2->is_directed) {
        return g1->is_directed ? 1 : -1;
    }
    
    /* Compare hashes first for quick rejection */
    uint32_t h1 = bliss_get_hash(g1);
    uint32_t h2 = bliss_get_hash(g2);
    if (h1 != h2) {
        return h1 < h2 ? -1 : 1;
    }
    
    /* Compare vertex colors */
    for (unsigned int i = 0; i < g1->num_vertices; i++) {
        if (g1->vertex_colors[i] != g2->vertex_colors[i]) {
            return g1->vertex_colors[i] < g2->vertex_colors[i] ? -1 : 1;
        }
    }
    
    /* Compare adjacency structure */
    for (unsigned int i = 0; i < g1->num_vertices; i++) {
        if (g1->adj_list_sizes[i] != g2->adj_list_sizes[i]) {
            return g1->adj_list_sizes[i] < g2->adj_list_sizes[i] ? -1 : 1;
        }
        
        for (unsigned int j = 0; j < g1->adj_list_sizes[i]; j++) {
            if (g1->adj_lists[i][j] != g2->adj_lists[i][j]) {
                return g1->adj_lists[i][j] < g2->adj_lists[i][j] ? -1 : 1;
            }
        }
    }
    
    return 0; /* Graphs are equal */
}

bliss_graph_t *bliss_permute(const bliss_graph_t *graph, const unsigned int *perm) {
    if (BLISS_UNLIKELY(!graph || !perm)) {
        return NULL;
    }
    
    bliss_graph_t *result = graph->is_directed ? 
        bliss_new_directed(graph->num_vertices) : 
        bliss_new(graph->num_vertices);
    
    /* Apply permutation to vertex colors */
    for (unsigned int i = 0; i < graph->num_vertices; i++) {
        result->vertex_colors[perm[i]] = graph->vertex_colors[i];
    }
    
    /* Apply permutation to edges */
    for (unsigned int i = 0; i < graph->num_vertices; i++) {
        for (unsigned int j = 0; j < graph->adj_list_sizes[i]; j++) {
            unsigned int neighbor = graph->adj_lists[i][j];
            bliss_add_edge(result, perm[i], perm[neighbor]);
        }
    }
    
    return result;
}

bool bliss_is_automorphism(const bliss_graph_t *graph, const unsigned int *perm) {
    if (BLISS_UNLIKELY(!graph || !perm)) {
        return false;
    }
    
    /* Check if permutation preserves vertex colors */
    for (unsigned int i = 0; i < graph->num_vertices; i++) {
        if (graph->vertex_colors[i] != graph->vertex_colors[perm[i]]) {
            return false;
        }
    }
    
    /* Check if permutation preserves edges */
    for (unsigned int i = 0; i < graph->num_vertices; i++) {
        for (unsigned int j = 0; j < graph->adj_list_sizes[i]; j++) {
            unsigned int neighbor = graph->adj_lists[i][j];
            unsigned int perm_i = perm[i];
            unsigned int perm_neighbor = perm[neighbor];
            
            /* Check if edge (perm_i, perm_neighbor) exists */
            bool found = false;
            for (unsigned int k = 0; k < graph->adj_list_sizes[perm_i]; k++) {
                if (graph->adj_lists[perm_i][k] == perm_neighbor) {
                    found = true;
                    break;
                }
            }
            
            if (!found) {
                return false;
            }
        }
    }
    
    return true;
}

/* ===================================================================
 * CONFIGURATION FUNCTIONS
 * =================================================================== */

void bliss_set_splitting_heuristic(bliss_graph_t *graph, bliss_splitting_heuristic_t heuristic) {
    if (graph) {
        graph->splitting_heuristic = heuristic;
        graph->canonical_labeling_valid = false;
    }
}

void bliss_set_component_recursion(bliss_graph_t *graph, bool enabled) {
    if (graph) {
        graph->use_component_recursion = enabled;
        graph->canonical_labeling_valid = false;
    }
}

void bliss_set_failure_recording(bliss_graph_t *graph, bool enabled) {
    if (graph) {
        graph->use_failure_recording = enabled;
    }
}

void bliss_set_long_prune(bliss_graph_t *graph, bool enabled) {
    if (graph) {
        graph->use_long_prune = enabled;
        graph->canonical_labeling_valid = false;
    }
}

/* ===================================================================
 * STATISTICS FUNCTIONS
 * =================================================================== */

bliss_stats_t *bliss_stats_new(void) {
    bliss_stats_t *stats = bliss_malloc(sizeof(bliss_stats_t));
    memset(stats, 0, sizeof(bliss_stats_t));
    return stats;
}

void bliss_stats_release(bliss_stats_t *stats) {
    if (stats) {
        bliss_free(stats->group_size_str);
        bliss_free(stats);
    }
}

const char *bliss_stats_get_group_size_string(const bliss_stats_t *stats) {
    return stats ? stats->group_size_str : NULL;
}

double bliss_stats_get_group_size_approx(const bliss_stats_t *stats) {
    return stats ? stats->group_size_approx : 0.0;
}

unsigned long bliss_stats_get_nof_nodes(const bliss_stats_t *stats) {
    return stats ? stats->nof_nodes : 0;
}

unsigned long bliss_stats_get_nof_leaf_nodes(const bliss_stats_t *stats) {
    return stats ? stats->nof_leaf_nodes : 0;
}

unsigned long bliss_stats_get_nof_generators(const bliss_stats_t *stats) {
    return stats ? stats->nof_generators : 0;
}

/* ===================================================================
 * UTILITY FUNCTIONS
 * =================================================================== */

void bliss_set_verbose_level(int level) {
    verbose_level = level;
}

void bliss_set_verbose_file(FILE *fp) {
    verbose_file = fp;
}

const char *bliss_version_string(void) {
    return "0.77.0-c";
}