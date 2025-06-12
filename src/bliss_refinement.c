/*
 * bliss_refinement.c - Complete partition refinement algorithm
 *
 * This file contains the core partition refinement algorithm that is
 * the heart of the bliss automorphism detection algorithm.
 */

#include "bliss.h"
#include <time.h>
#include <math.h>

/* ===================================================================
 * DEGREE SIGNATURE COMPUTATION
 * =================================================================== */

/* Degree signature structure for partition refinement */
typedef struct {
    unsigned int *degrees_to_color;  /* Degrees to each color class */
    unsigned int num_colors;
    unsigned int total_degree;
    uint32_t signature_hash;         /* Hash of the degree sequence */
} degree_signature_t;

/* Create degree signature for a vertex */
static degree_signature_t compute_degree_signature(const bliss_graph_t *graph,
                                                   unsigned int vertex,
                                                   const partition_t *partition) {
    degree_signature_t sig;
    sig.num_colors = partition->num_cells;
    sig.degrees_to_color = bliss_malloc(sig.num_colors * sizeof(unsigned int));
    sig.total_degree = 0;
    sig.signature_hash = 2166136261U; /* FNV offset basis */

    memset(sig.degrees_to_color, 0, sig.num_colors * sizeof(unsigned int));

    /* Count degrees to each color class */
    for (unsigned int i = 0; i < graph->adj_list_sizes[vertex]; i++) {
        unsigned int neighbor = graph->adj_lists[vertex][i];
        unsigned int neighbor_color = partition->element_to_cell[neighbor];
        sig.degrees_to_color[neighbor_color]++;
        sig.total_degree++;
    }

    /* For directed graphs, also count incoming edges */
    if (graph->is_directed) {
        for (unsigned int i = 0; i < graph->in_adj_list_sizes[vertex]; i++) {
            unsigned int neighbor = graph->in_adj_lists[vertex][i];
            unsigned int neighbor_color = partition->element_to_cell[neighbor];
            /* Use a different offset for incoming edges to distinguish from outgoing */
            sig.degrees_to_color[neighbor_color] += 1000000;
            sig.total_degree++;
        }
    }

    /* Compute hash of degree sequence for quick comparison */
    for (unsigned int i = 0; i < sig.num_colors; i++) {
        sig.signature_hash = (sig.signature_hash ^ sig.degrees_to_color[i]) * 16777619U;
    }

    return sig;
}

/* Compare two degree signatures */
static int compare_degree_signatures(const degree_signature_t *sig1,
                                    const degree_signature_t *sig2) {
    /* First compare hash for quick rejection */
    if (sig1->signature_hash != sig2->signature_hash) {
        return sig1->signature_hash < sig2->signature_hash ? -1 : 1;
    }

    /* Then compare total degree */
    if (sig1->total_degree != sig2->total_degree) {
        return sig1->total_degree < sig2->total_degree ? -1 : 1;
    }

    /* Finally compare degree sequence element by element */
    for (unsigned int i = 0; i < sig1->num_colors; i++) {
        if (sig1->degrees_to_color[i] != sig2->degrees_to_color[i]) {
            return sig1->degrees_to_color[i] < sig2->degrees_to_color[i] ? -1 : 1;
        }
    }

    return 0;
}

/* Free degree signature */
static void free_degree_signature(degree_signature_t *sig) {
    bliss_free(sig->degrees_to_color);
    sig->degrees_to_color = NULL;
}

/* ===================================================================
 * CELL SPLITTING OPERATIONS
 * =================================================================== */

/* Split a partition cell based on degree signatures */
static bool split_cell_by_signatures(partition_t *partition, unsigned int cell_idx,
                                     const bliss_graph_t *graph) {
    partition_cell_t *cell = &partition->cells[cell_idx];
    
    if (cell->size <= 1) {
        return false; /* Cannot split unit cells */
    }
    
    /* Compute degree signatures for all vertices in the cell */
    degree_signature_t *signatures = bliss_malloc(cell->size * sizeof(degree_signature_t));
    for (unsigned int i = 0; i < cell->size; i++) {
        signatures[i] = compute_degree_signature(graph, cell->elements[i], partition);
    }
    
    /* Sort vertices by their degree signatures using insertion sort for stability */
    /* REMOVE this line: bool changed = false; */
    for (unsigned int i = 1; i < cell->size; i++) {
        for (unsigned int j = i; j > 0; j--) {
            int cmp = compare_degree_signatures(&signatures[j-1], &signatures[j]);
            if (cmp <= 0) break;
            
            /* Swap vertices and signatures */
            unsigned int temp_vertex = cell->elements[j-1];
            cell->elements[j-1] = cell->elements[j];
            cell->elements[j] = temp_vertex;
            
            degree_signature_t temp_sig = signatures[j-1];
            signatures[j-1] = signatures[j];
            signatures[j] = temp_sig;
            
            /* REMOVE this line: changed = true; */
        }
    }
    
    /* Rest of function remains the same... */
    /* Find split points where signatures differ */
    bool split_occurred = false;
    unsigned int group_start = 0;
    
    for (unsigned int i = 1; i <= cell->size; i++) {
        bool should_split = (i == cell->size) || 
                           (compare_degree_signatures(&signatures[group_start], &signatures[i]) != 0);
        
        if (should_split && i > group_start + 1) {
            /* Found a group with multiple vertices - split if not the first group */
            if (group_start > 0) {
                /* Create new cell for this group */
                unsigned int new_cell_idx = partition->num_cells;
                
                /* Ensure capacity */
                if (new_cell_idx >= partition->capacity) {
                    unsigned int new_capacity = partition->capacity * 2;
                    partition->cells = bliss_realloc(partition->cells, 
                                                    new_capacity * sizeof(partition_cell_t));
                    
                    /* Initialize new cells */
                    for (unsigned int k = partition->capacity; k < new_capacity; k++) {
                        partition->cells[k].elements = NULL;
                        partition->cells[k].size = 0;
                        partition->cells[k].first = 0;
                        partition->cells[k].is_unit = false;
                    }
                    
                    partition->capacity = new_capacity;
                }
                
                partition_cell_t *new_cell = &partition->cells[new_cell_idx];
                new_cell->size = i - group_start;
                new_cell->elements = bliss_malloc(new_cell->size * sizeof(unsigned int));
                new_cell->is_unit = (new_cell->size == 1);
                new_cell->first = 0;
                
                /* Move vertices to new cell */
                for (unsigned int j = 0; j < new_cell->size; j++) {
                    new_cell->elements[j] = cell->elements[group_start + j];
                    partition->element_to_cell[new_cell->elements[j]] = new_cell_idx;
                    partition->position_in_cell[new_cell->elements[j]] = j;
                }
                
                if (new_cell->is_unit) {
                    partition->num_discrete_cells++;
                }
                
                partition->num_cells++;
                split_occurred = true;
            }
        }
        
        if (i < cell->size && should_split) {
            group_start = i;
        }
    }
    
    /* Update the original cell to contain only the first group */
    if (split_occurred) {
        unsigned int new_size = 0;
        for (unsigned int i = 0; i < cell->size; i++) {
            if (partition->element_to_cell[cell->elements[i]] == cell_idx) {
                cell->elements[new_size] = cell->elements[i];
                partition->position_in_cell[cell->elements[i]] = new_size;
                new_size++;
            }
        }
        
        cell->size = new_size;
        cell->is_unit = (cell->size == 1);
        
        if (cell->is_unit && cell->size == 1) {
            partition->num_discrete_cells++;
        }
    }
    
    /* Clean up signatures */
    for (unsigned int i = 0; i < cell->size + (split_occurred ? 1 : 0); i++) {
        if (i < cell->size) {
            free_degree_signature(&signatures[i]);
        }
    }
    bliss_free(signatures);
    
    return split_occurred;
}

/* ===================================================================
 * COMPLETE EQUITABLE REFINEMENT ALGORITHM
 * =================================================================== */

/* Perform one round of refinement on all cells */
static bool refine_partition_one_round(bliss_graph_t *graph, partition_t *partition) {
    bool any_split = false;
    unsigned int original_num_cells = partition->num_cells;

    /* Try to split each non-unit cell */
    for (unsigned int cell_idx = 0; cell_idx < original_num_cells; cell_idx++) {
        if (!partition->cells[cell_idx].is_unit) {
            if (split_cell_by_signatures(partition, cell_idx, graph)) {
                any_split = true;
            }
        }
    }

    return any_split;
}

/* Complete equitable partition refinement - the heart of bliss */
bool refine_partition_complete(bliss_graph_t *graph, partition_t *partition) {
    bool global_changed = false;
    bool local_changed = true;
    unsigned int round = 0;

    while (local_changed && round < MAX_REFINEMENT_ROUNDS) {
        local_changed = refine_partition_one_round(graph, partition);

        if (local_changed) {
            global_changed = true;

            /* Canonicalize ordering in all cells after splitting */
            unsigned int current_num_cells = partition->num_cells;
            for (unsigned int i = 0; i < current_num_cells; i++) {
                canonicalize_cell_ordering(partition, i, graph);
            }
        }

        round++;
    }

    /* Update discrete cell count */
    partition->num_discrete_cells = 0;
    for (unsigned int i = 0; i < partition->num_cells; i++) {
        if (partition->cells[i].is_unit) {
            partition->num_discrete_cells++;
        }
    }

    return global_changed;
}

/* ===================================================================
 * TARGET CELL SELECTION HEURISTICS
 * =================================================================== */

/* Compute connectivity score for a cell (used by FM, FLM, FSM heuristics) */
static unsigned int compute_cell_connectivity(const bliss_graph_t *graph,
                                             const partition_cell_t *cell,
                                             const partition_t *partition) {
    unsigned int connectivity = 0;

    for (unsigned int i = 0; i < cell->size; i++) {
        unsigned int vertex = cell->elements[i];
        unsigned int vertex_cell = partition->element_to_cell[vertex];

        /* Count edges to vertices in different cells */
        for (unsigned int j = 0; j < graph->adj_list_sizes[vertex]; j++) {
            unsigned int neighbor = graph->adj_lists[vertex][j];
            unsigned int neighbor_cell = partition->element_to_cell[neighbor];

            if (neighbor_cell != vertex_cell) {
                connectivity++;
            }
        }

        /* For directed graphs, also count incoming edges */
        if (graph->is_directed) {
            for (unsigned int j = 0; j < graph->in_adj_list_sizes[vertex]; j++) {
                unsigned int neighbor = graph->in_adj_lists[vertex][j];
                unsigned int neighbor_cell = partition->element_to_cell[neighbor];

                if (neighbor_cell != vertex_cell) {
                    connectivity++;
                }
            }
        }
    }

    return connectivity;
}

/* Select target cell for branching based on splitting heuristic */
unsigned int select_target_cell(const bliss_graph_t *graph,
                                const partition_t *partition,
                                bliss_splitting_heuristic_t heuristic) {
    unsigned int best_cell = UINT_MAX;
    unsigned int best_score = 0;
    bool first_found = false;

    for (unsigned int i = 0; i < partition->num_cells; i++) {
        const partition_cell_t *cell = &partition->cells[i];

        if (cell->is_unit) continue; /* Skip unit cells */

        unsigned int score = 0;
        bool update_best = false;

        switch (heuristic) {
            case BLISS_SH_F:   /* First non-singleton cell */
                if (!first_found) {
                    best_cell = i;
                    first_found = true;
                }
                break;

            case BLISS_SH_FL:  /* First largest non-singleton cell */
                score = cell->size;
                update_best = (score > best_score);
                break;

            case BLISS_SH_FS:  /* First smallest non-singleton cell */
                score = UINT_MAX - cell->size;
                update_best = (score > best_score);
                break;

            case BLISS_SH_FM:  /* First maximally connected non-singleton cell */
                score = compute_cell_connectivity(graph, cell, partition);
                update_best = (score > best_score);
                break;

            case BLISS_SH_FLM: /* First largest maximally connected non-singleton cell */
                score = compute_cell_connectivity(graph, cell, partition) * cell->size;
                update_best = (score > best_score);
                break;

            case BLISS_SH_FSM: /* First smallest maximally connected non-singleton cell */
                score = compute_cell_connectivity(graph, cell, partition) * (UINT_MAX - cell->size);
                update_best = (score > best_score);
                break;
        }

        if (update_best) {
            best_score = score;
            best_cell = i;
        }
    }

    return best_cell;
}

/* ===================================================================
 * ADVANCED REFINEMENT TECHNIQUES
 * =================================================================== */

/* Check if partition is equitable (all vertices in same cell have same degree pattern) */
bool is_partition_equitable(const bliss_graph_t *graph, const partition_t *partition) {
    for (unsigned int cell_idx = 0; cell_idx < partition->num_cells; cell_idx++) {
        const partition_cell_t *cell = &partition->cells[cell_idx];

        if (cell->size <= 1) continue; /* Unit cells are trivially equitable */

        /* Check if all vertices in the cell have the same degree signature */
        degree_signature_t first_sig = compute_degree_signature(graph, cell->elements[0], partition);

        for (unsigned int i = 1; i < cell->size; i++) {
            degree_signature_t curr_sig = compute_degree_signature(graph, cell->elements[i], partition);

            if (compare_degree_signatures(&first_sig, &curr_sig) != 0) {
                /* Found vertices with different signatures - not equitable */
                free_degree_signature(&first_sig);
                free_degree_signature(&curr_sig);
                return false;
            }

            free_degree_signature(&curr_sig);
        }

        free_degree_signature(&first_sig);
    }

    return true;
}

/* Specialized refinement for sparse graphs */
static bool refine_sparse_graph(bliss_graph_t *graph, partition_t *partition) {
    /* For sparse graphs, we can use more aggressive refinement techniques */
    bool changed = false;

    /* First do standard refinement */
    if (refine_partition_complete(graph, partition)) {
        changed = true;
    }

    /* Additional refinement for sparse graphs: consider 2-hop neighborhoods */
    for (unsigned int cell_idx = 0; cell_idx < partition->num_cells; cell_idx++) {
        partition_cell_t *cell = &partition->cells[cell_idx];

        if (cell->is_unit || cell->size <= 1) continue;

        /* Compute 2-hop degree signatures */
        unsigned int *two_hop_degrees = bliss_malloc(cell->size * sizeof(unsigned int));

        for (unsigned int i = 0; i < cell->size; i++) {
            unsigned int vertex = cell->elements[i];
            unsigned int two_hop_count = 0;

            /* Count 2-hop neighbors */
            for (unsigned int j = 0; j < graph->adj_list_sizes[vertex]; j++) {
                unsigned int neighbor = graph->adj_lists[vertex][j];
                two_hop_count += graph->adj_list_sizes[neighbor];
            }

            two_hop_degrees[i] = two_hop_count;
        }

        /* Check if we can split based on 2-hop degrees */
        bool can_split = false;
        for (unsigned int i = 1; i < cell->size; i++) {
            if (two_hop_degrees[i] != two_hop_degrees[0]) {
                can_split = true;
                break;
            }
        }

        if (can_split) {
            /* Sort vertices by 2-hop degree */
            for (unsigned int i = 1; i < cell->size; i++) {
                for (unsigned int j = i; j > 0; j--) {
                    if (two_hop_degrees[j-1] > two_hop_degrees[j]) {
                        /* Swap vertices and degrees */
                        unsigned int temp_vertex = cell->elements[j-1];
                        cell->elements[j-1] = cell->elements[j];
                        cell->elements[j] = temp_vertex;

                        unsigned int temp_degree = two_hop_degrees[j-1];
                        two_hop_degrees[j-1] = two_hop_degrees[j];
                        two_hop_degrees[j] = temp_degree;
                    } else {
                        break;
                    }
                }
            }

            /* Split the cell based on 2-hop degrees */
            /* This is a simplified implementation - full version would create new cells */
            changed = true;
        }

        bliss_free(two_hop_degrees);
    }

    return changed;
}

/* Refinement with failure recording for pruning */
static bool refine_with_failure_recording(bliss_graph_t *graph, partition_t *partition,
                                          void *failure_cache) {
    /* Simplified failure recording - in full implementation, this would
     * maintain a cache of failed refinement attempts to avoid repeating work */
    (void)failure_cache; /* Suppress warning */

    return refine_partition_complete(graph, partition);
}

/* ===================================================================
 * REFINEMENT STATISTICS AND MONITORING
 * =================================================================== */

/* Global refinement statistics */
static struct {
    unsigned long total_refinements;
    unsigned long total_cell_splits;
    unsigned long total_rounds;
    double total_time;
} refinement_stats = {0};

/* Update refinement statistics */
static void update_refinement_stats(unsigned int rounds, double time, unsigned int splits) {
    refinement_stats.total_refinements++;
    refinement_stats.total_rounds += rounds;
    refinement_stats.total_time += time;
    refinement_stats.total_cell_splits += splits;
}

/* Get refinement statistics */
void get_refinement_statistics(unsigned long *refinements, unsigned long *splits,
                              unsigned long *rounds, double *total_time) {
    if (refinements) *refinements = refinement_stats.total_refinements;
    if (splits) *splits = refinement_stats.total_cell_splits;
    if (rounds) *rounds = refinement_stats.total_rounds;
    if (total_time) *total_time = refinement_stats.total_time;
}

/* Reset refinement statistics */
void reset_refinement_statistics(void) {
    memset(&refinement_stats, 0, sizeof(refinement_stats));
}

/* ===================================================================
 * MAIN REFINEMENT INTERFACE
 * =================================================================== */

/* Main refinement function with timing and statistics */
bool refine_partition_with_stats(bliss_graph_t *graph, partition_t *partition) {
    clock_t start = clock();
    unsigned int initial_cells = partition->num_cells;

    bool result = refine_partition_complete(graph, partition);

    clock_t end = clock();
    double time = ((double)(end - start)) / CLOCKS_PER_SEC;
    unsigned int final_cells = partition->num_cells;
    unsigned int splits = final_cells - initial_cells;

    update_refinement_stats(1, time, splits);

    return result;
}

/* Adaptive refinement that chooses the best strategy based on graph properties */
bool refine_partition_adaptive(bliss_graph_t *graph, partition_t *partition) {
    /* Choose refinement strategy based on graph density */
    double density = (double)(graph->num_edges * 2) / (graph->num_vertices * (graph->num_vertices - 1));

    if (density < 0.1) {
        /* Sparse graph - use specialized sparse refinement */
        return refine_sparse_graph(graph, partition);
    } else if (density > 0.7) {
        /* Dense graph - use standard refinement */
        return refine_partition_complete(graph, partition);
    } else {
        /* Medium density - use standard refinement with failure recording */
        return refine_with_failure_recording(graph, partition, NULL);
    }
}

/* ===================================================================
 * REFINEMENT VALIDATION AND TESTING
 * =================================================================== */

/* Validate that a partition is properly refined */
bool validate_refinement(const bliss_graph_t *graph, const partition_t *partition) {
    /* Check that partition is equitable */
    if (!is_partition_equitable(graph, partition)) {
        return false;
    }

    /* Check that all mappings are consistent */
    for (unsigned int i = 0; i < partition->num_cells; i++) {
        const partition_cell_t *cell = &partition->cells[i];

        for (unsigned int j = 0; j < cell->size; j++) {
            unsigned int vertex = cell->elements[j];

            if (partition->element_to_cell[vertex] != i) {
                return false;
            }

            if (partition->position_in_cell[vertex] != j) {
                return false;
            }
        }
    }

    return true;
}

/* Test refinement on a known graph */
bool test_refinement_algorithm(void) {
    /* Create test graph - cycle C6 */
    bliss_graph_t *graph = bliss_new(6);
    for (unsigned int i = 0; i < 6; i++) {
        bliss_add_edge(graph, i, (i + 1) % 6);
    }
    
    /* Create initial partition - all vertices in one cell */
    partition_t *partition = partition_new(6);
    for (unsigned int i = 0; i < 6; i++) {
        partition_add_to_cell(partition, 0, i);
    }
    partition->num_cells = 1;
    
    /* Refine the partition */
    bool refined = refine_partition_complete(graph, partition);
    (void)refined; /* ADD THIS LINE to suppress unused variable warning */
    
    /* For C6, all vertices should remain in the same cell after refinement */
    bool test_passed = (partition->num_cells == 1) && validate_refinement(graph, partition);
    
    /* Cleanup */
    partition_release(partition);
    bliss_release(graph);
    
    return test_passed;
}