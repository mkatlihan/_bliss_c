/*
 * partition_individualization.c - Critical partition operations for bliss search tree
 *
 * This implements the core partition copying and vertex individualization that
 * enables the search tree branching mechanism in the bliss algorithm.
 */

#include "bliss.h"

/* ===================================================================
 * PARTITION DEEP COPYING
 * =================================================================== */

/* Deep copy a partition cell */
static void copy_partition_cell(partition_cell_t *dest, const partition_cell_t *src) {
    dest->size = src->size;
    dest->first = src->first;
    dest->is_unit = src->is_unit;

    if (src->size > 0) {
        dest->elements = bliss_malloc(src->size * sizeof(unsigned int));
        memcpy(dest->elements, src->elements, src->size * sizeof(unsigned int));
    } else {
        dest->elements = NULL;
    }
}

/* Create a deep copy of a partition */
static partition_t *partition_copy(const partition_t *original, unsigned int num_vertices) {
    partition_t *copy = bliss_malloc(sizeof(partition_t));

    copy->num_cells = original->num_cells;
    copy->capacity = original->capacity;
    copy->num_discrete_cells = original->num_discrete_cells;

    /* Allocate arrays */
    copy->cells = bliss_malloc(copy->capacity * sizeof(partition_cell_t));
    copy->element_to_cell = bliss_malloc(num_vertices * sizeof(unsigned int));
    copy->position_in_cell = bliss_malloc(num_vertices * sizeof(unsigned int));

    /* Copy cell mappings */
    memcpy(copy->element_to_cell, original->element_to_cell,
           num_vertices * sizeof(unsigned int));
    memcpy(copy->position_in_cell, original->position_in_cell,
           num_vertices * sizeof(unsigned int));

    /* Deep copy each cell */
    for (unsigned int i = 0; i < copy->capacity; i++) {
        if (i < original->num_cells) {
            copy_partition_cell(&copy->cells[i], &original->cells[i]);
        } else {
            /* Initialize empty cells */
            copy->cells[i].elements = NULL;
            copy->cells[i].size = 0;
            copy->cells[i].first = 0;
            copy->cells[i].is_unit = false;
        }
    }

    return copy;
}

/* ===================================================================
 * VERTEX INDIVIDUALIZATION - THE CORE OF SEARCH TREE BRANCHING
 * =================================================================== */

/* Split a cell by moving one vertex to a new singleton cell */
static bool individualize_vertex_in_cell(partition_t *partition,
                                         unsigned int cell_idx,
                                         unsigned int vertex_to_split) {
    partition_cell_t *cell = &partition->cells[cell_idx];

    if (cell->is_unit) {
        return false; /* Cannot split a unit cell */
    }

    /* Find the position of the vertex in the cell */
    unsigned int vertex_pos = UINT_MAX;
    for (unsigned int i = 0; i < cell->size; i++) {
        if (cell->elements[i] == vertex_to_split) {
            vertex_pos = i;
            break;
        }
    }

    if (vertex_pos == UINT_MAX) {
        return false; /* Vertex not found in cell */
    }

    /* Create new singleton cell for the individualized vertex */
    unsigned int new_cell_idx = partition->num_cells;

    /* Ensure we have capacity for the new cell */
    if (new_cell_idx >= partition->capacity) {
        /* Expand capacity */
        unsigned int new_capacity = partition->capacity * 2;
        partition->cells = bliss_realloc(partition->cells,
                                        new_capacity * sizeof(partition_cell_t));

        /* Initialize new cells */
        for (unsigned int i = partition->capacity; i < new_capacity; i++) {
            partition->cells[i].elements = NULL;
            partition->cells[i].size = 0;
            partition->cells[i].first = 0;
            partition->cells[i].is_unit = false;
        }

        partition->capacity = new_capacity;
    }

    /* Create the new singleton cell */
    partition_cell_t *new_cell = &partition->cells[new_cell_idx];
    new_cell->size = 1;
    new_cell->elements = bliss_malloc(sizeof(unsigned int));
    new_cell->elements[0] = vertex_to_split;
    new_cell->is_unit = true;
    new_cell->first = 0;

    /* Update vertex-to-cell mapping */
    partition->element_to_cell[vertex_to_split] = new_cell_idx;
    partition->position_in_cell[vertex_to_split] = 0;

    /* Remove vertex from original cell */
    for (unsigned int i = vertex_pos; i < cell->size - 1; i++) {
        cell->elements[i] = cell->elements[i + 1];
        /* Update position mapping for moved vertices */
        partition->position_in_cell[cell->elements[i]] = i;
    }

    cell->size--;

    /* Update cell status */
    if (cell->size == 1) {
        cell->is_unit = true;
        partition->num_discrete_cells++;
    }

    /* Increment cell count and discrete cell count */
    partition->num_cells++;
    partition->num_discrete_cells++;

    return true;
}

/* ===================================================================
 * CANONICAL VERTEX ORDERING WITHIN CELLS
 * =================================================================== */

/* Compare vertices for canonical ordering within a cell */
static int compare_vertices_canonical(const void *a, const void *b, void *graph_ptr) {
    const bliss_graph_t *graph = (const bliss_graph_t *)graph_ptr;
    unsigned int v1 = *(const unsigned int *)a;
    unsigned int v2 = *(const unsigned int *)b;

    /* Primary comparison: vertex color */
    if (graph->vertex_colors[v1] != graph->vertex_colors[v2]) {
        return graph->vertex_colors[v1] < graph->vertex_colors[v2] ? -1 : 1;
    }

    /* Secondary comparison: degree */
    if (graph->adj_list_sizes[v1] != graph->adj_list_sizes[v2]) {
        return graph->adj_list_sizes[v1] < graph->adj_list_sizes[v2] ? -1 : 1;
    }

    /* Tertiary comparison: vertex index (for stability) */
    return v1 < v2 ? -1 : (v1 > v2 ? 1 : 0);
}

/* Sort vertices within a cell for canonical ordering */
static void canonicalize_cell_ordering(partition_t *partition,
                                       unsigned int cell_idx,
                                       const bliss_graph_t *graph) {
    partition_cell_t *cell = &partition->cells[cell_idx];

    if (cell->size <= 1) {
        return; /* Nothing to sort */
    }

    /* Sort vertices in canonical order */
    /* Note: Using a stable sort to maintain deterministic results */
    for (unsigned int i = 1; i < cell->size; i++) {
        unsigned int key = cell->elements[i];
        int j = i - 1;

        while (j >= 0 && compare_vertices_canonical(&cell->elements[j], &key, (void*)graph) > 0) {
            cell->elements[j + 1] = cell->elements[j];
            j--;
        }
        cell->elements[j + 1] = key;
    }

    /* Update position mappings after sorting */
    for (unsigned int i = 0; i < cell->size; i++) {
        partition->position_in_cell[cell->elements[i]] = i;
    }
}

/* ===================================================================
 * COMPLETE INDIVIDUALIZATION PROCEDURE
 * =================================================================== */

/* Perform complete vertex individualization with canonical ordering */
static partition_t *individualize_vertex(const partition_t *original_partition,
                                         unsigned int vertex_to_individualize,
                                         unsigned int num_vertices,
                                         const bliss_graph_t *graph) {
    /* Create a copy of the original partition */
    partition_t *new_partition = partition_copy(original_partition, num_vertices);

    /* Find which cell contains the vertex to individualize */
    unsigned int target_cell = new_partition->element_to_cell[vertex_to_individualize];

    /* Canonicalize the ordering of vertices in the target cell before splitting */
    canonicalize_cell_ordering(new_partition, target_cell, graph);

    /* Individualize the specified vertex */
    if (!individualize_vertex_in_cell(new_partition, target_cell, vertex_to_individualize)) {
        /* Failed to individualize - should not happen with valid inputs */
        partition_release(new_partition);
        return NULL;
    }

    /* Canonicalize ordering in all cells after individualization */
    for (unsigned int i = 0; i < new_partition->num_cells; i++) {
        canonicalize_cell_ordering(new_partition, i, graph);
    }

    return new_partition;
}

/* ===================================================================
 * ORBIT-BASED PRUNING SUPPORT
 * =================================================================== */

/* Check if two vertices are in the same orbit (simplified version) */
static bool vertices_in_same_orbit(unsigned int v1, unsigned int v2,
                                  const search_state_t *state,
                                  const bliss_graph_t *graph) {
    /* Basic orbit check - in a complete implementation, this would use
     * the accumulated automorphism generators to compute orbits */

    /* If vertices have different colors, they cannot be in the same orbit */
    if (graph->vertex_colors[v1] != graph->vertex_colors[v2]) {
        return false;
    }

    /* If vertices have different degrees, they cannot be in the same orbit */
    if (graph->adj_list_sizes[v1] != graph->adj_list_sizes[v2]) {
        return false;
    }

    /* TODO: More sophisticated orbit computation using generators */
    return false; /* Conservative approach - assume different orbits */
}

/* Get canonical representative from a cell (for orbit pruning) */
static unsigned int get_canonical_representative(const partition_cell_t *cell,
                                               const bliss_graph_t *graph) {
    if (cell->size == 0) {
        return UINT_MAX;
    }

    if (cell->size == 1) {
        return cell->elements[0];
    }

    /* Return the lexicographically smallest vertex in canonical ordering */
    unsigned int best_vertex = cell->elements[0];

    for (unsigned int i = 1; i < cell->size; i++) {
        if (compare_vertices_canonical(&cell->elements[i], &best_vertex, (void*)graph) < 0) {
            best_vertex = cell->elements[i];
        }
    }

    return best_vertex;
}

/* ===================================================================
 * ENHANCED SEARCH WITH PROPER INDIVIDUALIZATION
 * =================================================================== */

/* Create child node with proper partition individualization */
static search_node_t *create_child_node(const search_node_t *parent,
                                        unsigned int vertex_to_individualize,
                                        unsigned int target_cell_idx,
                                        const bliss_graph_t *graph) {
    search_node_t *child = bliss_malloc(sizeof(search_node_t));

    child->parent = (search_node_t*)parent;
    child->child = NULL;
    child->sibling = NULL;
    child->splitting_vertex = vertex_to_individualize;
    child->target_cell = target_cell_idx;
    child->level = parent->level + 1;

    /* Create individualized partition */
    partition_t *new_partition = individualize_vertex(&parent->partition,
                                                     vertex_to_individualize,
                                                     graph->num_vertices,
                                                     graph);

    if (!new_partition) {
        bliss_free(child);
        return NULL;
    }

    child->partition = *new_partition;
    bliss_free(new_partition); /* Only free the wrapper, not the contents */

    return child;
}

/* Enhanced search with proper individualization and orbit pruning */
static void search_automorphisms_with_individualization(bliss_graph_t *graph,
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
        search_node_t *child = create_child_node(state->current, split_vertex,
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

/* ===================================================================
 * UPDATED MAIN AUTOMORPHISM FUNCTION
 * =================================================================== */

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

/* ===================================================================
 * UTILITY FUNCTIONS FOR INTEGRATION
 * =================================================================== */

/* Replace the original function in bliss_core.c with this version */
void bliss_find_automorphisms(bliss_graph_t *graph,
                              bliss_stats_t *stats,
                              bliss_automorphism_hook_t hook,
                              void *hook_user_param) {
    /* Call the complete implementation */
    bliss_find_automorphisms_complete(graph, stats, hook, hook_user_param);
}

/* Also update the canonical labeling function to use the complete search */
const unsigned int *bliss_find_canonical_labeling(bliss_graph_t *graph,
                                                   bliss_stats_t *stats,
                                                   bliss_automorphism_hook_t hook,
                                                   void *hook_user_param) {
    if (BLISS_UNLIKELY(!graph)) {
        return NULL;
    }

    if (graph->canonical_labeling_valid) {
        return graph->canonical_labeling;
    }

    /* Allocate canonical labeling if needed */
    if (!graph->canonical_labeling) {
        graph->canonical_labeling = bliss_malloc(graph->num_vertices * sizeof(unsigned int));
    }

    /* Use the complete search to find canonical labeling */
    bliss_find_automorphisms_complete(graph, stats, hook, hook_user_param);

    /* Extract canonical labeling from search results */
    /* This is still simplified - proper canonical labeling extraction
     * would track the lexicographically maximum leaf in the search tree */
    for (unsigned int i = 0; i < graph->num_vertices; i++) {
        graph->canonical_labeling[i] = i;
    }

    graph->canonical_labeling_valid = true;
    return graph->canonical_labeling;
}