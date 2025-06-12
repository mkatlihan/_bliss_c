/*
 * critical_fixes.c - Specific fixes for the identified crash-causing bugs
 */

#include "bliss.h"

/* ===================================================================
 * FIX 1: Partition Memory Management (Most Critical)
 * =================================================================== */

// REPLACE the create_child_node_unified function in bliss_unified.c with this:
static search_node_t *create_child_node_unified_fixed(const search_node_t *parent,
                                                      unsigned int vertex_to_individualize,
                                                      unsigned int target_cell_idx,
                                                      const bliss_graph_t *graph) {
    search_node_t *child = bliss_malloc(sizeof(search_node_t));
    if (!child) {
        return NULL;
    }

    child->parent = (search_node_t*)parent;
    child->child = NULL;
    child->sibling = NULL;
    child->splitting_vertex = vertex_to_individualize;
    child->target_cell = target_cell_idx;
    child->level = parent->level + 1;

    /* CRITICAL FIX: Proper partition individualization with error checking */
    partition_t *new_partition = individualize_vertex(&parent->partition,
                                                     vertex_to_individualize,
                                                     graph->num_vertices,
                                                     graph);

    if (!new_partition) {
        bliss_free(child);
        return NULL;
    }

    /* CORRECT: Transfer ownership of arrays, not shallow copy */
    child->partition.capacity = new_partition->capacity;
    child->partition.num_cells = new_partition->num_cells;
    child->partition.num_discrete_cells = new_partition->num_discrete_cells;

    /* Transfer ownership of allocated arrays */
    child->partition.cells = new_partition->cells;
    child->partition.element_to_cell = new_partition->element_to_cell;
    child->partition.position_in_cell = new_partition->position_in_cell;

    /* Clear the source pointers to prevent double-free */
    new_partition->cells = NULL;
    new_partition->element_to_cell = NULL;
    new_partition->position_in_cell = NULL;

    /* Now safe to free the wrapper */
    bliss_free(new_partition);

    return child;
}

/* ===================================================================
 * FIX 2: Safe Partition Copying with Error Handling
 * =================================================================== */

// REPLACE the copy_partition_cell function in bliss_partition.c:
static bool copy_partition_cell_safe(partition_cell_t *dest, const partition_cell_t *src) {
    dest->size = src->size;
    dest->first = src->first;
    dest->is_unit = src->is_unit;

    if (src->size > 0) {
        dest->elements = bliss_malloc(src->size * sizeof(unsigned int));
        if (!dest->elements) {
            return false; // PROPAGATE FAILURE
        }
        memcpy(dest->elements, src->elements, src->size * sizeof(unsigned int));
    } else {
        dest->elements = NULL;
    }

    return true; // SUCCESS
}

// REPLACE the partition_copy function in bliss_partition.c:
partition_t *partition_copy_safe(const partition_t *original, unsigned int num_vertices) {
    if (!original) {
        return NULL;
    }

    partition_t *copy = bliss_malloc(sizeof(partition_t));
    if (!copy) {
        return NULL;
    }

    copy->num_cells = original->num_cells;
    copy->capacity = original->capacity;
    copy->num_discrete_cells = original->num_discrete_cells;

    /* Allocate arrays */
    copy->cells = bliss_malloc(copy->capacity * sizeof(partition_cell_t));
    copy->element_to_cell = bliss_malloc(num_vertices * sizeof(unsigned int));
    copy->position_in_cell = bliss_malloc(num_vertices * sizeof(unsigned int));

    if (!copy->cells || !copy->element_to_cell || !copy->position_in_cell) {
        bliss_free(copy->cells);
        bliss_free(copy->element_to_cell);
        bliss_free(copy->position_in_cell);
        bliss_free(copy);
        return NULL;
    }

    /* Copy cell mappings */
    memcpy(copy->element_to_cell, original->element_to_cell,
           num_vertices * sizeof(unsigned int));
    memcpy(copy->position_in_cell, original->position_in_cell,
           num_vertices * sizeof(unsigned int));

    /* Deep copy each cell with error checking */
    for (unsigned int i = 0; i < copy->capacity; i++) {
        if (i < original->num_cells && original->cells[i].size > 0) {
            if (!copy_partition_cell_safe(&copy->cells[i], &original->cells[i])) {
                /* CLEANUP ON FAILURE */
                for (unsigned int j = 0; j < i; j++) {
                    bliss_free(copy->cells[j].elements);
                }
                bliss_free(copy->cells);
                bliss_free(copy->element_to_cell);
                bliss_free(copy->position_in_cell);
                bliss_free(copy);
                return NULL;
            }
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
 * FIX 3: Safe Partition Operations with Error Propagation
 * =================================================================== */

// REPLACE partition_add_to_cell in bliss_partition.c:
bool partition_add_to_cell_safe(partition_t *partition, unsigned int cell_idx, unsigned int vertex) {
    if (!partition || cell_idx >= partition->capacity) {
        return false; // PROPAGATE ERROR
    }

    partition_cell_t *cell = &partition->cells[cell_idx];

    if (cell->size == 0) {
        cell->elements = bliss_malloc(8 * sizeof(unsigned int));
        if (!cell->elements) {
            return false; // PROPAGATE ERROR
        }
    } else {
        unsigned int capacity = 8;
        while (capacity <= cell->size) capacity *= 2;

        if (cell->size > 0 && (cell->size & (cell->size - 1)) == 0 && cell->size >= 8) {
            unsigned int *new_elements = bliss_realloc(cell->elements,
                                                      capacity * sizeof(unsigned int));
            if (!new_elements) {
                return false; // PROPAGATE ERROR
            }
            cell->elements = new_elements;
        }
    }

    cell->elements[cell->size] = vertex;
    partition->element_to_cell[vertex] = cell_idx;
    partition->position_in_cell[vertex] = cell->size;
    cell->size++;

    cell->is_unit = (cell->size == 1);

    if (cell->is_unit && cell->size == 1) {
        partition->num_discrete_cells++;
    } else if (cell->size == 2) {
        partition->num_discrete_cells--;
    }

    return true; // SUCCESS
}

/* ===================================================================
 * FIX 4: Safe Refinement Loop
 * =================================================================== */

// REPLACE refine_partition_complete in bliss_refinement.c:
bool refine_partition_complete_safe(bliss_graph_t *graph, partition_t *partition) {
    if (!graph || !partition) {
        return false;
    }

    bool global_changed = false;
    bool local_changed = true;
    unsigned int round = 0;

    while (local_changed && round < MAX_REFINEMENT_ROUNDS) {
        local_changed = refine_partition_one_round(graph, partition);

        if (local_changed) {
            global_changed = true;

            /* FIXED: Snapshot num_cells before canonicalization */
            unsigned int cells_to_canonicalize = partition->num_cells;
            for (unsigned int i = 0; i < cells_to_canonicalize; i++) {
                canonicalize_cell_ordering(partition, i, graph);
            }
        }

        round++;
    }

    /* Update discrete cell count safely */
    partition->num_discrete_cells = 0;
    for (unsigned int i = 0; i < partition->num_cells; i++) {
        if (partition->cells[i].is_unit) {
            partition->num_discrete_cells++;
        }
    }

    return global_changed;
}

/* ===================================================================
 * FIX 5: Safe Generator Storage
 * =================================================================== */

// ADD this function to bliss_unified.c:
static bool expand_generator_storage(search_state_t *state) {
    if (!state) {
        return false;
    }

    unsigned int new_capacity = state->generator_capacity * 2;
    unsigned int **new_generators = bliss_realloc(state->generators,
                                                  new_capacity * sizeof(unsigned int*));

    if (!new_generators) {
        return false; // OUT OF MEMORY
    }

    state->generators = new_generators;
    state->generator_capacity = new_capacity;
    return true;
}

// REPLACE generator storage code in search_automorphisms_unified:
if (!is_identity) {
    /* Add to generator set with safe expansion */
    if (state->num_generators >= state->generator_capacity) {
        if (!expand_generator_storage(state)) {
            bliss_free(automorphism); // Clean up on failure
            return; // Can't store more generators
        }
    }

    state->generators[state->num_generators] = automorphism;
    state->num_generators++;
    state->stats->nof_generators++;
} else {
    bliss_free(automorphism);
}

/* ===================================================================
 * FIX 6: Partition Validation for Search Loop
 * =================================================================== */

// ADD this validation function:
static bool validate_partition_for_search(const partition_t *partition, unsigned int num_vertices) {
    if (!partition || !partition->cells || !partition->element_to_cell || !partition->position_in_cell) {
        return false;
    }

    // Quick sanity checks
    if (partition->num_cells > partition->capacity) {
        return false;
    }

    if (partition->num_discrete_cells > num_vertices) {
        return false;
    }

    // Check that used cells have valid elements arrays
    for (unsigned int i = 0; i < partition->num_cells; i++) {
        const partition_cell_t *cell = &partition->cells[i];
        if (cell->size > 0 && !cell->elements) {
            return false;
        }
    }

    return true;
}

// MODIFY the search loop in search_automorphisms_unified:
/* STEP 6: Create child node with FIXED individualization */
search_node_t *child = create_child_node_unified_fixed(state->current, split_vertex,
                                                      target_cell, graph);

if (!child) {
    continue; /* Skip if individualization failed */
}

/* CRITICAL: Validate child partition before proceeding */
if (!validate_partition_for_search(&child->partition, graph->num_vertices)) {
    printf("ERROR: Invalid child partition for vertex %u\n", split_vertex);
    partition_release(&child->partition);
    bliss_free(child);
    continue;
}

/* Update maximum level reached */
if (child->level > state->stats->max_level) {
    state->stats->max_level = child->level;
}

/* STEP 7: Recursively search the child */
search_node_t *prev_current = state->current;
state->current = child;
search_automorphisms_unified(graph, state, hook, hook_param);
state->current = prev_current;

/* STEP 8: Clean up child node safely */
partition_release(&child->partition);
bliss_free(child);

/* ===================================================================
 * INTEGRATION INSTRUCTIONS
 * =================================================================== */

/*
 * To apply these fixes:
 *
 * 1. Replace create_child_node_unified with create_child_node_unified_fixed
 * 2. Replace partition_copy with partition_copy_safe
 * 3. Replace partition_add_to_cell with partition_add_to_cell_safe
 * 4. Replace refine_partition_complete with refine_partition_complete_safe
 * 5. Add generator storage expansion and validation functions
 * 6. Modify search loop to use validation
 *
 * These fixes address the root causes of crashes in advanced automorphism tests.
 */