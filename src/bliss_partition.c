/*
 * bliss_partition.c - Partition copying and individualization operations
 *
 * This file contains the core partition manipulation functions that enable
 * the search tree branching mechanism in the bliss algorithm.
 */

#include "bliss.h"
#include <time.h>


/* ===================================================================
 * IMPLEMENTATION OF BASIC PARTITION FUNCTIONS
 * =================================================================== */

/* Remove the static inline versions and use these instead: */

partition_t *partition_new(unsigned int num_vertices) {
    partition_t *p = bliss_malloc(sizeof(partition_t));
    if (!p) return NULL;
    
    p->capacity = num_vertices + 1;
    p->cells = bliss_malloc(p->capacity * sizeof(partition_cell_t));
    p->element_to_cell = bliss_malloc(num_vertices * sizeof(unsigned int));
    p->position_in_cell = bliss_malloc(num_vertices * sizeof(unsigned int));
    p->num_cells = 0;
    p->num_discrete_cells = 0;
    
    if (!p->cells || !p->element_to_cell || !p->position_in_cell) {
        bliss_free(p->cells);
        bliss_free(p->element_to_cell);
        bliss_free(p->position_in_cell);
        bliss_free(p);
        return NULL;
    }
    
    /* Initialize arrays */
    for (unsigned int i = 0; i < p->capacity; i++) {
        p->cells[i].elements = NULL;
        p->cells[i].size = 0;
        p->cells[i].first = 0;
        p->cells[i].is_unit = false;
    }
    
    for (unsigned int i = 0; i < num_vertices; i++) {
        p->element_to_cell[i] = UINT_MAX;
        p->position_in_cell[i] = UINT_MAX;
    }
    
    return p;
}

void partition_release(partition_t *partition) {
    if (!partition) return;
    
    for (unsigned int i = 0; i < partition->capacity; i++) {
        bliss_free(partition->cells[i].elements);
    }
    
    bliss_free(partition->cells);
    bliss_free(partition->element_to_cell);
    bliss_free(partition->position_in_cell);
    bliss_free(partition);
}

void partition_add_to_cell(partition_t *partition, unsigned int cell_idx, unsigned int vertex) {
    if (cell_idx >= partition->capacity) {
        return;
    }
    
    partition_cell_t *cell = &partition->cells[cell_idx];
    
    if (cell->size == 0) {
        cell->elements = bliss_malloc(8 * sizeof(unsigned int));
        if (!cell->elements) return;
    } else {
        unsigned int capacity = 8;
        while (capacity <= cell->size) capacity *= 2;
        
        if (cell->size > 0 && (cell->size & (cell->size - 1)) == 0 && cell->size >= 8) {
            cell->elements = bliss_realloc(cell->elements, 
                                          capacity * sizeof(unsigned int));
            if (!cell->elements) return;
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
}

bool partition_is_discrete(const partition_t *partition, unsigned int num_vertices) {
    return partition && (partition->num_discrete_cells == num_vertices);
}

/* Get the cell index of a vertex */
unsigned int partition_get_cell_of_vertex(const partition_t *partition, unsigned int vertex) {
    if (!partition || vertex >= partition->capacity) {
        return UINT_MAX;
    }
    return partition->element_to_cell[vertex];
}

/* Get the position of a vertex within its cell */
unsigned int partition_get_position_in_cell(const partition_t *partition, unsigned int vertex) {
    if (!partition || vertex >= partition->capacity) {
        return UINT_MAX;
    }
    return partition->position_in_cell[vertex];
}

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
partition_t *partition_copy(const partition_t *original, unsigned int num_vertices) {
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
 * VERTEX INDIVIDUALIZATION
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
 * CANONICAL VERTEX ORDERING
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
void canonicalize_cell_ordering(partition_t *partition,
                                unsigned int cell_idx,
                                const bliss_graph_t *graph) {
    partition_cell_t *cell = &partition->cells[cell_idx];

    if (cell->size <= 1) {
        return; /* Nothing to sort */
    }

    /* Sort vertices in canonical order using insertion sort for stability */
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
partition_t *individualize_vertex(const partition_t *original_partition,
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
bool vertices_in_same_orbit(unsigned int v1, unsigned int v2,
                           const search_state_t *state,
                           const bliss_graph_t *graph) {
    /* Suppress unused parameter warning */
    (void)state; /* ADD THIS LINE */
    
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
    /* For now, conservatively assume different orbits */
    return false;
}


/* Get canonical representative from a cell (for orbit pruning) */
unsigned int get_canonical_representative(const partition_cell_t *cell,
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
 * PARTITION VALIDATION AND DEBUGGING
 * =================================================================== */

/* Validate partition consistency (useful for debugging) */
bool validate_partition(const partition_t *partition, unsigned int num_vertices) {
    if (!partition) return false;

    /* Check that all vertices are accounted for exactly once */
    bool *vertex_seen = bliss_malloc(num_vertices * sizeof(bool));
    memset(vertex_seen, 0, num_vertices * sizeof(bool));

    unsigned int total_vertices = 0;
    unsigned int discrete_count = 0;

    for (unsigned int cell_idx = 0; cell_idx < partition->num_cells; cell_idx++) {
        const partition_cell_t *cell = &partition->cells[cell_idx];

        total_vertices += cell->size;
        if (cell->is_unit) discrete_count++;

        /* Check each vertex in the cell */
        for (unsigned int i = 0; i < cell->size; i++) {
            unsigned int vertex = cell->elements[i];

            /* Check vertex is valid */
            if (vertex >= num_vertices) {
                bliss_free(vertex_seen);
                return false;
            }

            /* Check vertex hasn't been seen before */
            if (vertex_seen[vertex]) {
                bliss_free(vertex_seen);
                return false;
            }
            vertex_seen[vertex] = true;

            /* Check mappings are consistent */
            if (partition->element_to_cell[vertex] != cell_idx) {
                bliss_free(vertex_seen);
                return false;
            }

            if (partition->position_in_cell[vertex] != i) {
                bliss_free(vertex_seen);
                return false;
            }
        }
    }

    /* Check all vertices are accounted for */
    bool all_vertices_present = (total_vertices == num_vertices);
    for (unsigned int i = 0; i < num_vertices; i++) {
        if (!vertex_seen[i]) {
            all_vertices_present = false;
            break;
        }
    }

    /* Check discrete cell count */
    bool discrete_count_correct = (discrete_count == partition->num_discrete_cells);

    bliss_free(vertex_seen);

    return all_vertices_present && discrete_count_correct;
}

/* Print partition for debugging */
void print_partition(const partition_t *partition, FILE *fp) {
    if (!partition || !fp) return;

    fprintf(fp, "Partition with %u cells (%u discrete):\n",
            partition->num_cells, partition->num_discrete_cells);

    for (unsigned int i = 0; i < partition->num_cells; i++) {
        const partition_cell_t *cell = &partition->cells[i];
        fprintf(fp, "  Cell %u (%s): {", i, cell->is_unit ? "unit" : "non-unit");

        for (unsigned int j = 0; j < cell->size; j++) {
            fprintf(fp, "%u", cell->elements[j]);
            if (j < cell->size - 1) fprintf(fp, ", ");
        }
        fprintf(fp, "}\n");
    }
}

/* Extract labeling from a discrete partition */
void extract_labeling_from_partition(const partition_t *partition,
                                     unsigned int *labeling,
                                     unsigned int n) {
    /* Initialize labeling array */
    for (unsigned int i = 0; i < n; i++) {
        labeling[i] = UINT_MAX;
    }

    /* Extract vertices from unit cells in order */
    unsigned int position = 0;
    for (unsigned int cell_idx = 0; cell_idx < partition->num_cells && position < n; cell_idx++) {
        const partition_cell_t *cell = &partition->cells[cell_idx];
        if (cell->is_unit && cell->size > 0) {
            labeling[position] = cell->elements[0];
            position++;
        }
    }

    /* Handle non-discrete partitions by filling remaining positions */
    if (position < n) {
        for (unsigned int cell_idx = 0; cell_idx < partition->num_cells; cell_idx++) {
            const partition_cell_t *cell = &partition->cells[cell_idx];
            if (!cell->is_unit) {
                for (unsigned int i = 0; i < cell->size && position < n; i++) {
                    labeling[position] = cell->elements[i];
                    position++;
                }
            }
        }
    }
}