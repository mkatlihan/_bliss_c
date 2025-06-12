# Crash Analysis: Advanced Automorphism Tests

After examining your actual code files, here are the specific issues causing crashes in advanced tests:

## Critical Memory Corruption Issues

### 1. Partition Memory Management in `bliss_unified.c` (Lines 160-166)

**The Bug:**
```c
// create_child_node_unified() in bliss_unified.c:160-166
partition_t *new_partition = individualize_vertex(&parent->partition,
                                                 vertex_to_individualize,
                                                 graph->num_vertices,
                                                 graph);
if (!new_partition) {
    bliss_free(child);
    return NULL;
}

child->partition = *new_partition;        // SHALLOW COPY - DANGEROUS!
bliss_free(new_partition);               // FREES WRAPPER, LEAVES DANGLING POINTERS
```

**Why it crashes:** The `individualize_vertex()` returns a partition with allocated arrays (`cells`, `element_to_cell`, `position_in_cell`). The shallow copy `*new_partition` copies the pointer values, then `bliss_free(new_partition)` frees the wrapper struct but leaves the child partition pointing to the same memory. When the parent scope ends or another partition is created, these arrays get corrupted or freed multiple times.

**Fix:** Either deep copy the arrays or transfer ownership properly.

### 2. Partition Copying Bug in `bliss_partition.c` (Lines 124-153)

**The Bug:**
```c
// partition_copy() function doesn't handle failure cases properly
partition_t *copy = bliss_malloc(sizeof(partition_t));
// ... allocations ...
/* Deep copy each cell */
for (unsigned int i = 0; i < copy->capacity; i++) {
    if (i < original->num_cells) {
        copy_partition_cell(&copy->cells[i], &original->cells[i]);  // CAN FAIL
    } else {
        // Initialize empty cells
    }
}
```

**Why it crashes:** The `copy_partition_cell()` can fail (line 104 allocates memory), but there's no error checking. If allocation fails partway through, you have a partially constructed partition with some valid and some invalid cells.

### 3. Cell Element Array Management in `bliss_partition.c` (Lines 65-95)

**The Bug:**
```c
// partition_add_to_cell() has insufficient bounds checking
void partition_add_to_cell(partition_t *partition, unsigned int cell_idx, unsigned int vertex) {
    if (cell_idx >= partition->capacity) {
        return;  // SILENT FAILURE - NO ERROR PROPAGATION
    }
    
    partition_cell_t *cell = &partition->cells[cell_idx];
    
    if (cell->size == 0) {
        cell->elements = bliss_malloc(8 * sizeof(unsigned int));
        if (!cell->elements) return;  // SILENT FAILURE - PARTITION NOW INCONSISTENT
    }
    // ... growth logic that can fail silently
}
```

**Why it crashes:** Silent failures leave the partition in an inconsistent state. Later code assumes elements arrays are valid when they might be NULL.

## Algorithm Logic Errors

### 4. Incomplete Refinement in `bliss_refinement.c` (Lines 80-150)

**The Bug:**
```c
// split_cell_by_signatures() has logical issues
static bool split_cell_by_signatures(partition_t *partition, unsigned int cell_idx,
                                     const bliss_graph_t *graph) {
    // ... sorting code ...
    
    /* Find split points where signatures differ */
    bool split_occurred = false;
    unsigned int group_start = 0;
    
    for (unsigned int i = 1; i <= cell->size; i++) {
        bool should_split = (i == cell->size) || 
                           (compare_degree_signatures(&signatures[group_start], &signatures[i]) != 0);
        
        if (should_split && i > group_start + 1) {
            /* Found a group with multiple vertices - split if not the first group */
            if (group_start > 0) {  // ONLY SPLITS NON-FIRST GROUPS
                // Create new cell...
            }
        }
    }
}
```

**Why it crashes:** The splitting logic only creates new cells for non-first groups. This means the original cell retains mixed vertices, leading to incorrect partition states that confuse the search algorithm.

### 5. Missing Validation in Search Loop (`bliss_unified.c` Lines 200-250)

**The Bug:**
```c
// search_automorphisms_unified() doesn't validate partitions
for (unsigned int i = 0; i < split_cell->size; i++) {
    unsigned int split_vertex = split_cell->elements[i];
    
    // ... orbit pruning ...
    
    search_node_t *child = create_child_node_unified(state->current, split_vertex, 
                                                    target_cell, graph);
    
    if (!child) {
        continue; // Skip if individualization failed
    }
    
    // NO VALIDATION OF child->partition HERE
    
    search_node_t *prev_current = state->current;
    state->current = child;
    search_automorphisms_unified(graph, state, hook, hook_param);  // RECURSIVE CALL WITH BAD PARTITION
    state->current = prev_current;
}
```

**Why it crashes:** The recursive call proceeds even if the child partition is corrupted from the memory management bug above.

## Test-Specific Issues

### 6. Complex Graph Test Failures (`bliss_tests.c` Lines 250-300)

**The Bug:**
Your test for Petersen graph and complex automorphism detection:
```c
// test_complex_automorphisms() 
bliss_graph_t *petersen = bliss_create_petersen_graph();
bliss_find_automorphisms(petersen, stats, count_automorphisms, &counter);
```

**Why it crashes:** The Petersen graph has 120 automorphisms. The recursive search creates deep partition trees, amplifying the memory corruption issues. Each level has the potential partition copying bug, leading to cascade failures.

### 7. Generator Storage Overflow (`bliss_unified.c` Lines 120-140)

**The Bug:**
```c
// In handle automorphism storage
if (state->num_generators >= state->generator_capacity) {
    state->generator_capacity *= 2;
    state->generators = bliss_realloc(state->generators, 
        state->generator_capacity * sizeof(unsigned int*));
}

state->generators[state->num_generators] = automorphism;  // NO NULL CHECK AFTER REALLOC
state->num_generators++;
```

**Why it crashes:** If `bliss_realloc()` fails, `state->generators` becomes NULL, but the code proceeds to dereference it.

## Refinement Algorithm Issues

### 8. Infinite Refinement Loop (`bliss_refinement.c` Lines 180-200)

**The Bug:**
```c
bool refine_partition_complete(bliss_graph_t *graph, partition_t *partition) {
    bool global_changed = false;
    bool local_changed = true;
    unsigned int round = 0;

    while (local_changed && round < MAX_REFINEMENT_ROUNDS) {
        local_changed = refine_partition_one_round(graph, partition);
        
        if (local_changed) {
            global_changed = true;
            // Canonicalize ordering in all cells after splitting
            unsigned int current_num_cells = partition->num_cells;
            for (unsigned int i = 0; i < current_num_cells; i++) {
                canonicalize_cell_ordering(partition, i, graph);  // CAN MODIFY partition->num_cells
            }
        }
        round++;
    }
}
```

**Why it crashes:** `canonicalize_cell_ordering()` can modify `partition->num_cells` during iteration, causing the loop to access invalid cell indices.

## Specific Crash Scenarios

### Scenario 1: Petersen Graph Test
1. Create Petersen graph (10 vertices, complex symmetry)
2. Initial partition refinement works fine
3. First individualization creates corrupted child partition
4. Recursive search with corrupted partition accesses invalid memory
5. **Segmentation fault**

### Scenario 2: Complete Graph Test  
1. Create K5 (complete graph, 120 automorphisms)
2. Deep search tree due to high symmetry
3. Multiple partition copying operations
4. Memory corruption accumulates
5. Generator storage overflow
6. **Heap corruption crash**

### Scenario 3: Random Graph Test
1. Create random graph with specific seed
2. Refinement algorithm hits edge case
3. Infinite refinement loop due to cell modification during iteration
4. Eventually hits MAX_REFINEMENT_ROUNDS
5. Returns invalid partition state
6. **Assertion failure or infinite loop**

## Fix Priority Order

### Critical (Must Fix for Basic Functionality)
1. **Fix partition memory management** in `create_child_node_unified()`
2. **Add error checking** to partition operations
3. **Fix refinement loop** cell modification issue

### High (Needed for Advanced Tests)
4. **Fix cell splitting logic** in refinement
5. **Add partition validation** in search loop
6. **Fix generator storage** overflow handling

### Medium (Robustness)
7. **Improve error propagation** throughout
8. **Add memory leak prevention**
9. **Optimize for complex graphs**

The crashes in advanced tests are primarily due to the memory corruption cascade starting from the partition copying bug, amplified by the lack of validation and error checking throughout the algorithm.