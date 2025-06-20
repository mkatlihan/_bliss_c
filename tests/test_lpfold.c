#include "bliss.h"
#include <stdio.h>
#include <math.h>
#include <float.h>

// Structure to represent an LP problem
typedef struct {
  int num_variables;
  int num_constraints;
  int objective_sense;  // 0 for MIN, 1 for MAX
  double objective_constant;
  double* objective_coefficients;  // size: num_variables
  double* constraint_rhs;          // size: num_constraints
  char* constraint_types;          // size: num_constraints, e.g., 'L', 'E', 'G'

  // Sparse matrix for constraints A*x (op) b
  // Represented by Compressed Sparse Column (CSC)
  int A_num_nonzeros;
  int* A_col_starts;      // size: num_variables + 1
  int* A_row_indices;     // size: A_num_nonzeros
  double* A_coeffs;       // size: A_num_nonzeros
  int* A_col_indices;     // size: A_num_nonzeros, computed from A_row_indices

  double* lower_bounds;   // size: num_variables
  double* upper_bounds;   // size: num_variables
  // Additional fields for partitioning
  int is_ptr; // 1 if this LPData is a pointer to another LPData, 0 if it's a copy
} LPData;

// Unified structure for both variable and constraint symmetries
typedef struct {
    int *variable_partition;      // partition[i] = group number for variable i
    int *constraint_partition;    // partition[i] = group number for constraint i
    int num_variable_groups;
    int num_constraint_groups;
    int *variable_representatives;    // one representative per variable group
    int *constraint_representatives;  // one representative per constraint group
    double variable_reduction_ratio;
    double constraint_reduction_ratio;
    double total_reduction_ratio;
} lp_symmetry_result_t;

// Coefficient analysis for creating equivalence classes
typedef struct {
    double *unique_coeffs;
    int num_coeff_classes;
    double tolerance;
} coefficient_analysis_t;

// Context for automorphism processing
typedef struct {
    const LPData *lp;
    lp_symmetry_result_t *result;
    double tolerance;
} automorphism_context_t;

unsigned int compute_variable_color(const LPData *lp, int var_idx, double tolerance);
unsigned int compute_constraint_color(const LPData *lp, int con_idx, double tolerance);
void merge_groups(int *partition, int size, int elem1, int elem2);
void finalize_partitions(lp_symmetry_result_t *result, int num_variables, int num_constraints);
// ===================================================================
// STEP 1: COEFFICIENT ANALYSIS
// ===================================================================

coefficient_analysis_t *analyze_coefficients(const LPData *lp, double tolerance) {
    DPRINTF("=== Analyzing Coefficient Patterns ===\n");
    
    coefficient_analysis_t *analysis = bliss_malloc(sizeof(coefficient_analysis_t));
    analysis->tolerance = tolerance;
    
    // Collect all unique coefficient values
    double *temp_coeffs = bliss_malloc(lp->A_num_nonzeros * sizeof(double));
    int num_unique = 0;
    
    for (int idx = 0; idx < lp->A_num_nonzeros; idx++) {
        double coeff = fabs(lp->A_coeffs[idx]); // Use absolute value for symmetry
        
        // Check if this coefficient is already in our list (within tolerance)
        bool found = false;
        for (int i = 0; i < num_unique; i++) {
            if (fabs(coeff - temp_coeffs[i]) <= tolerance) {
                found = true;
                break;
            }
        }
        
        if (!found) {
            temp_coeffs[num_unique++] = coeff;
        }
    }
    
    analysis->num_coeff_classes = num_unique;
    analysis->unique_coeffs = bliss_malloc(num_unique * sizeof(double));
    memcpy(analysis->unique_coeffs, temp_coeffs, num_unique * sizeof(double));
    
    DPRINTF("Found %d unique coefficient classes:\n", num_unique);
    for (int i = 0; i < num_unique && i < 10; i++) {
        DPRINTF("  Class %d: %.6f\n", i, analysis->unique_coeffs[i]);
    }
    if (num_unique > 10) {
        DPRINTF("  ... and %d more classes\n", num_unique - 10);
    }
    
    bliss_free(temp_coeffs);
    return analysis;
}

int get_coefficient_class(const coefficient_analysis_t *analysis, double coeff) {
    double abs_coeff = fabs(coeff);
    for (int i = 0; i < analysis->num_coeff_classes; i++) {
        if (fabs(abs_coeff - analysis->unique_coeffs[i]) <= analysis->tolerance) {
            return i;
        }
    }
    return 0; // Default to class 0 if not found
}

void free_coefficient_analysis(coefficient_analysis_t *analysis) {
    if (analysis) {
        bliss_free(analysis->unique_coeffs);
        bliss_free(analysis);
    }
}

// ===================================================================
// STEP 2: GRAPH CONSTRUCTION FOR LP SYMMETRY DETECTION
// ===================================================================

bliss_graph_t *create_lp_symmetry_graph(const LPData *lp, double coeff_tolerance) {
    DPRINTF("=== Creating LP Symmetry Detection Graph ===\n");
    DPRINTF("LP dimensions: %d variables, %d constraints, %d nonzeros\n", 
            lp->num_variables, lp->num_constraints, lp->A_num_nonzeros);
    
    // Analyze coefficient patterns
    coefficient_analysis_t *coeff_analysis = analyze_coefficients(lp, coeff_tolerance);
    
    // Graph structure:
    // - Variables: vertices 0 to (num_variables-1)
    // - Constraints: vertices num_variables to (num_variables + num_constraints - 1)
    // - Coefficient classes: vertices (num_variables + num_constraints) to end
    int var_offset = 0;
    int constraint_offset = lp->num_variables;
    int coeff_offset = lp->num_variables + lp->num_constraints;
    int total_vertices = lp->num_variables + lp->num_constraints + coeff_analysis->num_coeff_classes;
    
    DPRINTF("Graph layout: %d variables [0-%d], %d constraints [%d-%d], %d coeff classes [%d-%d]\n",
            lp->num_variables, lp->num_variables-1,
            lp->num_constraints, constraint_offset, constraint_offset + lp->num_constraints - 1,
            coeff_analysis->num_coeff_classes, coeff_offset, total_vertices - 1);
    
    bliss_graph_t *graph = bliss_new(total_vertices);
    
    // Color variables by their structural properties
    DPRINTF("Coloring variables...\n");
    for (int var = 0; var < lp->num_variables; var++) {
        unsigned int color = compute_variable_color(lp, var, coeff_tolerance);
        bliss_change_color(graph, var_offset + var, color);
    }
    
    // Color constraints by their structural properties
    DPRINTF("Coloring constraints...\n");
    for (int con = 0; con < lp->num_constraints; con++) {
        unsigned int color = compute_constraint_color(lp, con, coeff_tolerance);
        bliss_change_color(graph, constraint_offset + con, color);
    }
    
    // Color coefficient class vertices
    DPRINTF("Coloring coefficient classes...\n");
    for (int cc = 0; cc < coeff_analysis->num_coeff_classes; cc++) {
        unsigned int color = 1000000 + cc; // Unique high-value colors for coeff classes
        bliss_change_color(graph, coeff_offset + cc, color);
    }
    
    // Add structural edges: variable -- coeff_class -- constraint
    DPRINTF("Adding structural edges...\n");
    int edges_added = 0;
    
    for (int var = 0; var < lp->num_variables; var++) {
        int start = lp->A_col_starts[var];
        int end = lp->A_col_starts[var + 1];
        
        for (int idx = start; idx < end; idx++) {
            int con = lp->A_row_indices[idx];
            double coeff = lp->A_coeffs[idx];
            
            if (fabs(coeff) > coeff_tolerance) {
                int coeff_class = get_coefficient_class(coeff_analysis, coeff);
                
                // Create path: variable -> coeff_class -> constraint
                bliss_add_edge(graph, var_offset + var, coeff_offset + coeff_class);
                bliss_add_edge(graph, coeff_offset + coeff_class, constraint_offset + con);
                edges_added += 2;
            }
        }
    }
    
    DPRINTF("Added %d structural edges\n", edges_added);
    
    free_coefficient_analysis(coeff_analysis);
    return graph;
}

// ===================================================================
// STEP 3: VERTEX COLORING FUNCTIONS
// ===================================================================

unsigned int compute_variable_color(const LPData *lp, int var_idx, double tolerance) {
    unsigned int color = 0;
    
    // Objective coefficient class
    double obj_coeff = lp->objective_coefficients[var_idx];
    int obj_class = (int)(obj_coeff / tolerance) + 10000; // Shift to avoid negatives
    color = color * 100000 + (obj_class % 100000);
    
    // Lower bound class
    double lb = lp->lower_bounds[var_idx];
    if (lb <= -1e20) {
        color = color * 10 + 1; // Unbounded below
    } else {
        int lb_class = (int)(lb / tolerance) + 1000;
        color = color * 10000 + (lb_class % 10000);
    }
    
    // Upper bound class
    double ub = lp->upper_bounds[var_idx];
    if (ub >= 1e20) {
        color = color * 10 + 1; // Unbounded above
    } else {
        int ub_class = (int)(ub / tolerance) + 1000;
        color = color * 10000 + (ub_class % 10000);
    }
    
    // Number of constraints this variable appears in
    int start = lp->A_col_starts[var_idx];
    int end = lp->A_col_starts[var_idx + 1];
    int num_constraints = end - start;
    color = color * 1000 + (num_constraints % 1000);
    
    return color;
}

unsigned int compute_constraint_color(const LPData *lp, int con_idx, double tolerance) {
    unsigned int color = 0;
    
    // Constraint type
    char type = lp->constraint_types[con_idx];
    switch (type) {
        case 'L': color = 100; break;  // <=
        case 'E': color = 200; break;  // =
        case 'G': color = 300; break;  // >=
        default: color = 0; break;
    }
    
    // RHS class
    double rhs = lp->constraint_rhs[con_idx];
    int rhs_class = (int)(rhs / tolerance) + 10000;
    color = color * 100000 + (rhs_class % 100000);
    
    // Count variables in this constraint
    int num_vars = 0;
    for (int var = 0; var < lp->num_variables; var++) {
        int start = lp->A_col_starts[var];
        int end = lp->A_col_starts[var + 1];
        
        for (int idx = start; idx < end; idx++) {
            if (lp->A_row_indices[idx] == con_idx) {
                num_vars++;
                break;
            }
        }
    }
    color = color * 1000 + (num_vars % 1000);
    
    return color;
}

// ===================================================================
// STEP 4: AUTOMORPHISM PROCESSING
// ===================================================================

void process_lp_automorphism(void *user_param, unsigned int n, const unsigned int *aut) {
    automorphism_context_t *ctx = (automorphism_context_t*)user_param;
    const LPData *lp = ctx->lp;
    lp_symmetry_result_t *result = ctx->result;
    
    // Process variable-to-variable mappings
    for (int i = 0; i < lp->num_variables; i++) {
        int mapped_to = aut[i];
        
        if (mapped_to < lp->num_variables && i != mapped_to) {
            DPRINTF("Variable symmetry: %d <-> %d\n", i, mapped_to);
            merge_groups(result->variable_partition, lp->num_variables, i, mapped_to);
        }
    }
    
    // Process constraint-to-constraint mappings
    int constraint_offset = lp->num_variables;
    for (int i = 0; i < lp->num_constraints; i++) {
        int vertex_i = constraint_offset + i;
        int mapped_to = aut[vertex_i];
        
        if (mapped_to >= constraint_offset && 
            mapped_to < constraint_offset + lp->num_constraints && 
            vertex_i != mapped_to) {
            
            int con_j = mapped_to - constraint_offset;
            DPRINTF("Constraint symmetry: %d <-> %d\n", i, con_j);
            merge_groups(result->constraint_partition, lp->num_constraints, i, con_j);
        }
    }
}

// Union-find style group merging
void merge_groups(int *partition, int size, int elem1, int elem2) {
    int group1 = partition[elem1];
    int group2 = partition[elem2];
    
    if (group1 != group2) {
        int min_group = (group1 < group2) ? group1 : group2;
        int max_group = (group1 < group2) ? group2 : group1;
        
        for (int i = 0; i < size; i++) {
            if (partition[i] == max_group) {
                partition[i] = min_group;
            }
        }
    }
}

// ===================================================================
// STEP 5: MAIN SYMMETRY DETECTION FUNCTION
// ===================================================================

lp_symmetry_result_t *find_lp_equitable_partitions(const LPData *lp, double tolerance) {
    DPRINTF("=== Finding LP Equitable Partitions ===\n");
    DPRINTF("LP: %d variables, %d constraints, tolerance: %.2e\n", 
            lp->num_variables, lp->num_constraints, tolerance);
    
    // Create symmetry detection graph
    bliss_graph_t *graph = create_lp_symmetry_graph(lp, tolerance);
    
    // Initialize result structure
    lp_symmetry_result_t *result = bliss_malloc(sizeof(lp_symmetry_result_t));
    result->variable_partition = bliss_malloc(lp->num_variables * sizeof(int));
    result->constraint_partition = bliss_malloc(lp->num_constraints * sizeof(int));
    
    // Initialize partitions (each element in its own group)
    for (int i = 0; i < lp->num_variables; i++) {
        result->variable_partition[i] = i;
    }
    for (int i = 0; i < lp->num_constraints; i++) {
        result->constraint_partition[i] = i;
    }
    
    // Find automorphisms
    bliss_stats_t *stats = bliss_stats_new();
    automorphism_context_t context = {
        .lp = lp,
        .result = result,
        .tolerance = tolerance
    };
    
    DPRINTF("Running automorphism detection...\n");
    bliss_find_automorphisms(graph, stats, process_lp_automorphism, &context);
    
    DPRINTF("Automorphism search complete: %lu generators, %lu nodes\n",
            bliss_stats_get_nof_generators(stats), bliss_stats_get_nof_nodes(stats));
    
    // Finalize partitions
    finalize_partitions(result, lp->num_variables, lp->num_constraints);
    
    DPRINTF("Final results:\n");
    DPRINTF("  Variable groups: %d (%.2f%% reduction)\n", 
            result->num_variable_groups,
            100.0 * (lp->num_variables - result->num_variable_groups) / lp->num_variables);
    DPRINTF("  Constraint groups: %d (%.2f%% reduction)\n", 
            result->num_constraint_groups,
            100.0 * (lp->num_constraints - result->num_constraint_groups) / lp->num_constraints);
    
    bliss_stats_release(stats);
    bliss_release(graph);
    
    return result;
}

// ===================================================================
// STEP 6: PARTITION FINALIZATION
// ===================================================================

void finalize_partitions(lp_symmetry_result_t *result, int num_variables, int num_constraints) {
    DPRINTF("=== Finalizing Partitions ===\n");
    
    // Compact variable groups
    int *var_group_map = bliss_malloc(num_variables * sizeof(int));
    memset(var_group_map, -1, num_variables * sizeof(int));
    
    int next_var_group = 0;
    for (int i = 0; i < num_variables; i++) {
        int old_group = result->variable_partition[i];
        if (var_group_map[old_group] == -1) {
            var_group_map[old_group] = next_var_group++;
        }
        result->variable_partition[i] = var_group_map[old_group];
    }
    
    result->num_variable_groups = next_var_group;
    
    // Compact constraint groups
    int *con_group_map = bliss_malloc(num_constraints * sizeof(int));
    memset(con_group_map, -1, num_constraints * sizeof(int));
    
    int next_con_group = 0;
    for (int i = 0; i < num_constraints; i++) {
        int old_group = result->constraint_partition[i];
        if (con_group_map[old_group] == -1) {
            con_group_map[old_group] = next_con_group++;
        }
        result->constraint_partition[i] = con_group_map[old_group];
    }
    
    result->num_constraint_groups = next_con_group;
    
    // Select representatives
    result->variable_representatives = bliss_malloc(result->num_variable_groups * sizeof(int));
    result->constraint_representatives = bliss_malloc(result->num_constraint_groups * sizeof(int));
    
    memset(result->variable_representatives, -1, result->num_variable_groups * sizeof(int));
    memset(result->constraint_representatives, -1, result->num_constraint_groups * sizeof(int));
    
    for (int i = 0; i < num_variables; i++) {
        int group = result->variable_partition[i];
        if (result->variable_representatives[group] == -1) {
            result->variable_representatives[group] = i;
        }
    }
    
    for (int i = 0; i < num_constraints; i++) {
        int group = result->constraint_partition[i];
        if (result->constraint_representatives[group] == -1) {
            result->constraint_representatives[group] = i;
        }
    }
    
    // Calculate reduction ratios
    result->variable_reduction_ratio = (double)result->num_variable_groups / num_variables;
    result->constraint_reduction_ratio = (double)result->num_constraint_groups / num_constraints;
    result->total_reduction_ratio = result->variable_reduction_ratio * result->constraint_reduction_ratio;
    
    DPRINTF("Variable groups finalized: %d groups, representatives: ", result->num_variable_groups);
    for (int i = 0; i < result->num_variable_groups && i < 10; i++) {
        DPRINTF("x_%d ", result->variable_representatives[i]);
    }
    if (result->num_variable_groups > 10) DPRINTF("...");
    DPRINTF("\n");
    
    DPRINTF("Constraint groups finalized: %d groups, representatives: ", result->num_constraint_groups);
    for (int i = 0; i < result->num_constraint_groups && i < 10; i++) {
        DPRINTF("c_%d ", result->constraint_representatives[i]);
    }
    if (result->num_constraint_groups > 10) DPRINTF("...");
    DPRINTF("\n");
    
    bliss_free(var_group_map);
    bliss_free(con_group_map);
}

// ===================================================================
// CLEANUP FUNCTIONS
// ===================================================================

void free_lp_symmetry_result(lp_symmetry_result_t *result) {
    if (result) {
        bliss_free(result->variable_partition);
        bliss_free(result->constraint_partition);
        bliss_free(result->variable_representatives);
        bliss_free(result->constraint_representatives);
        bliss_free(result);
    }
}

// Create doubly-reduced LP from both variable and constraint equitable partitions
LPData *create_reduced_lp(const LPData *original_lp, const lp_symmetry_result_t *symmetries) {
    DPRINTF("=== Creating Doubly-Reduced LP (Variables + Constraints) ===\n");
    DPRINTF("Original LP: %d variables, %d constraints, %d nonzeros\n", 
            original_lp->num_variables, original_lp->num_constraints, original_lp->A_num_nonzeros);
    DPRINTF("Variable groups: %d (reduction: %.2f%%)\n", 
            symmetries->num_variable_groups,
            100.0 * (original_lp->num_variables - symmetries->num_variable_groups) / original_lp->num_variables);
    DPRINTF("Constraint groups: %d (reduction: %.2f%%)\n", 
            symmetries->num_constraint_groups,
            100.0 * (original_lp->num_constraints - symmetries->num_constraint_groups) / original_lp->num_constraints);
    
    // Allocate reduced LP structure
    LPData *reduced_lp = (LPData*)bliss_malloc(sizeof(LPData));
    
    // Basic dimensions - both variables AND constraints are reduced
    reduced_lp->num_variables = symmetries->num_variable_groups;
    reduced_lp->num_constraints = symmetries->num_constraint_groups;
    reduced_lp->objective_sense = original_lp->objective_sense;
    reduced_lp->objective_constant = original_lp->objective_constant;
    reduced_lp->is_ptr = 0;
    
    DPRINTF("Reduced LP dimensions: %d variables (was %d), %d constraints (was %d)\n", 
            reduced_lp->num_variables, original_lp->num_variables,
            reduced_lp->num_constraints, original_lp->num_constraints);
    
    // Allocate arrays for reduced LP
    reduced_lp->objective_coefficients = (double*)bliss_malloc(reduced_lp->num_variables * sizeof(double));
    reduced_lp->lower_bounds = (double*)bliss_malloc(reduced_lp->num_variables * sizeof(double));
    reduced_lp->upper_bounds = (double*)bliss_malloc(reduced_lp->num_variables * sizeof(double));
    
    reduced_lp->constraint_rhs = (double*)bliss_malloc(reduced_lp->num_constraints * sizeof(double));
    reduced_lp->constraint_types = (char*)bliss_malloc(reduced_lp->num_constraints * sizeof(char));
    
    // STEP 1: Build reduced objective and variable bounds
    DPRINTF("=== STEP 1: Building reduced objective and variable bounds ===\n");
    for (int var_group = 0; var_group < symmetries->num_variable_groups; var_group++) {
        int rep_var = symmetries->variable_representatives[var_group];
        
        // Count variables in this group
        int group_size = 0;
        for (int var = 0; var < original_lp->num_variables; var++) {
            if (symmetries->variable_partition[var] == var_group) {
                group_size++;
            }
        }
        
        DPRINTF("Variable group %d: representative var_%d, size %d\n", var_group, rep_var, group_size);
        
        // Objective coefficient: sum of all variables in the group
        reduced_lp->objective_coefficients[var_group] = 
            original_lp->objective_coefficients[rep_var] * group_size;
        
        // Bounds: scaled by group size (assuming all variables in group have same bounds)
        reduced_lp->lower_bounds[var_group] = original_lp->lower_bounds[rep_var] * group_size;
        reduced_lp->upper_bounds[var_group] = original_lp->upper_bounds[rep_var] * group_size;
        
        DPRINTF("  Objective coeff: %.6f (%.6f * %d)\n", 
                reduced_lp->objective_coefficients[var_group], 
                original_lp->objective_coefficients[rep_var], group_size);
        DPRINTF("  Bounds: [%.6f, %.6f]\n", 
                reduced_lp->lower_bounds[var_group], reduced_lp->upper_bounds[var_group]);
    }
    
    // STEP 2: Build reduced constraint RHS and types
    DPRINTF("=== STEP 2: Building reduced constraint RHS and types ===\n");
    for (int con_group = 0; con_group < symmetries->num_constraint_groups; con_group++) {
        int rep_con = symmetries->constraint_representatives[con_group];
        
        // Count constraints in this group
        int group_size = 0;
        for (int con = 0; con < original_lp->num_constraints; con++) {
            if (symmetries->constraint_partition[con] == con_group) {
                group_size++;
            }
        }
        
        DPRINTF("Constraint group %d: representative con_%d, size %d\n", con_group, rep_con, group_size);
        
        // RHS: use representative's RHS (should be same for all in group)
        // NOTE: For constraint grouping, we typically don't scale RHS - they should be identical
        reduced_lp->constraint_rhs[con_group] = original_lp->constraint_rhs[rep_con];
        
        // Type: use representative's type (should be same for all in group)
        reduced_lp->constraint_types[con_group] = original_lp->constraint_types[rep_con];
        
        DPRINTF("  RHS: %.6f, Type: %c\n", 
                reduced_lp->constraint_rhs[con_group], reduced_lp->constraint_types[con_group]);
        
        // Verify that all constraints in group have same RHS and type
        for (int con = 0; con < original_lp->num_constraints; con++) {
            if (symmetries->constraint_partition[con] == con_group) {
                if (fabs(original_lp->constraint_rhs[con] - reduced_lp->constraint_rhs[con_group]) > 1e-10) {
                    DPRINTF("  WARNING: Constraint %d has different RHS %.6f vs %.6f\n", 
                            con, original_lp->constraint_rhs[con], reduced_lp->constraint_rhs[con_group]);
                }
                if (original_lp->constraint_types[con] != reduced_lp->constraint_types[con_group]) {
                    DPRINTF("  WARNING: Constraint %d has different type %c vs %c\n", 
                            con, original_lp->constraint_types[con], reduced_lp->constraint_types[con_group]);
                }
            }
        }
    }
    
    // STEP 3: Build reduced constraint matrix (most complex part)
    DPRINTF("=== STEP 3: Building reduced constraint matrix ===\n");
    
    // Create coefficient accumulator: [constraint_group][variable_group] -> coefficient
    double **coeff_matrix = (double**)bliss_malloc(reduced_lp->num_constraints * sizeof(double*));
    for (int i = 0; i < reduced_lp->num_constraints; i++) {
        coeff_matrix[i] = (double*)bliss_malloc(reduced_lp->num_variables * sizeof(double));
        memset(coeff_matrix[i], 0, reduced_lp->num_variables * sizeof(double));
    }
    
    DPRINTF("Accumulating coefficients from original matrix...\n");
    
    // Accumulate coefficients: for each original (constraint, variable) -> add to (constraint_group, variable_group)
    for (int orig_var = 0; orig_var < original_lp->num_variables; orig_var++) {
        int var_group = symmetries->variable_partition[orig_var];
        int start = original_lp->A_col_starts[orig_var];
        int end = original_lp->A_col_starts[orig_var + 1];
        
        for (int idx = start; idx < end; idx++) {
            int orig_con = original_lp->A_row_indices[idx];
            int con_group = symmetries->constraint_partition[orig_con];
            double coeff = original_lp->A_coeffs[idx];
            
            coeff_matrix[con_group][var_group] += coeff;
        }
    }
    
    // Count nonzeros in reduced matrix
    int reduced_nonzeros = 0;
    for (int con_group = 0; con_group < reduced_lp->num_constraints; con_group++) {
        for (int var_group = 0; var_group < reduced_lp->num_variables; var_group++) {
            if (fabs(coeff_matrix[con_group][var_group]) > 1e-12) {
                reduced_nonzeros++;
            }
        }
    }
    
    reduced_lp->A_num_nonzeros = reduced_nonzeros;
    
    DPRINTF("Reduced matrix: %d nonzeros (was %d), reduction: %.2f%%\n", 
            reduced_nonzeros, original_lp->A_num_nonzeros,
            100.0 * (original_lp->A_num_nonzeros - reduced_nonzeros) / original_lp->A_num_nonzeros);
    
    // Allocate CSC arrays for reduced matrix
    reduced_lp->A_col_starts = (int*)bliss_malloc((reduced_lp->num_variables + 1) * sizeof(int));
    reduced_lp->A_row_indices = (int*)bliss_malloc(reduced_nonzeros * sizeof(int));
    reduced_lp->A_coeffs = (double*)bliss_malloc(reduced_nonzeros * sizeof(double));
    reduced_lp->A_col_indices = (int*)bliss_malloc(reduced_nonzeros * sizeof(int));
    
    // Build CSC structure for reduced matrix
    DPRINTF("Building CSC structure for reduced matrix...\n");
    
    // First pass: count nonzeros per column
    int *col_counts = (int*)bliss_malloc(reduced_lp->num_variables * sizeof(int));
    memset(col_counts, 0, reduced_lp->num_variables * sizeof(int));
    
    for (int var_group = 0; var_group < reduced_lp->num_variables; var_group++) {
        for (int con_group = 0; con_group < reduced_lp->num_constraints; con_group++) {
            if (fabs(coeff_matrix[con_group][var_group]) > 1e-12) {
                col_counts[var_group]++;
            }
        }
    }
    
    // Build column starts
    reduced_lp->A_col_starts[0] = 0;
    for (int var_group = 0; var_group < reduced_lp->num_variables; var_group++) {
        reduced_lp->A_col_starts[var_group + 1] = reduced_lp->A_col_starts[var_group] + col_counts[var_group];
        DPRINTF("Column %d: %d nonzeros, starts at index %d\n", 
                var_group, col_counts[var_group], reduced_lp->A_col_starts[var_group]);
    }
    
    // Second pass: fill in the actual values
    int *current_pos = (int*)bliss_malloc(reduced_lp->num_variables * sizeof(int));
    memcpy(current_pos, reduced_lp->A_col_starts, reduced_lp->num_variables * sizeof(int));
    
    for (int var_group = 0; var_group < reduced_lp->num_variables; var_group++) {
        DPRINTF("Filling column %d (variable group %d, rep: var_%d):\n", 
                var_group, var_group, symmetries->variable_representatives[var_group]);
        
        int entries_added = 0;
        for (int con_group = 0; con_group < reduced_lp->num_constraints; con_group++) {
            double coeff = coeff_matrix[con_group][var_group];
            
            if (fabs(coeff) > 1e-12) {
                int pos = current_pos[var_group]++;
                reduced_lp->A_row_indices[pos] = con_group;
                reduced_lp->A_coeffs[pos] = coeff;
                reduced_lp->A_col_indices[pos] = var_group;
                
                if (entries_added < 5) {  // Log first few entries
                    DPRINTF("  A[%d,%d] = %.6f (constraint group %d, rep: con_%d)\n", 
                            con_group, var_group, coeff, con_group,
                            symmetries->constraint_representatives[con_group]);
                }
                entries_added++;
            }
        }
        
        if (entries_added > 5) {
            DPRINTF("  ... and %d more entries\n", entries_added - 5);
        }
        DPRINTF("  Total entries in column %d: %d\n", var_group, entries_added);
    }
    
    // STEP 4: Verify matrix construction
    DPRINTF("=== STEP 4: Verifying reduced matrix construction ===\n");
    bool matrix_valid = true;
    int total_entries = 0;
    
    for (int var_group = 0; var_group < reduced_lp->num_variables; var_group++) {
        int expected_count = col_counts[var_group];
        int actual_count = reduced_lp->A_col_starts[var_group + 1] - reduced_lp->A_col_starts[var_group];
        total_entries += actual_count;
        
        if (expected_count != actual_count) {
            DPRINTF("ERROR: Column %d count mismatch: expected %d, got %d\n", 
                    var_group, expected_count, actual_count);
            matrix_valid = false;
        }
    }
    
    if (total_entries != reduced_nonzeros) {
        DPRINTF("ERROR: Total entries mismatch: expected %d, got %d\n", 
                reduced_nonzeros, total_entries);
        matrix_valid = false;
    }
    
    if (matrix_valid) {
        DPRINTF("Matrix construction verified successfully\n");
    } else {
        DPRINTF("ERROR: Matrix construction has inconsistencies\n");
    }
    
    // STEP 5: Print comprehensive summary statistics
    DPRINTF("\n=== DOUBLY-REDUCED LP SUMMARY ===\n");
    DPRINTF("Original dimensions: %d variables × %d constraints\n", 
            original_lp->num_variables, original_lp->num_constraints);
    DPRINTF("Reduced dimensions:  %d variables × %d constraints\n", 
            reduced_lp->num_variables, reduced_lp->num_constraints);
    DPRINTF("Variable reduction:  %.2f%% (%d → %d)\n", 
            100.0 * (original_lp->num_variables - reduced_lp->num_variables) / original_lp->num_variables,
            original_lp->num_variables, reduced_lp->num_variables);
    DPRINTF("Constraint reduction: %.2f%% (%d → %d)\n", 
            100.0 * (original_lp->num_constraints - reduced_lp->num_constraints) / original_lp->num_constraints,
            original_lp->num_constraints, reduced_lp->num_constraints);
    
    double total_size_reduction = 100.0 * (1.0 - (double)(reduced_lp->num_variables * reduced_lp->num_constraints) / 
                                                 (original_lp->num_variables * original_lp->num_constraints));
    DPRINTF("Total problem size reduction: %.2f%%\n", total_size_reduction);
    
    DPRINTF("Matrix sparsity: %d → %d nonzeros (%.2f%% reduction)\n", 
            original_lp->A_num_nonzeros, reduced_lp->A_num_nonzeros,
            100.0 * (original_lp->A_num_nonzeros - reduced_lp->A_num_nonzeros) / original_lp->A_num_nonzeros);
    
    // Calculate density changes
    double orig_density = (double)original_lp->A_num_nonzeros / (original_lp->num_variables * original_lp->num_constraints);
    double reduced_density = (double)reduced_lp->A_num_nonzeros / (reduced_lp->num_variables * reduced_lp->num_constraints);
    DPRINTF("Matrix density: %.4f → %.4f (%.2f%% change)\n", 
            orig_density, reduced_density, 100.0 * (reduced_density - orig_density) / orig_density);
    
    // Objective function statistics
    double orig_obj_norm = 0.0, reduced_obj_norm = 0.0;
    for (int i = 0; i < original_lp->num_variables; i++) {
        orig_obj_norm += original_lp->objective_coefficients[i] * original_lp->objective_coefficients[i];
    }
    for (int i = 0; i < reduced_lp->num_variables; i++) {
        reduced_obj_norm += reduced_lp->objective_coefficients[i] * reduced_lp->objective_coefficients[i];
    }
    DPRINTF("Objective ||c||₂: %.6f → %.6f\n", sqrt(orig_obj_norm), sqrt(reduced_obj_norm));
    
    // Clean up temporary arrays
    for (int i = 0; i < reduced_lp->num_constraints; i++) {
        bliss_free(coeff_matrix[i]);
    }
    bliss_free(coeff_matrix);
    bliss_free(col_counts);
    bliss_free(current_pos);
    
    DPRINTF("=== DOUBLY-REDUCED LP CREATION COMPLETE ===\n\n");
    
    return reduced_lp;
}

// Enhanced solution expansion for both variable and constraint symmetries
typedef struct {
    double *original_solution;     // Solution for original variables
    double *original_duals;        // Dual solution for original constraints
} expanded_solution_t;

expanded_solution_t *expand_solution_complete(const double *reduced_solution, 
                                             const double *reduced_duals,
                                             const lp_symmetry_result_t *symmetries,
                                             int original_num_variables,
                                             int original_num_constraints) {
    DPRINTF("=== Expanding Complete Solution ===\n");
    DPRINTF("Expanding %d variable groups to %d original variables\n", 
            symmetries->num_variable_groups, original_num_variables);
    DPRINTF("Expanding %d constraint groups to %d original constraints\n", 
            symmetries->num_constraint_groups, original_num_constraints);
    
    expanded_solution_t *expanded = bliss_malloc(sizeof(expanded_solution_t));
    
    // Expand primal solution (variables)
    expanded->original_solution = (double*)bliss_malloc(original_num_variables * sizeof(double));
    for (int var = 0; var < original_num_variables; var++) {
        int var_group = symmetries->variable_partition[var];
        expanded->original_solution[var] = reduced_solution[var_group];
        
        if (var < 10) {
            DPRINTF("  x_%d = group_%d_value = %.6f\n", var, var_group, reduced_solution[var_group]);
        }
    }
    
    // Expand dual solution (constraints)
    if (reduced_duals) {
        expanded->original_duals = (double*)bliss_malloc(original_num_constraints * sizeof(double));
        for (int con = 0; con < original_num_constraints; con++) {
            int con_group = symmetries->constraint_partition[con];
            expanded->original_duals[con] = reduced_duals[con_group];
            
            if (con < 10) {
                DPRINTF("  π_%d = group_%d_dual = %.6f\n", con, con_group, reduced_duals[con_group]);
            }
        }
    } else {
        expanded->original_duals = NULL;
    }
    
    if (original_num_variables > 10) {
        DPRINTF("  ... and %d more variable mappings\n", original_num_variables - 10);
    }
    if (original_num_constraints > 10) {
        DPRINTF("  ... and %d more dual mappings\n", original_num_constraints - 10);
    }
    
    DPRINTF("=== Complete Solution Expansion Complete ===\n");
    
    return expanded;
}

void free_expanded_solution(expanded_solution_t *expanded) {
    if (expanded) {
        bliss_free(expanded->original_solution);
        bliss_free(expanded->original_duals);
        bliss_free(expanded);
    }
}

// Function to free memory allocated for LPData members
void free_lp_data_members(LPData* lp) {
  if (!lp) return;

  free(lp->objective_coefficients);
  free(lp->constraint_rhs);
  free(lp->constraint_types);
  free(lp->A_col_starts);
  free(lp->A_row_indices);
  free(lp->A_coeffs);
  free(lp->lower_bounds);
  free(lp->upper_bounds);

  lp->objective_coefficients = NULL;
  lp->constraint_rhs = NULL;
  lp->constraint_types = NULL;
  lp->A_col_starts = NULL;
  lp->A_row_indices = NULL;
  lp->A_coeffs = NULL;
  lp->lower_bounds = NULL;
  lp->upper_bounds = NULL;
}

/**
 * Creates an LPData instance from the matrix example provided in the Grohe et al paper
 *
 * @return Pointer to newly allocated LPData structure (caller must free)
 */
LPData* create_example_lp6() {
  LPData* lp = (LPData*)calloc(1, sizeof(LPData));
  if (!lp) return NULL;

  // Set dimensions
  lp->num_variables = 13;
  lp->num_constraints = 7;
  lp->objective_sense = 1; // MAX
  lp->objective_constant = 0.0;

  // Allocate memory for LP components
  lp->objective_coefficients = (double*)malloc(lp->num_variables * sizeof(double));
  lp->constraint_rhs = (double*)malloc(lp->num_constraints * sizeof(double));
  lp->constraint_types = (char*)malloc(lp->num_constraints * sizeof(char));
  lp->lower_bounds = (double*)malloc(lp->num_variables * sizeof(double));
  lp->upper_bounds = (double*)malloc(lp->num_variables * sizeof(double));

  if (!lp->objective_coefficients || !lp->constraint_rhs || !lp->constraint_types ||
    !lp->lower_bounds || !lp->upper_bounds) {
    free_lp_data_members(lp);
    free(lp);
    return NULL;
  }

  // Set objective coefficients (cost function) - the green row at the bottom
  double obj_coeffs[13] = {
      2., 2., 2., // First 3 variables
      3. / 2, 3. / 2, 3. / 2, 3. / 2, // Next 4 variables
      1., 1., // Next 2 variables
      1. / 2, 1. / 2, 1. / 2, 1. / 2  // Last 4 variables
  };
  memcpy(lp->objective_coefficients, obj_coeffs, lp->num_variables * sizeof(double));

  // Set RHS values - all 6 from the green column on the right
  for (int i = 0; i < lp->num_constraints; i++) {
    lp->constraint_rhs[i] = 1.0;
    lp->constraint_types[i] = 'E'; // Equality constraints
  }

  // Set bounds - assuming all variables are non-negative
  for (int j = 0; j < lp->num_variables; j++) {
    lp->lower_bounds[j] = 0.0;
    lp->upper_bounds[j] = DBL_MAX; // Unbounded above
  }

  // Count non-zeros in the constraint matrix
  int nnz = 0;
  double A[7][13] = {
    // First 3 rows
    {3.0, -1.0, 1.0, 0.25, 0.25, 0.25, 0.25, 0.0, 0.0, 3.0, -2.0, 0.5, 0.5},
    {-1.0, 1.0, 3.0, 0.25, 0.25, 0.25, 0.25, 0.0, 0.0, -2.0, 3.0, 0.5, 0.5},
    {1.0, 3.0, -1.0, 0.25, 0.25, 0.25, 0.25, 0.0, 0.0, 0.5, 0.5, 0.5, 0.5},

    // Next 4 rows
    {0.0, 1.0 / 3.0, 2.0 / 3.0, 0.0, 1.5, 0.0, 1.5, 2.0, 0.0, 1.0, 0.0, -1.0, 0.0},
    {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0, 1.5, 0.0, 1.5, 0.0, 2.0, 0.0, 0.0, 1.0, 0.0, -1.0},
    {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0, 0.0, 1.5, 0.0, 1.5, 0.0, 2.0, -1.0, 0.0, 1.0, 0.0},
    {2.0 / 3.0, 1.0 / 3.0, 0.0, 1.5, 0.0, 1.5, 0.0, 0.0, 2.0, 0.0, -1.0, 0.0, 1.0}
  };

  // Count non-zeros
  for (int i = 0; i < lp->num_constraints; i++) {
    for (int j = 0; j < lp->num_variables; j++) {
      if (A[i][j] != 0.0) {
        nnz++;
      }
    }
  }

  lp->A_num_nonzeros = nnz;

  // Allocate memory for CSC format
  lp->A_col_starts = (int*)malloc((lp->num_variables + 1) * sizeof(int));
  lp->A_row_indices = (int*)malloc(nnz * sizeof(int));
  lp->A_coeffs = (double*)malloc(nnz * sizeof(double));

  if (!lp->A_col_starts || !lp->A_row_indices || !lp->A_coeffs) {
    free_lp_data_members(lp);
    free(lp);
    return NULL;
  }

  // Fill CSC format
  int idx = 0;
  lp->A_col_starts[0] = 0;

  for (int j = 0; j < lp->num_variables; j++) {
    for (int i = 0; i < lp->num_constraints; i++) {
      if (A[i][j] != 0.0) {
        lp->A_row_indices[idx] = i;
        lp->A_coeffs[idx] = A[i][j];
        idx++;
      }
    }
    lp->A_col_starts[j + 1] = idx;
  }

  return lp;
}

// Helper function to print LPData (simplified)
void print_lp_data_summary(const char* title, const LPData* lp) {
  if (!lp) {
    printf("%s: LPData is NULL\n", title);
    return;
  }
  printf("--- %s ---\n", title);
  printf("Variables: %d, Constraints: %d\n", lp->num_variables, lp->num_constraints);
  if (lp->num_variables > 0 && lp->objective_coefficients) {
    printf("Obj Coeffs (first few): ");
    for (int i = 0; i < lp->num_variables && i < 5; ++i) printf("%.2f ", lp->objective_coefficients[i]);
    printf("\n");
  }
  if (lp->num_constraints > 0 && lp->constraint_rhs) {
    printf("RHS (first few): ");
    for (int i = 0; i < lp->num_constraints && i < 5; ++i) printf("%.2f ", lp->constraint_rhs[i]);
    printf("\n");
  }
  printf("A non-zeros: %d\n", lp->A_num_nonzeros);
  if (lp->A_num_nonzeros > 0 && lp->A_col_starts && lp->A_row_indices && lp->A_coeffs) {
    printf("A_col_starts (first few): ");
    for (int i = 0; i < lp->num_variables + 1 && i < 6; ++i) printf("%d ", lp->A_col_starts[i]);
    printf("\n");
    printf("A_row_indices (first few): ");
    for (int i = 0; i < lp->A_num_nonzeros && i < 10; ++i) printf("%d ", lp->A_row_indices[i]);
    printf("\n");
    printf("A_coeffs (first few): ");
    for (int i = 0; i < lp->A_num_nonzeros && i < 10; ++i) printf("%.2f ", lp->A_coeffs[i]);
    printf("\n");
  }
  printf("--------------------\n");
}

void demonstrate_lp_folding(LPData *lp) {
    double tolerance = 1e-6;  // Coefficient tolerance
    
    // Find equitable partitions
    lp_symmetry_result_t *symmetries = find_lp_equitable_partitions(lp, tolerance);
    
    printf("Original LP: %d variables, %d constraints\n", 
           lp->num_variables, lp->num_constraints);

    
    // TODO: Create reduced LP using the equitable partition
    LPData *reduced_lp = create_reduced_lp(lp, symmetries);
    
    free_lp_symmetry_result(symmetries);
}

int main(void) {
    LPData *example_lp = create_example_lp6();
    if (!example_lp) {
        fprintf(stderr, "Failed to create example LP\n");
        return 1;
    }
    print_lp_data_summary("Example LP6", example_lp);
    demonstrate_lp_folding(example_lp);
    free_lp_data_members(example_lp);
    return 0;
}