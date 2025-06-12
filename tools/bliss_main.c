/*
 * bliss_main.c - Command line interface for bliss C library
 *
 * Provides a command-line tool similar to the original bliss executable
 * for computing automorphisms and canonical forms of graphs.
 */

#include "../bliss.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <time.h>

/* Command line options */
typedef struct {
    char *input_file;
    char *output_file;
    bool directed;
    bool canonical_labeling;
    bool automorphisms_only;
    bool quiet;
    bool verbose;
    bool statistics;
    bool write_dot;
    bliss_splitting_heuristic_t splitting_heuristic;
    bool use_component_recursion;
    bool use_failure_recording;
    bool use_long_prune;
    int verbose_level;
} options_t;

/* Global variables */
static unsigned int automorphism_count = 0;
static bool output_automorphisms = true;

/* Print usage information */
static void print_usage(const char *program_name) {
    printf("Usage: %s [OPTIONS] [INPUT_FILE]\n\n", program_name);
    printf("Compute automorphism groups and canonical forms of graphs.\n\n");
    printf("Options:\n");
    printf("  -d, --directed           Treat input as directed graph\n");
    printf("  -c, --canonical          Output canonical labeling\n");
    printf("  -a, --automorphisms      Output automorphisms only\n");
    printf("  -o, --output FILE        Output file (default: stdout)\n");
    printf("  -q, --quiet              Suppress automorphism output\n");
    printf("  -v, --verbose            Verbose output\n");
    printf("  -s, --statistics         Print search statistics\n");
    printf("  --dot                    Output in DOT format\n");
    printf("  --splitting HEURISTIC    Splitting heuristic (f, fl, fs, fm, flm, fsm)\n");
    printf("  --no-component-recursion Disable component recursion\n");
    printf("  --no-failure-recording   Disable failure recording\n");
    printf("  --no-long-prune          Disable long prune\n");
    printf("  --verbose-level LEVEL    Set verbose level (0-3)\n");
    printf("  -h, --help               Show this help message\n");
    printf("  --version                Show version information\n\n");
    printf("Input format: DIMACS graph format\n");
    printf("If no input file is specified, reads from stdin.\n\n");
    printf("Examples:\n");
    printf("  %s graph.dimacs                    # Find automorphisms\n", program_name);
    printf("  %s -c graph.dimacs                 # Find canonical labeling\n", program_name);
    printf("  %s -s --splitting fl graph.dimacs  # Use FL heuristic with stats\n", program_name);
    printf("  %s --dot -o graph.dot graph.dimacs # Convert to DOT format\n", program_name);
}

/* Print version information */
static void print_version(void) {
    printf("bliss %s\n", bliss_version_string());
    printf("Graph automorphism and canonical labeling tool\n");
    printf("C implementation based on the original bliss library\n");
    printf("Copyright (C) 2025\n");
}

/* Parse splitting heuristic string */
static bliss_splitting_heuristic_t parse_splitting_heuristic(const char *str) {
    if (strcmp(str, "f") == 0) return BLISS_SH_F;
    if (strcmp(str, "fl") == 0) return BLISS_SH_FL;
    if (strcmp(str, "fs") == 0) return BLISS_SH_FS;
    if (strcmp(str, "fm") == 0) return BLISS_SH_FM;
    if (strcmp(str, "flm") == 0) return BLISS_SH_FLM;
    if (strcmp(str, "fsm") == 0) return BLISS_SH_FSM;

    fprintf(stderr, "Error: Unknown splitting heuristic '%s'\n", str);
    fprintf(stderr, "Valid options: f, fl, fs, fm, flm, fsm\n");
    exit(1);
}

/* Automorphism output callback */
static void print_automorphism(void *user_param, unsigned int n, const unsigned int *aut) {
    FILE *output = (FILE *)user_param;

    if (!output_automorphisms) {
        automorphism_count++;
        return;
    }

    fprintf(output, "(");
    for (unsigned int i = 0; i < n; i++) {
        if (i > 0) fprintf(output, " ");
        fprintf(output, "%u", aut[i] + 1); /* Convert to 1-based indexing */
    }
    fprintf(output, ")\n");

    automorphism_count++;
}

/* Print canonical labeling */
static void print_canonical_labeling(FILE *output, const unsigned int *labeling, unsigned int n) {
    fprintf(output, "Canonical labeling: (");
    for (unsigned int i = 0; i < n; i++) {
        if (i > 0) fprintf(output, " ");
        fprintf(output, "%u", labeling[i] + 1); /* Convert to 1-based indexing */
    }
    fprintf(output, ")\n");
}

/* Print search statistics */
static void print_statistics(FILE *output, const bliss_stats_t *stats) {
    fprintf(output, "\nSearch statistics:\n");
    fprintf(output, "  Search tree nodes: %lu\n", bliss_stats_get_nof_nodes(stats));
    fprintf(output, "  Leaf nodes: %lu\n", bliss_stats_get_nof_leaf_nodes(stats));
    fprintf(output, "  Generators found: %lu\n", bliss_stats_get_nof_generators(stats));
    fprintf(output, "  Automorphisms: %u\n", automorphism_count);

    const char *group_size_str = bliss_stats_get_group_size_string(stats);
    if (group_size_str) {
        fprintf(output, "  Group size: %s\n", group_size_str);
    } else {
        fprintf(output, "  Group size (approx): %.2e\n", bliss_stats_get_group_size_approx(stats));
    }
}

/* Main program */
int main(int argc, char *argv[]) {
    /* Default options */
    options_t opts = {
        .input_file = NULL,
        .output_file = NULL,
        .directed = false,
        .canonical_labeling = false,
        .automorphisms_only = false,
        .quiet = false,
        .verbose = false,
        .statistics = false,
        .write_dot = false,
        .splitting_heuristic = BLISS_SH_F,
        .use_component_recursion = true,
        .use_failure_recording = true,
        .use_long_prune = true,
        .verbose_level = 0
    };

    /* Long options */
    static struct option long_options[] = {
        {"directed", no_argument, 0, 'd'},
        {"canonical", no_argument, 0, 'c'},
        {"automorphisms", no_argument, 0, 'a'},
        {"output", required_argument, 0, 'o'},
        {"quiet", no_argument, 0, 'q'},
        {"verbose", no_argument, 0, 'v'},
        {"statistics", no_argument, 0, 's'},
        {"dot", no_argument, 0, 1001},
        {"splitting", required_argument, 0, 1002},
        {"no-component-recursion", no_argument, 0, 1003},
        {"no-failure-recording", no_argument, 0, 1004},
        {"no-long-prune", no_argument, 0, 1005},
        {"verbose-level", required_argument, 0, 1006},
        {"help", no_argument, 0, 'h'},
        {"version", no_argument, 0, 1007},
        {0, 0, 0, 0}
    };

    /* Parse command line options */
    int c;
    while ((c = getopt_long(argc, argv, "dcao:qvsh", long_options, NULL)) != -1) {
        switch (c) {
        case 'd':
            opts.directed = true;
            break;
        case 'c':
            opts.canonical_labeling = true;
            break;
        case 'a':
            opts.automorphisms_only = true;
            break;
        case 'o':
            opts.output_file = optarg;
            break;
        case 'q':
            opts.quiet = true;
            break;
        case 'v':
            opts.verbose = true;
            break;
        case 's':
            opts.statistics = true;
            break;
        case 1001: /* --dot */
            opts.write_dot = true;
            break;
        case 1002: /* --splitting */
            opts.splitting_heuristic = parse_splitting_heuristic(optarg);
            break;
        case 1003: /* --no-component-recursion */
            opts.use_component_recursion = false;
            break;
        case 1004: /* --no-failure-recording */
            opts.use_failure_recording = false;
            break;
        case 1005: /* --no-long-prune */
            opts.use_long_prune = false;
            break;
        case 1006: /* --verbose-level */
            opts.verbose_level = atoi(optarg);
            break;
        case 'h':
            print_usage(argv[0]);
            return 0;
        case 1007: /* --version */
            print_version();
            return 0;
        case '?':
            return 1;
        default:
            fprintf(stderr, "Internal error: unhandled option %c\n", c);
            return 1;
        }
    }

    /* Get input file from remaining arguments */
    if (optind < argc) {
        opts.input_file = argv[optind];
    }

    /* Set verbose output options */
    if (opts.verbose) {
        bliss_set_verbose_level(opts.verbose_level > 0 ? opts.verbose_level : 1);
    }

    if (opts.quiet) {
        output_automorphisms = false;
    }

    /* Open input file */
    FILE *input = stdin;
    if (opts.input_file) {
        input = fopen(opts.input_file, "r");
        if (!input) {
            fprintf(stderr, "Error: Cannot open input file '%s'\n", opts.input_file);
            return 1;
        }
    }

    /* Open output file */
    FILE *output = stdout;
    if (opts.output_file) {
        output = fopen(opts.output_file, "w");
        if (!output) {
            fprintf(stderr, "Error: Cannot open output file '%s'\n", opts.output_file);
            if (input != stdin) fclose(input);
            return 1;
        }
    }

    /* Set verbose output file */
    if (opts.verbose) {
        bliss_set_verbose_file(stderr);
    }

    /* Read graph from input */
    if (opts.verbose) {
        fprintf(stderr, "Reading graph from %s...\n",
                opts.input_file ? opts.input_file : "stdin");
    }

    bliss_graph_t *graph = bliss_read_dimacs(input);
    if (input != stdin) fclose(input);

    if (!graph) {
        fprintf(stderr, "Error: Failed to read graph from input\n");
        if (output != stdout) fclose(output);
        return 1;
    }

    /* Override direction if specified */
    if (opts.directed && !bliss_is_directed(graph)) {
        fprintf(stderr, "Warning: Graph read as undirected but --directed specified\n");
    }

    /* Configure graph options */
    bliss_set_splitting_heuristic(graph, opts.splitting_heuristic);
    bliss_set_component_recursion(graph, opts.use_component_recursion);
    bliss_set_failure_recording(graph, opts.use_failure_recording);
    bliss_set_long_prune(graph, opts.use_long_prune);

    /* Print graph information */
    if (opts.verbose) {
        bliss_print_graph_info(graph, stderr);
    }

    /* Handle DOT output */
    if (opts.write_dot) {
        if (opts.verbose) {
            fprintf(stderr, "Writing graph in DOT format...\n");
        }
        bliss_write_dot(graph, output);
        if (output != stdout) fclose(output);
        bliss_release(graph);
        return 0;
    }

    /* Prepare for automorphism search */
    clock_t start_time = clock();
    bliss_stats_t *stats = bliss_stats_new();
    automorphism_count = 0;

    if (opts.verbose) {
        fprintf(stderr, "Starting automorphism search...\n");
    }

    /* Find automorphisms and/or canonical labeling */
    const unsigned int *canonical_labeling = NULL;

    if (opts.canonical_labeling) {
        canonical_labeling = bliss_find_canonical_labeling(graph, stats,
                                                          opts.quiet ? NULL : print_automorphism,
                                                          output);
    } else {
        bliss_find_automorphisms(graph, stats,
                                opts.quiet ? NULL : print_automorphism,
                                output);
    }

    clock_t end_time = clock();
    double elapsed_time = ((double)(end_time - start_time)) / CLOCKS_PER_SEC;

    /* Output results */
    if (!opts.quiet && !opts.automorphisms_only) {
        if (canonical_labeling) {
            print_canonical_labeling(output, canonical_labeling, bliss_get_nof_vertices(graph));
        }

        fprintf(output, "Automorphism group has %u generators\n", automorphism_count);
    }

    /* Print statistics */
    if (opts.statistics) {
        print_statistics(output, stats);
        fprintf(output, "Computation time: %.3f seconds\n", elapsed_time);
    }

    if (opts.verbose) {
        fprintf(stderr, "Search completed in %.3f seconds\n", elapsed_time);
        fprintf(stderr, "Found %u automorphisms\n", automorphism_count);
    }

    /* Cleanup */
    bliss_stats_release(stats);
    bliss_release(graph);
    if (output != stdout) fclose(output);

    return 0;
}