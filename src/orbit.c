#include "bliss.h"

void orbit_init(search_state_t *state, unsigned int n) {
    state->orbit = bliss_malloc(n * sizeof(unsigned int));
    state->orbit_size = bliss_malloc(n * sizeof(unsigned int));
    for (unsigned int i = 0; i < n; i++) {
        state->orbit[i] = i;
        state->orbit_size[i] = 1;
    }
}

static unsigned int orbit_find_internal(search_state_t *state, unsigned int v) {
    while (state->orbit[v] != v) {
        state->orbit[v] = state->orbit[state->orbit[v]];
        v = state->orbit[v];
    }
    return v;
}

unsigned int orbit_find(const search_state_t *state, unsigned int v) {
    return orbit_find_internal((search_state_t *)state, v);
}

void orbit_union(search_state_t *state, unsigned int a, unsigned int b) {
    unsigned int ra = orbit_find_internal(state, a);
    unsigned int rb = orbit_find_internal(state, b);
    if (ra == rb) return;
    if (state->orbit_size[ra] < state->orbit_size[rb]) {
        unsigned int tmp = ra;
        ra = rb;
        rb = tmp;
    }
    state->orbit[rb] = ra;
    state->orbit_size[ra] += state->orbit_size[rb];
}

void orbit_free(search_state_t *state) {
    bliss_free(state->orbit);
    bliss_free(state->orbit_size);
    state->orbit = NULL;
    state->orbit_size = NULL;
}
