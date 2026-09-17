/* Fixed logical work; dynamic chunks balance heterogeneous physical cores. */
#include <errno.h>
#include <inttypes.h>
#include <omp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

static uint64_t mix(uint64_t x) {
    x ^= x << 13;
    x ^= x >> 7;
    x ^= x << 17;
    return x;
}

int main(int argc, char **argv) {
    if (argc != 3) return 2;
    char *end;
    errno = 0;
    uint64_t iterations = strtoull(argv[1], &end, 10);
    if (errno || *end || !iterations || iterations > UINT64_C(10000000000) || argv[1][0] == '-') return 2;
    long workers = strtol(argv[2], &end, 10);
    if (errno || *end || workers < 1 || workers > 256) return 2;
    uint64_t chunk_size = UINT64_C(10000000);
    uint64_t chunks_per_worker = (iterations + chunk_size - 1) / chunk_size;
    uint64_t chunks = chunks_per_worker * (uint64_t)workers;
    uint64_t *results = calloc((size_t)chunks, sizeof(uint64_t));
    if (!results) return 3;
    int failed = 0;
    omp_set_dynamic(0);
    #pragma omp parallel num_threads(workers) reduction(|:failed)
    {
        size_t length = 262144;
        volatile uint64_t *arena = calloc(length, sizeof(uint64_t));
        if (!arena || omp_get_num_threads() != workers) {
            failed = 1;
        }
        #pragma omp for schedule(dynamic,1)
        for (uint64_t chunk = 0; chunk < chunks; ++chunk) {
            if (!arena) continue;
            uint64_t state = UINT64_C(20260920) + chunk;
            for (size_t j = 0; j < length; ++j) {
                state = mix(state);
                arena[j] = state;
            }
            uint64_t offset = (chunk % chunks_per_worker) * chunk_size;
            uint64_t count = iterations - offset;
            if (count > chunk_size) count = chunk_size;
            for (uint64_t j = 0; j < count; ++j) {
                size_t position = (size_t)(state & (length - 1));
                state = mix(state ^ arena[position] ^ j);
                arena[position] = state;
            }
            results[chunk] = state;
        }
        free((void *)arena);
    }
    if (failed) {
        free(results);
        return 3;
    }
    printf("{\"workers\":%ld,\"iterations_per_worker\":%" PRIu64 ",\"checksums\":[", workers, iterations);
    for (long i = 0; i < workers; ++i) {
        uint64_t checksum = 0;
        for (uint64_t j = 0; j < chunks_per_worker; ++j)
            checksum = mix(checksum ^ results[(uint64_t)i * chunks_per_worker + j]);
        printf("%s\"%016" PRIx64 "\"", i ? "," : "", checksum);
    }
    puts("]}");
    free(results);
    return 0;
}
