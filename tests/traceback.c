#include "config.h"

#include <assert.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <errno.h>
#include <unistd.h>

#include "kseq.h"
KSEQ_INIT(int, read)

#include "parasail.h"
#include "parasail/memory.h"

static inline void parse_sequences(
        const char *filename,
        char ***strings_,
        unsigned long **sizes_,
        unsigned long *count_)
{
    FILE* fp;
    kseq_t *seq = NULL;
    int l = 0;
    char **strings = NULL;
    unsigned long *sizes = NULL;
    unsigned long count = 0;
    unsigned long memory = 1000;

    errno = 0;
    fp = fopen(filename, "r");
    if(fp == NULL) {
        perror("fopen");
        exit(1);
    }
    strings = malloc(sizeof(char*) * memory);
    sizes = malloc(sizeof(unsigned long) * memory);
    seq = kseq_init(fileno(fp));
    while ((l = kseq_read(seq)) >= 0) {
        errno = 0;
        strings[count] = strdup(seq->seq.s);
        if (NULL == strings[count]) {
            perror("strdup");
            exit(1);
        }
        sizes[count] = seq->seq.l;
        ++count;
        if (count >= memory) {
            char **new_strings = NULL;
            unsigned long *new_sizes = NULL;
            memory *= 2;
            errno = 0;
            new_strings = realloc(strings, sizeof(char*) * memory);
            if (NULL == new_strings) {
                perror("realloc");
                exit(1);
            }
            strings = new_strings;
            errno = 0;
            new_sizes = realloc(sizes, sizeof(unsigned long) * memory);
            if (NULL == new_sizes) {
                perror("realloc");
                exit(1);
            }
            sizes = new_sizes;
        }
    }
    kseq_destroy(seq);
    fclose(fp);

    *strings_ = strings;
    *sizes_ = sizes;
    *count_ = count;
}


int main(int argc, char **argv)
{
    long seqA_index = 0;
    long seqB_index = 1;
    const char *seqA = NULL;
    const char *seqB = NULL;
    int lena = 0;
    int lenb = 0;
    parasail_result_t *result = NULL;
    int c = 0;
    char *filename = NULL;
    char **sequences = NULL;
    unsigned long *sizes = NULL;
    unsigned long seq_count = 0;
    char *endptr = NULL;
    const char *funcname = "sw_trace";
    parasail_function_t *function = NULL;
    parasail_pfunction_t *pfunction = NULL;
    parasail_pcreator_t *pcreator = NULL;
    parasail_profile_t *profile = NULL;
    char *matrixname = NULL;
    const parasail_matrix_t *matrix = NULL;
    int open = 10;
    int extend = 1;
    int match = 1;
    int mismatch = 0;
    int use_dna = 0;

    while ((c = getopt(argc, argv, "a:x:y:f:m:o:e:dX:Y:")) != -1) {
        switch (c) {
            case 'a':
                funcname = optarg;
                break;
            case 'x':
                errno = 0;
                seqA_index = strtol(optarg, &endptr, 10);
                if (errno) {
                    perror("strtol seqA_index");
                    fprintf(stderr, "invalid seqA index\n");
                    exit(1);
                }
                break;
            case 'y':
                errno = 0;
                seqB_index = strtol(optarg, &endptr, 10);
                if (errno) {
                    perror("strtol seqB_index");
                    fprintf(stderr, "invalid seqB index\n");
                    exit(1);
                }
                break;
            case 'f':
                filename = optarg;
                break;
            case 'd':
                use_dna = 1;
                break;
            case 'X':
                match = atoi(optarg);
                if (match < 0) {
                    fprintf(stderr, "match must be >= 0\n");
                    exit(1);
                }
                break;
            case 'Y':
                mismatch = atoi(optarg);
                if (mismatch < 0) {
                    fprintf(stderr, "mismatch must be >= 0\n");
                    exit(1);
                }
                break;
            case 'm':
                matrixname = optarg;
                break;
            case 'o':
                errno = 0;
                open = strtol(optarg, &endptr, 10);
                if (errno) {
                    perror("strtol open");
                    exit(1);
                }
                break;
            case 'e':
                errno = 0;
                extend = strtol(optarg, &endptr, 10);
                if (errno) {
                    perror("strtol extend");
                    exit(1);
                }
                break;
            case '?':
                if (optopt == 'a'
                        || optopt == 'b'
                        || optopt == 'f'
                        || optopt == 'n') {
                    fprintf(stderr,
                            "Option -%c requires an argument.\n",
                            optopt);
                }
                else if (isprint(optopt)) {
                    fprintf(stderr, "Unknown option `-%c'.\n",
                            optopt);
                }
                else {
                    fprintf(stderr,
                            "Unknown option character `\\x%x'.\n",
                            optopt);
                }
                exit(1);
            default:
                fprintf(stderr, "default case in getopt\n");
                exit(1);
        }
    }

    if (filename) {
        parse_sequences(filename, &sequences, &sizes, &seq_count);
    }
    else {
        fprintf(stderr, "no filename specified\n");
        exit(1);
    }

    /* select the substitution matrix */
    if (NULL != matrixname && use_dna) {
        fprintf(stderr, "Cannot specify matrix name for DNA alignments.\n");
        exit(EXIT_FAILURE);
    }
    if (use_dna) {
        matrix = parasail_matrix_create("ACGT", match, -mismatch);
    }
    else {
        if (NULL == matrixname) {
            matrixname = "blosum62";
        }
        matrix = parasail_matrix_lookup(matrixname);
        if (NULL == matrix) {
            fprintf(stderr, "Specified substitution matrix not found.\n");
            exit(EXIT_FAILURE);
        }
    }


    if ((unsigned long)seqA_index >= seq_count) {
        fprintf(stderr, "seqA index out of bounds\n");
        exit(1);
    }

    if ((unsigned long)seqB_index >= seq_count) {
        fprintf(stderr, "seqB index out of bounds\n");
        exit(1);
    }

    /* make sure we were given a function that traces */
    if (NULL ==  strstr(funcname, "trace")) {
        fprintf(stderr, "function must enable backtraces\n");
        exit(1);
    }

    /* select the function */
    if (funcname) {
        if (NULL != strstr(funcname, "profile")) {
            pfunction = parasail_lookup_pfunction(funcname);
            if (NULL == pfunction) {
                fprintf(stderr, "Specified profile function not found.\n");
                exit(EXIT_FAILURE);
            }
            pcreator = parasail_lookup_pcreator(funcname);
            if (NULL == pcreator) {
                fprintf(stderr, "Specified profile creator function not found.\n");
                exit(EXIT_FAILURE);
            }

        }
        else {
            function = parasail_lookup_function(funcname);
            if (NULL == function) {
                fprintf(stderr, "Specified function not found.\n");
                exit(EXIT_FAILURE);
            }
        }
    }
    else {
        fprintf(stderr, "No alignment function specified.\n");
        exit(EXIT_FAILURE);
    }

    seqA = sequences[seqA_index];
    seqB = sequences[seqB_index];
    lena = strlen(seqA);
    lenb = strlen(seqB);

    printf("        file: %s\n", filename);
    printf("    function: %s\n", funcname);
    printf("      matrix: %s\n", matrixname);
    printf("    gap open: %d\n", open);
    printf("  gap extend: %d\n", extend);
    printf("    seq pair: %lu,%lu\n", seqA_index, seqB_index);
    printf("query length: %d\n", lena);
    printf("   db length: %d\n", lenb);

    if (pfunction) {
        profile = pcreator(seqA, lena, matrix);
        result = pfunction(profile, seqB, lenb, open, extend);
        parasail_profile_free(profile);
    }
    else {
        result = function(seqA, lena, seqB, lenb, open, extend, matrix);
    }
    assert(result->trace_table);
    assert(result->trace_ins_table);
    assert(result->trace_del_table);

    /* do the traceback */
    parasail_traceback(seqA, lena, seqB, lenb, matrix, result);

    printf("Score:         %d\n", result->score);
    printf("end_query:     %d\n", result->end_query);
    printf("end_ref:       %d\n", result->end_ref);

    parasail_result_free(result);

    /* let's see how the stats version of the same trace function
     * behaves */
    {
        char *dup = strdup(funcname);
        char *loc = strstr(dup, "trace");
        if (NULL != loc) {
            /* this works because lengths are same */
            (void)strcpy(loc, "stats");
        }
        printf("stats function '%s'\n", dup);
        function = parasail_lookup_function(dup);
        if (NULL == function) {
            fprintf(stderr, "corresponding stats function not found.\n");
            exit(EXIT_FAILURE);
        }
        
        result = function(seqA, lena, seqB, lenb, open, extend, matrix);
        printf("Length:        %d\n", result->length);
        printf("Identity:   %d/%d\n", result->matches, result->length);
        printf("Similarity: %d/%d\n", result->similar, result->length);
        printf("Gaps:       %d/%d\n", -1, result->length);
        printf("Score:         %d\n", result->score);
        printf("end_query:     %d\n", result->end_query);
        printf("end_ref:       %d\n", result->end_ref);

        parasail_result_free(result);

        free(dup);
    }
    /* cleanup parsed sequences */

    {
        unsigned long s;
        for (s=0; s<seq_count; ++s) {
            free(sequences[s]);
        }
    }
    free(sequences);
    free(sizes);

    return 0;
}

