#include <stdio.h>
#include <stdlib.h>

#include "parasail.h"
#include "parasail/io.h"

int main(int argc, char **argv)
{
    parasail_file_t *pf = NULL;
    char *packed = NULL;
    long packed_size = 0;
    parasail_file_labels_t *labels = NULL;
    unsigned long i;
    
    if (argc != 2) {
        fprintf(stderr, "program only takes file name as arg\n");
        exit(EXIT_FAILURE);
    }

    pf = parasail_open(argv[1]);

    packed = parasail_pack_with_labels(pf, &packed_size, &labels);

    for (i=0; i<labels->sequences; ++i) {
        printf("%s\n", labels->labels[i]);
    }

    free(packed);
    parasail_free_file_labels(labels);
    parasail_close(pf);

    return EXIT_SUCCESS;
}

