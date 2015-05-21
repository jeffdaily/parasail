/**
 * @file
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright (c) 2015 Battelle Memorial Institute.
 */
#include "config.h"

#include <stdlib.h>
#include <string.h>

#include "parasail.h"
#include "parasail/function_lookup.h"

parasail_function_t * parasail_lookup_function(const char *funcname)
{
    parasail_function_t * function = NULL;

    if (funcname) {
        int index = 0;
        parasail_function_info_t f;
        f = functions[index++];
        while (f.pointer) {
            if (0 == strcmp(funcname, f.name)) {
                function = f.pointer;
                break;
            }
            f = functions[index++];
        }
        if (!function) {
            /* perhaps caller forgot "parasail_" prefix? */
            const char *prefix = "parasail_";
            char *newname = (char*)malloc(strlen(prefix)+strlen(funcname)+1);
            strcpy(newname, prefix);
            strcat(newname, funcname);
            index = 0;
            f = functions[index++];
            while (f.pointer) {
                if (0 == strcmp(newname, f.name)) {
                    function = f.pointer;
                    break;
                }
                f = functions[index++];
            }
            free(newname);
        }
    }

    return function;
}

