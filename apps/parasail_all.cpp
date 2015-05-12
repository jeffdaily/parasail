/**
 * @file parasail_all
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 *
 * Reads packed fasta file of N sequences.
 * Indexes input to learn length and end locations of each sequence.
 * Creates SA, LCP, and BWT.
 * Runs maximal pairs algorithm with the given minimum cutoff length.
 * For each sequence pair, performs semi-global alignment.
 *
 * Note about the input file. It is expected to be a packed fasta file
 * with each sequence delimited by the '$' sentinal. For example,
 * "banana$mississippi$foo$bar$".
 */
#include "config.h"

#include <unistd.h>

#include <cctype>
#include <cfloat>
#include <climits>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <set>
#include <stack>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "parasail.h"
#include "parasail/memory.h"
#include "parasail/function_lookup.h"
#include "parasail/matrix_lookup.h"

#include "sais.h"

using ::std::make_pair;
using ::std::pair;
using ::std::set;
using ::std::size_t;
using ::std::stack;
using ::std::vector;

typedef pair<int,int> Pair;

typedef set<Pair> PairSet;

struct quad {
    int lcp;
    int lb;
    int rb;
    vector<quad> children;

    quad()
        : lcp(0), lb(0), rb(INT_MAX), children() {}
    quad(int lcp, int lb, int rb)
        : lcp(lcp), lb(lb), rb(rb), children() {}
    quad(int lcp, int lb, int rb, vector<quad> children)
        : lcp(lcp), lb(lb), rb(rb), children(children) {}

    bool empty() { return rb == INT_MAX; }
};

inline static void pair_check(
        unsigned long &count_generated,
        PairSet &pairs,
        const int &i,
        const int &j,
        const int * const restrict SA,
        const unsigned char * const restrict BWT,
        const int * const restrict SID,
        const char &sentinal);

inline static void process(
        unsigned long &count_generated,
        PairSet &pairs,
        const quad &q,
        const int * const restrict SA,
        const unsigned char * const restrict BWT,
        const int * const restrict SID,
        const char &sentinal,
        const int &cutoff);

inline static int self_score(
        const char * const restrict seq,
        int len,
        const parasail_matrix_t *matrix);

static void print_help(const char *progname, int status) {
    fprintf(stderr, "\nusage: %s "
            "[-a funcname] "
            "[-c cutoff] "
            "[-e gap_extend] "
            "[-o gap_open] "
            "[-m matrix] "
            "[-l AOL] "
            "[-s SIM] "
            "[-i OS] "
            "-f file "
            "[-g output] "
            "\n\n",
            progname);
    fprintf(stderr, "Defaults:\n"
            "  funcname: sg_stats_striped_16\n"
            "    cutoff: 7, must be >= 1, exact match length cutoff\n"
            "gap_extend: 1, must be >= 0\n"
            "  gap_open: 10, must be >= 0\n"
            "    matrix: blosum62\n"
            "       AOL: 80, must be 0 <= AOL <= 100, percent alignment length\n"
            "       SIM: 40, must be 0 <= SIM <= 100, percent exact matches\n"
            "        OS: 30, must be 0 <= OS <= 100, percent optimal score over self score\n"
            "      file: no default, must be in FASTA format\n"
            "    output: edges.csv\n"
            );
    exit(status);
}

int main(int argc, char **argv) {
    FILE *fip = NULL;
    FILE *fop = NULL;
    const char *fname = NULL;
    const char *oname = "edges.csv";
    unsigned char *T = NULL;
#ifdef _OPENMP
    int num_threads = 1;
#endif
    int *SA = NULL;
    int *LCP = NULL;
    unsigned char *BWT = NULL;
    int *SID = NULL;
    vector<int> BEG;
    vector<int> END;
    long n = 0;
    double start = 0;
    double finish = 0;
    int i = 0;
    int longest = 0;
    int sid = 0;
    char sentinal = 0;
    int cutoff = 7;
    PairSet pairs;
    unsigned long count_possible = 0;
    unsigned long count_generated = 0;
    unsigned long work = 0;
    int c = 0;
    const char *funcname = "sg_stats_striped_16";
    parasail_function_t *function = NULL;
    const char *matrixname = "blosum62";
    const parasail_matrix_t *matrix = NULL;
    int gap_open = 10;
    int gap_extend = 1;
    const char *progname = "parasail_all";
    int AOL = 80;
    int SIM = 40;
    int OS = 30;

    /* Check arguments. */
    while ((c = getopt(argc, argv, "a:c:e:f:g:hm:o:l:s:i:")) != -1) {
        switch (c) {
            case 'a':
                funcname = optarg;
                break;
            case 'c':
                cutoff = atoi(optarg);
                if (cutoff <= 0) {
                    print_help(progname, EXIT_FAILURE);
                }
                break;
            case 'e':
                gap_extend = atoi(optarg);
                if (gap_extend < 0) {
                    print_help(progname, EXIT_FAILURE);
                }
                break;
            case 'f':
                fname = optarg;
                break;
            case 'g':
                oname = optarg;
                break;
            case 'h':
                print_help(progname, EXIT_FAILURE);
                break;
            case 'm':
                matrixname = optarg;
                break;
            case 'o':
                gap_open = atoi(optarg);
                if (gap_open < 0) {
                    print_help(progname, EXIT_FAILURE);
                }
                break;
            case 'l':
                AOL = atoi(optarg);
                if (AOL < 0 || AOL > 100) {
                    print_help(progname, EXIT_FAILURE);
                }
                break;
            case 's':
                SIM = atoi(optarg);
                if (SIM < 0 || SIM > 100) {
                    print_help(progname, EXIT_FAILURE);
                }
                break;
            case 'i':
                OS = atoi(optarg);
                if (OS < 0 || OS > 100) {
                    print_help(progname, EXIT_FAILURE);
                }
                break;
            case '?':
                if (optopt == 'a'
                        || optopt == 'c'
                        || optopt == 'e'
                        || optopt == 'f'
                        || optopt == 'g'
                        || optopt == 'm'
                        || optopt == 'o'
                        || optopt == 'l'
                        || optopt == 's'
                        || optopt == 'i'
                        ) {
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
                exit(EXIT_FAILURE);
            default:
                fprintf(stderr, "default case in getopt\n");
                exit(EXIT_FAILURE);
        }
    }

    /* select the function */
    if (funcname) {
        if (NULL == strstr(funcname, "stats") ) {
            fprintf(stderr, "Specified function does not calculate stats.\n");
            fprintf(stderr, "Use a 'stats' function, e.g, sg_stats_striped_16.\n");
            exit(EXIT_FAILURE);
        }
        function = parasail_lookup_function(funcname);
        if (NULL == function) {
            fprintf(stderr, "Specified function not found.\n");
            exit(EXIT_FAILURE);
        }
    }
    else {
        fprintf(stderr, "No alignment function specified.\n");
        exit(EXIT_FAILURE);
    }

    /* select the substitution matrix */
    if (matrixname) {
        matrix = parasail_matrix_lookup(matrixname);
        if (NULL == matrix) {
            fprintf(stderr, "Specified substitution matrix not found.\n");
            exit(EXIT_FAILURE);
        }
    }

    if (fname == NULL) {
        fprintf(stderr, "missing input file\n");
        print_help(progname, EXIT_FAILURE);
    }

    /* Best to know early whether we can open the output file. */
    if((fop = fopen(oname, "w")) == NULL) {
        fprintf(stderr, "%s: Cannot open output file `%s': ", progname, oname);
        perror("fopen");
        exit(EXIT_FAILURE);
    }
    
    /* Open a file for reading. */
    if((fip = fopen(fname, "r")) == NULL) {
        fprintf(stderr, "%s: Cannot open input file `%s': ", progname, fname);
        perror("fopen");
        exit(EXIT_FAILURE);
    }

    /* Get the file size. */
    if(fseek(fip, 0, SEEK_END) == 0) {
        n = ftell(fip);
        if(n < 0) {
            fprintf(stderr, "%s: Cannot ftell `%s': ", progname, fname);
            perror("ftell");
            exit(EXIT_FAILURE);
        }
        rewind(fip);
    } else {
        fprintf(stderr, "%s: Cannot fseek `%s': ", progname, fname);
        perror("fseek");
        exit(EXIT_FAILURE);
    }

    /* Allocate file buffer, read the entire file, then pack it. */
    T = (unsigned char *)malloc((size_t)(n+1) * sizeof(unsigned char));
    if(fread(T, sizeof(unsigned char), (size_t)n, fip) != (size_t)n) {
        fprintf(stderr, "%s: %s `%s': ",
                progname,
                (ferror(fip) || !feof(fip)) ? "Cannot read from" : "Unexpected EOF in",
                fname);
        perror("fread");
        exit(EXIT_FAILURE);
    }
    fclose(fip);
    T[n]='\0'; /* so we can print it */
    fprintf(stdout, "%20s: %s\n", "filename", fname);
    fprintf(stdout, "%20s: %ld bytes\n", "original size", n);
    /* Pack the buffer so it's ready for sais function. */
    {
        long i;
        char first = 1;
        long w = 0;
        long save = 0;
        for (i=0; i<n; ++i) { 
            if (T[i] == '>') {
                if (first) {
                    first = 0;
                }
                else {
                    T[w++] = '$';
                }
                /* skip rest of this line */
                while (T[i] != '\n') {
                    ++i;
                }
            }
            else if (isalpha(T[i])) {
                T[w++] = T[i];
            }
            else if (T[i] == '\n') {
                /* ignore newline */
            }
            else {
                fprintf(stderr, "uh oh at T[%ld]='%c'\n", i, T[i]);
                exit(EXIT_FAILURE);
            }
        }
        T[w++] = '$';
        save = w;
        /* nullifiy rest of buffer */
        while (w < n) {
            T[w++] = '\0';
        }
        /* new length */
        n = save;
    }
    fprintf(stdout, "%20s: %ld bytes\n", "packed size", n);

    /* Allocate memory. */
    SA = (int *)malloc((size_t)(n+1) * sizeof(int)); /* +1 for computing LCP */
    LCP = (int *)malloc((size_t)(n+1) * sizeof(int)); /* +1 for lcp tree */
    BWT = (unsigned char *)malloc((size_t)(n+1) * sizeof(unsigned char));
    SID = (int *)malloc((size_t)n * sizeof(int));
    if((T == NULL)
            || (SA == NULL)
            || (LCP == NULL)
            || (BWT == NULL)
            || (SID == NULL))
    {
        fprintf(stderr, "%s: Cannot allocate memory.\n", progname);
        exit(EXIT_FAILURE);
    }

    /* determine sentinal */
    if (sentinal == 0) {
        int off = 0;
        while (!isgraph(T[n-off])) {
            ++off;
        }
        sentinal = T[n-off];
    }

    /* determine actual end of file (last char) */
    {
        int off = 0;
        while (!isgraph(T[n-off])) {
            ++off;
        }
        n = n - off + 1;
    }

    /* scan T from left to build sequence ID and end index */
    longest = 0;
    sid = 0;
    BEG.push_back(0);
    for (i=0; i<n; ++i) {
        SID[i] = sid;
        if (T[i] == sentinal) {
            int len = i-BEG[sid];
            longest = (len>longest) ? len : longest;
            END.push_back(i);
            BEG.push_back(i+1);
            ++sid;
        }
    }
    longest += 1;
    if (0 == sid) { /* no sentinal found */
        fprintf(stderr, "no sentinal(%c) found in input\n", sentinal);
        exit(EXIT_FAILURE);
    }

    /* Construct the suffix and LCP arrays.
     * The following sais routine is from Fischer, with bugs fixed. */
    start = parasail_time();
    if(sais(T, SA, LCP, (int)n) != 0) {
        fprintf(stderr, "%s: Cannot allocate memory.\n", progname);
        exit(EXIT_FAILURE);
    }
    finish = parasail_time();
    fprintf(stdout, "%20s: %.4f seconds\n", "induced SA time", finish-start);

    /* construct naive BWT: */
    start = parasail_time();
    for (i = 0; i < n; ++i) {
        BWT[i] = (SA[i] > 0) ? T[SA[i]-1] : sentinal;
    }
    finish = parasail_time();
    fprintf(stdout, "%20s: %.4f seconds\n", "naive BWT time", finish-start);

    /* "fix" the LCP array to clamp LCP's that are too long */
    start = parasail_time();
    for (i = 0; i < n; ++i) {
        int len = END[SID[SA[i]]] - SA[i]; /* don't include sentinal */
        if (LCP[i] > len) LCP[i] = len;
    }
    finish = parasail_time();
    fprintf(stdout, "%20s: %.4f seconds\n", "clamp LCP time", finish-start);

    /* The GSA we create will put all sentinals either at the beginning
     * or end of the SA. We don't want to count all of the terminals,
     * nor do we want to process them in our bottom-up traversal. */
    /* do the sentinals appear at the beginning or end of SA? */
    int bup_start = 1;
    int bup_stop = n;
    if (T[SA[0]] == sentinal) {
        /* sentinals at beginning */
        bup_start = sid+1;
        bup_stop = n;
    }
    else if (T[SA[n-1]] == sentinal) {
        /* sentinals at end */
        bup_start = 1;
        bup_stop = n-sid;
    }
    else {
        fprintf(stderr, "sentinals not found at beginning or end of SA\n");
        exit(EXIT_FAILURE);
    }

    /* DFS of enhanced SA, from Abouelhoda et al */
    start = parasail_time();
    count_generated = 0;
    LCP[n] = 0; /* doesn't really exist, but for the root */
    {
        stack<quad> the_stack;
        quad last_interval;
        the_stack.push(quad());
        for (i = bup_start; i <= bup_stop; ++i) {
            int lb = i - 1;
            while (LCP[i] < the_stack.top().lcp) {
                the_stack.top().rb = i - 1;
                last_interval = the_stack.top();
                the_stack.pop();
                process(count_generated, pairs, last_interval, SA, BWT, SID, sentinal, cutoff);
                lb = last_interval.lb;
                if (LCP[i] <= the_stack.top().lcp) {
                    last_interval.children.clear();
                    the_stack.top().children.push_back(last_interval);
                    last_interval = quad();
                }
            }
            if (LCP[i] > the_stack.top().lcp) {
                if (!last_interval.empty()) {
                    last_interval.children.clear();
                    the_stack.push(quad(LCP[i],lb,INT_MAX,vector<quad>(1, last_interval)));
                    last_interval = quad();
                }
                else {
                    the_stack.push(quad(LCP[i],lb,INT_MAX));
                }
            }
        }
        the_stack.top().rb = bup_stop - 1;
        process(count_generated, pairs, the_stack.top(), SA, BWT, SID, sentinal, cutoff);
    }
    finish = parasail_time();
    count_possible = ((unsigned long)sid)*((unsigned long)sid-1)/2;
    fprintf(stdout, "%20s: %.4f seconds\n", "processing time", finish-start);
    fprintf(stdout, "%20s: %d\n", "number of sequences", sid);
    fprintf(stdout, "%20s: %lu\n", "possible pairs", count_possible);
    fprintf(stdout, "%20s: %lu\n", "generated pairs", count_generated);
    fprintf(stdout, "%20s: %zu\n", "unique pairs", pairs.size());

    /* Deallocate memory. */
    free(SA);
    free(LCP);
    free(BWT);
    free(SID);

#ifdef _OPENMP
    num_threads = omp_get_max_threads();
    fprintf(stdout, "%20s: %d\n", "omp num threads", num_threads);
#endif

    /* OpenMP can't iterate over an STL set. Convert to STL vector. */
    start = parasail_time();
    vector<Pair> vpairs(pairs.begin(), pairs.end());
    vector<pair<size_t,parasail_result_t*> > results;
    results.reserve(vpairs.size());
    finish = parasail_time();
    fprintf(stdout, "%20s: %.4f seconds\n", "openmp prep time", finish-start);

    /* align pairs */
    start = parasail_time();
#pragma omp parallel
    {
#pragma omp for schedule(dynamic) nowait
        for (size_t index=0; index<vpairs.size(); ++index) {
            int i = vpairs[index].first;
            int j = vpairs[index].second;
            int i_beg = BEG[i];
            int i_end = END[i];
            int i_len = i_end-i_beg;
            int j_beg = BEG[j];
            int j_end = END[j];
            int j_len = j_end-j_beg;
            unsigned long local_work = i_len * j_len;
            parasail_result_t *result = function(
                    (const char*)&T[i_beg], i_len,
                    (const char*)&T[j_beg], j_len,
                    gap_open, gap_extend, matrix);
#pragma omp critical
            {
                work += local_work;
                results.push_back(make_pair(index,result));
            }
        }
    }
    finish = parasail_time();
    fprintf(stdout, "%20s: %lu cells\n", "work", work);
    fprintf(stdout, "%20s: %.4f seconds\n", "alignment time", finish-start);
    fprintf(stdout, "%20s: %.4f \n", "gcups", double(work)/(finish-start)/1000000);

    /* Output results. */
    unsigned long edge_count = 0;
    for (size_t result_index=0; result_index<results.size(); ++result_index) {
        size_t index = results[result_index].first;
        parasail_result_t *result = results[result_index].second;
        int i = vpairs[index].first;
        int j = vpairs[index].second;
        int i_beg = BEG[i];
        int i_end = END[i];
        int i_len = i_end-i_beg;
        int j_beg = BEG[j];
        int j_end = END[j];
        int j_len = j_end-j_beg;
        int i_self_score = 0;
        int j_self_score = 0;
        int max_len = 0;
        int self_score_ = 0;

        if (result->score <= 0) {
            continue; /* skip negative scores */
        }

        i_self_score = self_score((const char*)&T[i_beg], i_len, matrix);
        j_self_score = self_score((const char*)&T[j_beg], j_len, matrix);

        if (i_len > j_len) {
            max_len = i_len;
            self_score_ = i_self_score;
        }
        else {
            max_len = j_len;
            self_score_ = j_self_score;
        }

        if ((result->length * 100 >= AOL * int(max_len))
                && (result->matches * 100 >= SIM * result->length)
                && (result->score * 100 >= OS * self_score_)) {
            ++edge_count;
            fprintf(fop, "%d,%d,%f,%f,%f\n",
                    i, j,
                    1.0*result->length/max_len,
                    1.0*result->matches/result->length,
                    1.0*result->score/self_score_);
        }
    }
    fclose(fop);

    fprintf(stdout, "%20s: %lu\n", "edges count", edge_count);

    /* Done with input text. */
    free(T);

    return 0;
}


inline static void pair_check(
        unsigned long &count_generated,
        PairSet &pairs,
        const int &i,
        const int &j,
        const int * const restrict SA,
        const unsigned char * const restrict BWT,
        const int * const restrict SID,
        const char &sentinal)
{
    const int &sidi = SID[SA[i]];
    const int &sidj = SID[SA[j]];
    if (BWT[i] != BWT[j] || BWT[i] == sentinal) {
        if (sidi != sidj) {
            ++count_generated;
            if (sidi < sidj) {
                pairs.insert(make_pair(sidi,sidj));
            }
            else {
                pairs.insert(make_pair(sidj,sidi));
            }
        }
    }
}

/* try to reduce number of duplicate pairs generated */
/* we observe that l-intervals (i.e. internal nodes) always have at
 * least two children, but these children could be singleton
 * l-intervals, e.g., [i..j]=[1..1], in addition to l-intervals with
 * non-singleton ranges/quads. For each l-interval, we take the cross
 * product of its child l-intervals. Naively, we could take the cross
 * product of the entire lb/rb range of the l-interval, but this
 * generates too many duplicate pairs. Instead, the complexity should be
 * bounded by the number of exact matches...
 */
inline static void process(
        unsigned long &count_generated,
        PairSet &pairs,
        const quad &q,
        const int * const restrict SA,
        const unsigned char * const restrict BWT,
        const int * const restrict SID,
        const char &sentinal,
        const int &cutoff)
{
    const int n_children = q.children.size();
    int child_index = 0;

    if (q.lcp < cutoff) return;

    if (n_children) {
        for (int i=q.lb; i<=q.rb; ++i) {
            int j = i+1;
            if (child_index < n_children) {
                if (i >= q.children[child_index].lb) {
                    j = q.children[child_index].rb+1;
                    if (i >= q.children[child_index].rb) {
                        ++child_index;
                    }
                }
            }
            for (/*nope*/; j<=q.rb; ++j) {
                pair_check(count_generated, pairs, i, j, SA, BWT, SID, sentinal);
            }
        }
    }
    else {
        for (int i=q.lb; i<=q.rb; ++i) {
            for (int j=i+1; j<=q.rb; ++j) {
                pair_check(count_generated, pairs, i, j, SA, BWT, SID, sentinal);
            }
        }
    }
}

inline static int self_score(
        const char * const restrict seq,
        int len,
        const parasail_matrix_t *matrix)
{
    int score = 0;
    for (int i=0; i<len; ++i) {
        unsigned char mapped = matrix->mapper[(unsigned char)seq[i]];
        score += matrix->matrix[matrix->size*mapped+mapped];
    }
    return score;
}

