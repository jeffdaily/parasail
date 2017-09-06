/**
 * @file parasail_query
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright 2012 Pacific Northwest National Laboratory. All rights reserved.
 *
 * Reads fasta file of database sequences.
 * Optionally reads fasta file of query sequences.
 * Optionally filters inputs through suffix array.
 * Indexes input to learn length and end locations of each sequence.
 * Creates SA, LCP, and BWT.
 * Runs maximal pairs algorithm with the given minimum cutoff length.
 * For each sequence pair, performs alignment of user choice.
 * Output is csv of alignments between sequences.
 */
#include "config.h"

#include <errno.h>
#include <sys/types.h>
#if defined(_MSC_VER)
#include "wingetopt/src/getopt.h"
#else
#include <pwd.h>
#include <unistd.h>
#endif

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

#ifdef USE_CILK
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include <cilk/reducer_opadd.h> 
#endif

#include "parasail.h"
#include "parasail/io.h"

#include "sais.h"

#if HAVE_VARIADIC_MACROS
#define eprintf(STREAM, ...) fprintf(STREAM, __VA_ARGS__); fflush(STREAM)
#else
#define eprintf fprintf
#endif

using ::std::bad_alloc;
using ::std::make_pair;
using ::std::pair;
using ::std::set;
using ::std::size_t;
using ::std::stack;
using ::std::vector;

typedef pair<int,int> Pair;

typedef set<Pair> PairSet;
typedef vector<Pair> PairVec;

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
        const vector<int> &DB,
        const char &sentinal);

inline static void process(
        unsigned long &count_generated,
        PairSet &pairs,
        const quad &q,
        const int * const restrict SA,
        const unsigned char * const restrict BWT,
        const int * const restrict SID,
        const vector<int> &DB,
        const char &sentinal,
        const int &cutoff);

inline static void print_array(
        const char * filename_,
        const int * const restrict array,
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len);

inline static int self_score(
        const char * const restrict seq,
        int len,
        const parasail_matrix_t *matrix);

inline static void output_edges(
        FILE *fop,
        bool has_query,
        long sid_crossover,
        unsigned char *T,
        int AOL,
        int SIM,
        int OS,
        const parasail_matrix_t *matrix,
        const vector<long> &BEG,
        const vector<long> &END,
        const PairVec &vpairs,
        const vector<parasail_result_t*> &results);

inline static void output_graph(
        FILE *fop,
        int which,
        unsigned long count,
        unsigned char *T,
        int AOL,
        int SIM,
        int OS,
        const parasail_matrix_t *matrix,
        const vector<long> &BEG,
        const vector<long> &END,
        const PairVec &vpairs,
        const vector<parasail_result_t*> &results);

inline static void output_stats(
        FILE *fop,
        bool has_query,
        long sid_crossover,
        const vector<long> &BEG,
        const vector<long> &END,
        const PairVec &vpairs,
        const vector<parasail_result_t*> &results);

inline static void output_basic(
        FILE *fop,
        bool has_query,
        long sid_crossover,
        const vector<long> &BEG,
        const vector<long> &END,
        const PairVec &vpairs,
        const vector<parasail_result_t*> &results);

inline static void output_tables(
        bool has_query,
        long sid_crossover,
        unsigned char *T,
        const vector<long> &BEG,
        const vector<long> &END,
        const PairVec &vpairs,
        const vector<parasail_result_t*> &results);

static void print_help(const char *progname, int status) {
    eprintf(stderr, "\nusage: %s "
            "[-a funcname] "
            "[-c cutoff] "
            "[-x] "
            "[-e gap_extend] "
            "[-o gap_open] "
            "[-m matrix] "
            "[-t threads] "
            "[-d] "
            "[-M match] "
            "[-X mismatch] "
            "[-k band size (for nw_banded)] "
            "[-l AOL] "
            "[-s SIM] "
            "[-i OS] "
            "-f file "
            "[-q query_file] "
            "[-g output_file] "
            "\n\n",
            progname);
    eprintf(stderr, "Defaults:\n"
            "   funcname: sw_stats_striped_16\n"
            "     cutoff: 7, must be >= 1, exact match length cutoff\n"
            "         -x: if present, don't use suffix array filter\n"
            " gap_extend: 1, must be >= 0\n"
            "   gap_open: 10, must be >= 0\n"
            "     matrix: blosum62\n"
            "         -d: if present, assume DNA alphabet\n"
            "      match: 1, must be >= 0\n"
            "   mismatch: 0, must be >= 0\n"
#ifdef _OPENMP
            "    threads: system-specific default, must be >= 1\n"
#else
            "    threads: Warning: ignored; OpenMP was not supported by your compiler\n"
#endif
            "        AOL: 80, must be 0 <= AOL <= 100, percent alignment length\n"
            "        SIM: 40, must be 0 <= SIM <= 100, percent exact matches\n"
            "         OS: 30, must be 0 <= OS <= 100, percent optimal score over self score\n"
            "       file: no default, must be in FASTA format\n"
            " query_file: no default, must be in FASTA format\n"
            "output_file: parasail.csv\n"
            );
    exit(status);
}

int main(int argc, char **argv) {
    FILE *fop = NULL;
    const char *fname = NULL;
    const char *qname = NULL;
    const char *oname = "parasail.csv";
    unsigned char *T = NULL;
    unsigned char *Q = NULL;
    int num_threads = -1;
    int *SA = NULL;
    int *LCP = NULL;
    unsigned char *BWT = NULL;
    int *SID = NULL;
    vector<long> BEG;
    vector<long> END;
    vector<int> DB;
    long n = 0;
    long t = 0;
    long q = 0;
    double start = 0;
    double finish = 0;
    long i = 0;
    long sid = 0;
    long sid_crossover = -1;
    char sentinal = 0;
    int cutoff = 7;
    bool use_filter = true;
    PairSet pairs;
    PairVec vpairs;
    unsigned long count_possible = 0;
    unsigned long count_generated = 0;
#ifdef USE_CILK
    cilk::reducer_opadd<unsigned long> work;
#else
    unsigned long work = 0;
#endif
    int c = 0;
    const char *funcname = "sw_stats_striped_16";
    parasail_function_t *function = NULL;
    parasail_pfunction_t *pfunction = NULL;
    parasail_pcreator_t *pcreator = NULL;
    int banded = 0;
    int kbandsize = 3;
    const char *matrixname = NULL;
    const parasail_matrix_t *matrix = NULL;
    int gap_open = 10;
    int gap_extend = 1;
    int match = 1;
    int mismatch = 0;
    bool use_dna = false;
    bool pairs_only = false;
    bool edge_output = false;
    bool graph_output = false;
    bool fpack = false;
    bool qpack = false;
    bool is_stats = true;
    bool is_table = false;
    bool has_query = false;
    const char *progname = "parasail_aligner";
    int AOL = 80;
    int SIM = 40;
    int OS = 30;

    /* Check arguments. */
    while ((c = getopt(argc, argv, "a:c:de:Ef:F:g:Ghi:k:l:m:M:o:pq:Q:s:t:xX:")) != -1) {
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
            case 'd':
                use_dna = true;
                break;
            case 'e':
                gap_extend = atoi(optarg);
                if (gap_extend < 0) {
                    print_help(progname, EXIT_FAILURE);
                }
                break;
            case 'E':
                edge_output = true;
                break;
            case 'f':
                fname = optarg;
                break;
            case 'F':
                fname = optarg;
                fpack = true;
                break;
            case 'g':
                oname = optarg;
                break;
            case 'G':
                graph_output = true;
                break;
            case 'h':
                print_help(progname, EXIT_FAILURE);
                break;
            case 'i':
                OS = atoi(optarg);
                if (OS < 0 || OS > 100) {
                    print_help(progname, EXIT_FAILURE);
                }
                break;
            case 'k':
                kbandsize = atoi(optarg);
                if (kbandsize <= 0) {
                    print_help(progname, EXIT_FAILURE);
                }
                break;
            case 'l':
                AOL = atoi(optarg);
                if (AOL < 0 || AOL > 100) {
                    print_help(progname, EXIT_FAILURE);
                }
                break;
            case 'm':
                matrixname = optarg;
                break;
            case 'M':
                match = atoi(optarg);
                if (match < 0) {
                    print_help(progname, EXIT_FAILURE);
                }
                break;
            case 'o':
                gap_open = atoi(optarg);
                if (gap_open < 0) {
                    print_help(progname, EXIT_FAILURE);
                }
                break;
            case 'p':
                pairs_only = true;
                break;
            case 'q':
                qname = optarg;
                break;
            case 'Q':
                qname = optarg;
                qpack = true;
                break;
            case 't':
                num_threads = atoi(optarg);
#ifdef _OPENMP
#else
                printf("-t number of threads requested, but OpenMP was not found during configuration. Running without threads.");
#endif
                break;
            case 'x':
                use_filter = false;
                break;
            case 'X':
                mismatch = atoi(optarg);
                if (mismatch < 0) {
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
                        || optopt == 'M'
                        || optopt == 'o'
                        || optopt == 'q'
                        || optopt == 'X'
                        || optopt == 'E'
                        || optopt == 'l'
                        || optopt == 's'
                        || optopt == 'i'
                        ) {
                    eprintf(stderr,
                            "Option -%c requires an argument.\n",
                            optopt);
                }
                else if (isprint(optopt)) {
                    eprintf(stderr, "Unknown option `-%c'.\n",
                            optopt);
                }
                else {
                    eprintf(stderr,
                            "Unknown option character `\\x%x'.\n",
                            optopt);
                }
                print_help(progname, EXIT_FAILURE);
            default:
                eprintf(stderr, "default case in getopt\n");
                exit(EXIT_FAILURE);
        }
    }

    /* select the function */
    if (funcname) {
        if (NULL != strstr(funcname, "profile")) {
            pfunction = parasail_lookup_pfunction(funcname);
            if (NULL == pfunction) {
                eprintf(stderr, "Specified profile function not found.\n");
                exit(EXIT_FAILURE);
            }
            pcreator = parasail_lookup_pcreator(funcname);
            if (NULL == pcreator) {
                eprintf(stderr, "Specified profile creator function not found.\n");
                exit(EXIT_FAILURE);
            }

        }
        else {
            function = parasail_lookup_function(funcname);
            if (NULL == function && NULL != strstr(funcname, "nw_banded")) {
                banded = 1;
            }
            if (NULL == function && 0 == banded) {
                eprintf(stderr, "Specified function not found.\n");
                exit(EXIT_FAILURE);
            }
        }
    }
    else {
        eprintf(stderr, "No alignment function specified.\n");
        exit(EXIT_FAILURE);
    }

    is_stats = (NULL != strstr(funcname, "stats"));
    is_table = (NULL != strstr(funcname, "table"));
    has_query = (NULL != qname);

    if (edge_output && graph_output) {
        eprintf(stderr, "Can only request one of edge or graph output.\n");
        exit(EXIT_FAILURE);
    }
    if ((edge_output || graph_output) && !is_stats) {
        eprintf(stderr, "Edge or graph output requested, but alignment function does not return statistics.\n");
        exit(EXIT_FAILURE);
    }
    if (graph_output && has_query) {
        eprintf(stderr, "Cannot specify a query file and output as a graph.\n");
        exit(EXIT_FAILURE);
    }

    /* select the substitution matrix */
    if (NULL == matrixname && use_dna) {
        matrixname = "ACGT";
        matrix = parasail_matrix_create("ACGT", match, -mismatch);
    }
    else {
        if (NULL == matrixname) {
            matrixname = "blosum62";
        }
        matrix = parasail_matrix_lookup(matrixname);
        if (NULL == matrix) {
            /* try as a filename */
            matrix = parasail_matrix_from_file(matrixname);
        }
        if (NULL == matrix) {
            eprintf(stderr, "Specified substitution matrix not found.\n");
            exit(EXIT_FAILURE);
        }
    }

    if (fname == NULL) {
        eprintf(stderr, "missing input file\n");
        print_help(progname, EXIT_FAILURE);
    }

    /* print the parameters for reference */
    eprintf(stdout,
            "%20s: %s\n"
            "%20s: %d\n"
            "%20s: %s\n"
            "%20s: %d\n"
            "%20s: %d\n"
            "%20s: %s\n"
            "%20s: %d\n"
            "%20s: %d\n"
            "%20s: %d\n"
            "%20s: %s\n"
            "%20s: %s\n"
            "%20s: %s\n",
            "funcname", funcname,
            "cutoff", cutoff,
            "use filter", use_filter ? "yes" : "no",
            "gap_extend", gap_extend,
            "gap_open", gap_open,
            "matrix", matrixname,
            "AOL", AOL,
            "SIM", SIM,
            "OS", OS,
            "file", fname,
            "query", (NULL == qname) ? "<no query>" : qname,
            "output", oname
            );
    if (use_dna) {
        eprintf(stdout,
                "%20s: %d\n"
                "%20s: %d\n",
                "match", match,
                "mismatch", mismatch);
    }

    /* Best to know early whether we can open the output file. */
    if((fop = fopen(oname, "w")) == NULL) {
        eprintf(stderr, "%s: Cannot open output file `%s': ", progname, oname);
        perror("fopen");
        exit(EXIT_FAILURE);
    }
    
    start = parasail_time();
    if (!has_query) {
        parasail_file_t *pf = parasail_open(fname);
        if (fpack) {
            T = (unsigned char*)parasail_read(pf, &n);
        }
        else {
            T = (unsigned char*)parasail_pack(pf, &n);
        }
        parasail_close(pf);
    }
    else {
        parasail_file_t *pf = NULL;
        pf = parasail_open(fname);
        if (fpack) {
            T = (unsigned char*)parasail_read(pf, &t);
        }
        else {
            T = (unsigned char*)parasail_pack(pf, &t);
        }
        parasail_close(pf);
        pf = parasail_open(qname);
        if (qpack) {
            Q = (unsigned char*)parasail_read(pf, &q);
        }
        else {
            Q = (unsigned char*)parasail_pack(pf, &q);
        }
        parasail_close(pf);
        n = t+q;
        /* realloc T and copy Q into it */
        T = (unsigned char*)realloc(T, (n+1)*sizeof(unsigned char));
        if (T == NULL) {
            eprintf(stderr, "%s: Cannot reallocate memory.\n", progname);
            perror("realloc");
            exit(EXIT_FAILURE);
        }
        (void)memcpy(T+t, Q, q);
        free(Q);
    }
    T[n] = '\0';
    finish = parasail_time();
    eprintf(stdout, "%20s: %.4f seconds\n", "read and pack time", finish-start);

    /* Allocate memory for sequence ID array. */
    if (use_filter) {
        SID = (int *)malloc((size_t)n * sizeof(int));
        if(SID == NULL) {
            eprintf(stderr, "%s: Cannot allocate memory.\n", progname);
            perror("malloc");
            exit(EXIT_FAILURE);
        }
    }

    /* determine sentinal */
    if (sentinal == 0) {
        long off = 0;
        while (!isgraph(T[n-off])) {
            ++off;
        }
        sentinal = T[n-off];
    }
    eprintf(stdout, "%20s: %c\n", "sentinal", sentinal);

    /* determine actual end of file (last char) */
    {
        long off = 0;
        while (!isgraph(T[n-off])) {
            ++off;
        }
        n = n - off + 1;
    }
    eprintf(stdout, "%20s: %ld\n", "end of packed buffer", n);

    /* scan T from left to count number of sequences */
    sid = 0;
    for (i=0; i<n; ++i) {
        if (T[i] == sentinal) {
            ++sid;
        }
    }
    if (0 == sid) { /* no sentinal found */
        eprintf(stderr, "no sentinal(%c) found in input\n", sentinal);
        exit(EXIT_FAILURE);
    }
    eprintf(stdout, "%20s: %ld\n", "number of sequences", sid);

    /* scan T from left to build sequence ID and end index */
    /* allocate vectors now that number of sequences is known */
    try {
        BEG.reserve(sid+1);
        END.reserve(sid+1);
        if (use_filter) {
            DB.reserve(sid+1);
        }
    } catch (const bad_alloc&) {
        eprintf(stderr, "Cannot allocate memory for vectors\n");
        exit(EXIT_FAILURE);
    }
    sid = 0;
    BEG.push_back(0);
    if (use_filter) {
        for (i=0; i<n; ++i) {
            SID[i] = sid;
            if (T[i] == sentinal) {
                END.push_back(i);
                BEG.push_back(i+1);
                DB.push_back(i<t);
                if (-1 == sid_crossover && i>=t) {
                    sid_crossover = sid;
                }
                ++sid;
            }
        }
    }
    else {
        for (i=0; i<n; ++i) {
            if (T[i] == sentinal) {
                END.push_back(i);
                BEG.push_back(i+1);
                if (-1 == sid_crossover && i>=t) {
                    sid_crossover = sid;
                }
                ++sid;
            }
        }
    }

    /* if we don't have a query file, clear the DB flags */
    if (!has_query) {
        DB.clear();
        sid_crossover = -1;
    }
    else {
        eprintf(stdout, "%20s: %ld\n", "number of queries", sid - sid_crossover);
        eprintf(stdout, "%20s: %ld\n", "number of db seqs", sid_crossover);
    }

    /* use the enhanced SA filter */
    if (use_filter) {
        /* Allocate memory for enhanced SA. */
        SA = (int *)malloc((size_t)(n+1) * sizeof(int)); /* +1 for LCP */
        LCP = (int *)malloc((size_t)(n+1) * sizeof(int)); /* +1 for lcp tree */
        BWT = (unsigned char *)malloc((size_t)(n+1) * sizeof(unsigned char));
        if((SA == NULL) || (LCP == NULL) || (BWT == NULL))
        {
            eprintf(stderr, "%s: Cannot allocate ESA memory.\n", progname);
            perror("malloc");
            exit(EXIT_FAILURE);
        }

        /* Construct the suffix and LCP arrays.
         * The following sais routine is from Fischer, with bugs fixed. */
        start = parasail_time();
        if(sais(T, SA, LCP, (int)n) != 0) {
            eprintf(stderr, "%s: Cannot allocate memory.\n", progname);
            exit(EXIT_FAILURE);
        }
        finish = parasail_time();
        eprintf(stdout,"%20s: %.4f seconds\n", "induced SA time", finish-start);

        /* construct naive BWT: */
        start = parasail_time();
        for (i = 0; i < n; ++i) {
            BWT[i] = (SA[i] > 0) ? T[SA[i]-1] : sentinal;
        }
        finish = parasail_time();
        eprintf(stdout, "%20s: %.4f seconds\n", "naive BWT time", finish-start);

        /* "fix" the LCP array to clamp LCP's that are too long */
        start = parasail_time();
        for (i = 0; i < n; ++i) {
            int len = END[SID[SA[i]]] - SA[i]; /* don't include sentinal */
            if (LCP[i] > len) LCP[i] = len;
        }
        finish = parasail_time();
        eprintf(stdout, "%20s: %.4f seconds\n", "clamp LCP time", finish-start);

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
            eprintf(stderr, "sentinals not found at beginning or end of SA\n");
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
                    process(count_generated, pairs, last_interval, SA, BWT, SID, DB, sentinal, cutoff);
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
            process(count_generated, pairs, the_stack.top(), SA, BWT, SID, DB, sentinal, cutoff);
        }
        finish = parasail_time();
        if (!has_query) {
            count_possible = ((unsigned long)sid)*((unsigned long)sid-1)/2;
        } else {
            count_possible = (sid-sid_crossover)*sid_crossover;
        }
        eprintf(stdout, "%20s: %.4f seconds\n", "ESA time", finish-start);
        eprintf(stdout, "%20s: %lu\n", "possible pairs", count_possible);
        eprintf(stdout, "%20s: %lu\n", "generated pairs", count_generated);

        /* Deallocate memory. */
        free(SID);
        free(SA);
        free(LCP);
        free(BWT);

        eprintf(stdout, "%20s: %zu\n", "unique pairs", pairs.size());
    }
    else {
        /* don't use enhanced SA filter -- generate all pairs */
        start = parasail_time();
        if (!has_query) {
            /* no query file, so all against all comparison */
            for (int i=0; i<sid; ++i) {
                for (int j=i+1; j<sid; ++j) {
                    vpairs.push_back(make_pair(i,j));
                }
            }
        }
        else {
            /* query given, so only compare query against database */
            for (int i=sid_crossover; i<sid; ++i) {
                for (int j=0; j<sid_crossover; ++j) {
                    vpairs.push_back(make_pair(i,j));
                }
            }
        }
        finish = parasail_time();
        eprintf(stdout, "%20s: %.4f seconds\n", "enumerate time", finish-start);
        eprintf(stdout, "%20s: %zu\n", "unique pairs", vpairs.size());
    }

    if (pairs_only) {
        /* Done with input text. */
        free(T);
        if (vpairs.empty() && !pairs.empty()) {
            for (PairSet::iterator it=pairs.begin(); it!=pairs.end(); ++it) {
                int i = it->first;
                int j = it->second;
                eprintf(fop, "%d,%d\n", i, j);
            }
        }
        else if (!vpairs.empty() && pairs.empty()) {
            for (PairVec::iterator it=vpairs.begin(); it!=vpairs.end(); ++it) {
                int i = it->first;
                int j = it->second;
                eprintf(fop, "%d,%d\n", i, j);
            }
        }
        else {
            eprintf(stderr, "pairs and vpairs were empty\n");
            exit(EXIT_FAILURE);
        }
        fclose(fop);
        return 0;
    }

#ifdef _OPENMP
    if (-1 == num_threads) {
        num_threads = omp_get_max_threads();
    }
    else if (num_threads >= 1) {
        omp_set_num_threads(num_threads);
    }
    else {
        eprintf(stderr, "invalid number of threads chosen (%d)\n", num_threads);
        exit(EXIT_FAILURE);
    }
    eprintf(stdout, "%20s: %d\n", "omp num threads", num_threads);
#endif
#ifdef USE_CILK
    if (-1 == num_threads) {
        /* use defaults */
    }
    else if (num_threads >= 1) {
        char num_threads_str[256];
        sprintf(num_threads_str, "%d", num_threads);
        __cilkrts_set_param("nworkers", num_threads_str);
    }
    else {
        eprintf(stderr, "invalid number of threads chosen (%d)\n", num_threads);
        exit(EXIT_FAILURE);
    }
    eprintf(stdout, "%20s: %d\n", "omp num threads", num_threads);
#endif

    /* OpenMP can't iterate over an STL set. Convert to STL vector. */
    start = parasail_time();
    if (vpairs.empty()) {
        if (pairs.empty()) {
            eprintf(stderr, "pairs and vpairs were empty\n");
            exit(EXIT_FAILURE);
        }
        vpairs.assign(pairs.begin(), pairs.end());
        pairs.clear();
    }
    vector<parasail_result_t*> results(vpairs.size(), NULL);
    finish = parasail_time();
    eprintf(stdout, "%20s: %.4f seconds\n", "openmp prep time", finish-start);

    /* create profiles, if necessary */
    vector<parasail_profile_t*> profiles(sid, (parasail_profile_t*)NULL);
    if (pfunction) {
        start = parasail_time();
        set<int> profile_indices_set;
        for (size_t index=0; index<vpairs.size(); ++index) {
            profile_indices_set.insert(vpairs[index].first);
        }
        vector<int> profile_indices(
                profile_indices_set.begin(),
                profile_indices_set.end());
        //profiles.assign(sid, NULL);
        finish = parasail_time();
        eprintf(stdout, "%20s: %.4f seconds\n", "profile init", finish-start);
        start = parasail_time();
#ifdef USE_CILK
        cilk_for (size_t index=0; index<profile_indices.size(); ++index)
#else
#pragma omp parallel
        {
#pragma omp for schedule(guided)
            for (long long index=0; index<(long long)profile_indices.size(); ++index)
#endif
            {
                int i = profile_indices[index];
                long i_beg = BEG[i];
                long i_end = END[i];
                long i_len = i_end-i_beg;
                profiles[i] = pcreator((const char*)&T[i_beg], i_len, matrix);
            }
#ifdef USE_CILK
#else
        }
#endif
        finish = parasail_time();
        eprintf(stdout, "%20s: %.4f seconds\n", "profile creation", finish-start);
    }

    /* align pairs */
    start = parasail_time();
    if (function) {
#ifdef USE_CILK
            cilk_for (size_t index=0; index<vpairs.size(); ++index)
#else
#pragma omp parallel
            {
#pragma omp for schedule(guided)
            for (long long index=0; index<(long long)vpairs.size(); ++index)
#endif
            {
                int i = vpairs[index].first;
                int j = vpairs[index].second;
                long i_beg = BEG[i];
                long i_end = END[i];
                long i_len = i_end-i_beg;
                long j_beg = BEG[j];
                long j_end = END[j];
                long j_len = j_end-j_beg;
                unsigned long local_work = i_len * j_len;
                parasail_result_t *result = function(
                        (const char*)&T[i_beg], i_len,
                        (const char*)&T[j_beg], j_len,
                        gap_open, gap_extend, matrix);
#ifdef USE_CILK
                work += local_work;
#else
#pragma omp atomic
                work += local_work;
#endif
                results[index] = result;
            }
#ifdef USE_CILK
#else
        }
#endif
    }
    else if (banded) {
#ifdef USE_CILK
            cilk_for (size_t index=0; index<vpairs.size(); ++index)
#else
#pragma omp parallel
            {
#pragma omp for schedule(guided)
            for (long long index=0; index<(long long)vpairs.size(); ++index)
#endif
            {
                int i = vpairs[index].first;
                int j = vpairs[index].second;
                long i_beg = BEG[i];
                long i_end = END[i];
                long i_len = i_end-i_beg;
                long j_beg = BEG[j];
                long j_end = END[j];
                long j_len = j_end-j_beg;
                unsigned long local_work = i_len * j_len;
                parasail_result_t *result = parasail_nw_banded(
                        (const char*)&T[i_beg], i_len,
                        (const char*)&T[j_beg], j_len,
                        gap_open, gap_extend, kbandsize, matrix);
#ifdef USE_CILK
                work += local_work;
#else
#pragma omp atomic
                work += local_work;
#endif
                results[index] = result;
            }
#ifdef USE_CILK
#else
        }
#endif
    }
    else if (pfunction) {
#ifdef USE_CILK
            cilk_for (size_t index=0; index<vpairs.size(); ++index)
#else
#pragma omp parallel
        {
#pragma omp for schedule(guided)
            for (long long index=0; index<(long long)vpairs.size(); ++index)
#endif
            {
                int i = vpairs[index].first;
                int j = vpairs[index].second;
                long j_beg = BEG[j];
                long j_end = END[j];
                long j_len = j_end-j_beg;
                parasail_profile_t *profile = profiles[i];
                if (NULL == profile) {
                    eprintf(stderr, "BAD PROFILE %d\n", i);
                    exit(EXIT_FAILURE);
                }
                unsigned long local_work = profile->s1Len * j_len;
                parasail_result_t *result = pfunction(
                        profile, (const char*)&T[j_beg], j_len,
                        gap_open, gap_extend);
#ifdef USE_CILK
                work += local_work;
#else
#pragma omp atomic
                work += local_work;
#endif
                results[index] = result;
            }
#ifdef USE_CILK
#else
        }
#endif
    }
    else {
        /* shouldn't get here */
        eprintf(stderr, "alignment function was not properly set (shouldn't happen)\n");
        exit(EXIT_FAILURE);
    }
    finish = parasail_time();
#ifdef USE_CILK
    eprintf(stdout, "%20s: %lu cells\n", "work", work.get_value());
#else
    eprintf(stdout, "%20s: %lu cells\n", "work", work);
#endif
    eprintf(stdout, "%20s: %.4f seconds\n", "alignment time", finish-start);
#ifdef USE_CILK
    eprintf(stdout, "%20s: %.4f \n", "gcups", double(work.get_value())/(finish-start)/1000000000);
#else
    eprintf(stdout, "%20s: %.4f \n", "gcups", double(work)/(finish-start)/1000000000);
#endif

    if (pfunction) {
        start = parasail_time();
#ifdef USE_CILK
            cilk_for (size_t index=0; index<profiles.size(); ++index)
#else
#pragma omp parallel
        {
#pragma omp for schedule(guided)
            for (long long index=0; index<(long long)profiles.size(); ++index)
#endif
            {
                if (NULL != profiles[index]) {
                    parasail_profile_free(profiles[index]);
                }
            }
            profiles.clear();
#ifdef USE_CILK
#else
        }
#endif
        finish = parasail_time();
        eprintf(stdout, "%20s: %.4f seconds\n", "profile cleanup", finish-start);
    }

    /* Output results. */
    if (is_stats) {
        if (edge_output) {
            output_edges(fop, has_query, sid_crossover, T, AOL, SIM, OS, matrix, BEG, END, vpairs, results);
        }
        else if (graph_output) {
            output_graph(fop, 0, sid, T, AOL, SIM, OS, matrix, BEG, END, vpairs, results);
        }
        else {
            output_stats(fop, has_query, sid_crossover, BEG, END, vpairs, results);
        }
    }
    else {
        output_basic(fop, has_query, sid_crossover, BEG, END, vpairs, results);
    }
    if (is_table) {
        output_tables(has_query, sid_crossover, T, BEG, END, vpairs, results);
    }

    /* free results */
    for (size_t index=0; index<results.size(); ++index) {
        parasail_result_t *result = results[index];
        parasail_result_free(result);
    }

    /* close output file */
    fclose(fop);

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
        const vector<int> &DB,
        const char &sentinal)
{
    const int &sidi = SID[SA[i]];
    const int &sidj = SID[SA[j]];
    if (BWT[i] != BWT[j] || BWT[i] == sentinal) {
        if (DB.empty()) {
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
        else {
            if (sidi != sidj && DB[sidi] != DB[sidj]) {
                ++count_generated;
                if (sidi > sidj) {
                    pairs.insert(make_pair(sidi,sidj));
                }
                else {
                    pairs.insert(make_pair(sidj,sidi));
                }
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
        const vector<int> &DB,
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
                pair_check(count_generated, pairs, i, j, SA, BWT, SID, DB, sentinal);
            }
        }
    }
    else {
        for (int i=q.lb; i<=q.rb; ++i) {
            for (int j=i+1; j<=q.rb; ++j) {
                pair_check(count_generated, pairs, i, j, SA, BWT, SID, DB, sentinal);
            }
        }
    }
}

#ifdef __MIC__
static const char *get_user_name()
{
    uid_t uid = geteuid();
    struct passwd *pw = getpwuid(uid);
    if (pw) {
        return pw->pw_name;
    }
    return "";
}
#endif

inline static void print_array(
        const char * filename_,
        const int * const restrict array,
        const char * const restrict s1, const int s1Len,
        const char * const restrict s2, const int s2Len)
{
    int i;
    int j;
    FILE *f = NULL;
#ifdef __MIC__
    const char *username = get_user_name();
    char filename[4096] = {0};
    strcat(filename, "/tmp/");
    if (username[0] != '\0') {
        strcat(filename, username);
        strcat(filename, "/");
    }
    strcat(filename, filename_);
#else
    const char *filename = filename_;
#endif
    f = fopen(filename, "w");
    if (NULL == f) {
        printf("fopen(\"%s\") error: %s\n", filename, strerror(errno));
        exit(-1);
    }
    fprintf(f, " ");
    for (j=0; j<s2Len; ++j) {
        fprintf(f, "%4c", s2[j]);
    }
    fprintf(f, "\n");
    for (i=0; i<s1Len; ++i) {
        fprintf(f, "%c", s1[i]);
        for (j=0; j<s2Len; ++j) {
            fprintf(f, "%4d", array[i*s2Len + j]);
        }
        fprintf(f, "\n");
    }
    fclose(f);
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

inline static void output_edges(
        FILE *fop,
        bool has_query,
        long sid_crossover,
        unsigned char *T,
        int AOL,
        int SIM,
        int OS,
        const parasail_matrix_t *matrix,
        const vector<long> &BEG,
        const vector<long> &END,
        const PairVec &vpairs,
        const vector<parasail_result_t*> &results)
{
    unsigned long edge_count = 0;
    for (size_t index=0; index<results.size(); ++index) {
        parasail_result_t *result = results[index];
        int i = vpairs[index].first;
        int j = vpairs[index].second;
        long i_beg = BEG[i];
        long i_end = END[i];
        long i_len = i_end-i_beg;
        long j_beg = BEG[j];
        long j_end = END[j];
        long j_len = j_end-j_beg;

        if (has_query) {
            i = i - sid_crossover;
        }

        int self_score_ = 0;
        int max_len = 0;
        int i_self_score = self_score(
                (const char*)&T[i_beg], i_len, matrix);
        int j_self_score = self_score(
                (const char*)&T[j_beg], j_len, matrix);

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

    fprintf(stdout, "%20s: %lu\n", "edges count", edge_count);
}

inline static void output_graph(
        FILE *fop,
        int which,
        unsigned long sid,
        unsigned char *T,
        int AOL,
        int SIM,
        int OS,
        const parasail_matrix_t *matrix,
        const vector<long> &BEG,
        const vector<long> &END,
        const PairVec &vpairs,
        const vector<parasail_result_t*> &results)
{
    vector<vector<pair<int,float> > > graph(sid);
    unsigned long edge_count = 0;
    for (size_t index=0; index<results.size(); ++index) {
        parasail_result_t *result = results[index];
        int i = vpairs[index].first;
        int j = vpairs[index].second;
        long i_beg = BEG[i];
        long i_end = END[i];
        long i_len = i_end-i_beg;
        long j_beg = BEG[j];
        long j_end = END[j];
        long j_len = j_end-j_beg;

        int self_score_ = 0;
        int max_len = 0;
        int i_self_score = self_score(
                (const char*)&T[i_beg], i_len, matrix);
        int j_self_score = self_score(
                (const char*)&T[j_beg], j_len, matrix);

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
            float value;
            ++edge_count;
            switch (which) {
                case 0:
                    value = 1.0*result->length/max_len;
                    break;
                case 1:
                    value = 1.0*result->matches/result->length;
                    break;
                case 2:
                    value = 1.0*result->score/self_score_;
                    break;
            }
            graph[i].push_back(make_pair(j,value));
            graph[j].push_back(make_pair(i,value));
        }
    }

    fprintf(fop, "%lu %lu 1\n", (unsigned long)graph.size(), edge_count);
    for (size_t i=0; i<graph.size(); ++i) {
        if (graph[i].size() > 0) {
            fprintf(fop, "%d %f", graph[i][0].first+1, graph[i][0].second);
            for (size_t j=1; j<graph[i].size(); ++j) {
                fprintf(fop, "  %d %f", graph[i][j].first+1, graph[i][j].second);
            }
            fprintf(fop, "\n");
        } else {
            fprintf(fop, "\n");
        }
    }

    fprintf(stdout, "%20s: %lu\n", "edges count", edge_count);
}

inline static void output_stats(
        FILE *fop,
        bool has_query,
        long sid_crossover,
        const vector<long> &BEG,
        const vector<long> &END,
        const PairVec &vpairs,
        const vector<parasail_result_t*> &results)
{
    for (size_t index=0; index<results.size(); ++index) {
        parasail_result_t *result = results[index];
        int i = vpairs[index].first;
        int j = vpairs[index].second;
        long i_beg = BEG[i];
        long i_end = END[i];
        long i_len = i_end-i_beg;
        long j_beg = BEG[j];
        long j_end = END[j];
        long j_len = j_end-j_beg;

        if (has_query) {
            i = i - sid_crossover;
        }

        eprintf(fop, "%d,%d,%ld,%ld,%d,%d,%d,%d,%d,%d\n",
                i,
                j,
                i_len,
                j_len,
                result->score,
                result->end_query,
                result->end_ref,
                result->matches,
                result->similar,
                result->length);
    }
}

inline static void output_basic(
        FILE *fop,
        bool has_query,
        long sid_crossover,
        const vector<long> &BEG,
        const vector<long> &END,
        const PairVec &vpairs,
        const vector<parasail_result_t*> &results)
{
    for (size_t index=0; index<results.size(); ++index) {
        parasail_result_t *result = results[index];
        int i = vpairs[index].first;
        int j = vpairs[index].second;
        long i_beg = BEG[i];
        long i_end = END[i];
        long i_len = i_end-i_beg;
        long j_beg = BEG[j];
        long j_end = END[j];
        long j_len = j_end-j_beg;

        if (has_query) {
            i = i - sid_crossover;
        }

        eprintf(fop, "%d,%d,%ld,%ld,%d,%d,%d\n",
                i,
                j,
                i_len,
                j_len,
                result->score,
                result->end_query,
                result->end_ref);
    }
}

inline static void output_tables(
        bool has_query,
        long sid_crossover,
        unsigned char *T,
        const vector<long> &BEG,
        const vector<long> &END,
        const PairVec &vpairs,
        const vector<parasail_result_t*> &results)
{
    for (size_t index=0; index<results.size(); ++index) {
        parasail_result_t *result = results[index];
        int i = vpairs[index].first;
        int j = vpairs[index].second;
        long i_beg = BEG[i];
        long i_end = END[i];
        long i_len = i_end-i_beg;
        long j_beg = BEG[j];
        long j_end = END[j];
        long j_len = j_end-j_beg;

        if (has_query) {
            i = i - sid_crossover;
        }

        char filename[256] = {'\0'};
        sprintf(filename, "parasail_%d_%d.txt", i, j);
        print_array(filename, result->score_table,
                (const char*)&T[i_beg], i_len,
                (const char*)&T[j_beg], j_len);
    }
}

