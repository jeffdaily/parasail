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

inline static int self_score(
        const char * const restrict seq,
        int len,
        const parasail_matrix_t *matrix);

inline static void read_and_pack_file(
        const char *fname, 
        unsigned char * &T,
        long &n,
        const char *progname);

static void print_help(const char *progname, int status) {
    fprintf(stderr, "\nusage: %s "
            "[-a funcname] "
            "[-c cutoff] "
            "[-x] "
            "[-e gap_extend] "
            "[-o gap_open] "
            "[-m matrix] "
            "-f file "
            "[-q query_file] "
            "[-g output_file] "
            "\n\n",
            progname);
    fprintf(stderr, "Defaults:\n"
            "   funcname: sw_stats_striped_16\n"
            "     cutoff: 7, must be >= 1, exact match length cutoff\n"
            "         -x: if present, don't use suffix array filter\n"
            " gap_extend: 1, must be >= 0\n"
            "   gap_open: 10, must be >= 0\n"
            "     matrix: blosum62\n"
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
#ifdef _OPENMP
    int num_threads = -1;
#endif
    int *SA = NULL;
    int *LCP = NULL;
    unsigned char *BWT = NULL;
    int *SID = NULL;
    vector<int> BEG;
    vector<int> END;
    vector<int> DB;
    long n = 0;
    long t = 0;
    long q = 0;
    double start = 0;
    double finish = 0;
    int i = 0;
    int sid = 0;
    int sid_crossover = -1;
    char sentinal = 0;
    int cutoff = 7;
    bool use_filter = true;
    PairSet pairs;
    unsigned long count_possible = 0;
    unsigned long count_generated = 0;
    unsigned long work = 0;
    int c = 0;
    const char *funcname = "sw_stats_striped_16";
    parasail_function_t *function = NULL;
    parasail_pfunction_t *pfunction = NULL;
    parasail_pcreator_t *pcreator = NULL;
    const char *matrixname = "blosum62";
    const parasail_matrix_t *matrix = NULL;
    int gap_open = 10;
    int gap_extend = 1;
    const char *progname = "parasail_aligner";

    /* Check arguments. */
    while ((c = getopt(argc, argv, "a:c:e:f:g:hm:o:q:t:x")) != -1) {
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
            case 'q':
                qname = optarg;
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
            case 't':
                num_threads = atoi(optarg);
                break;
            case 'x':
                use_filter = false;
                break;
            case '?':
                if (optopt == 'a'
                        || optopt == 'c'
                        || optopt == 'e'
                        || optopt == 'f'
                        || optopt == 'g'
                        || optopt == 'm'
                        || optopt == 'o'
                        || optopt == 'q'
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
                print_help(progname, EXIT_FAILURE);
            default:
                fprintf(stderr, "default case in getopt\n");
                exit(EXIT_FAILURE);
        }
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

    /* print the parameters for reference */
    fprintf(stdout,
            "%20s: %s\n"
            "%20s: %d\n"
            "%20s: %s\n"
            "%20s: %d\n"
            "%20s: %d\n"
            "%20s: %s\n"
            "%20s: %s\n"
            "%20s: %s\n"
            "%20s: %s\n",
            "funcname", funcname,
            "cutoff", cutoff,
            "use filter", use_filter ? "yes" : "no",
            "gap_extend", gap_extend,
            "gap_open", gap_open,
            "matrix", matrixname,
            "file", fname,
            "query", (NULL == qname) ? "<no query>" : qname,
            "output", oname
            );

    /* Best to know early whether we can open the output file. */
    if((fop = fopen(oname, "w")) == NULL) {
        fprintf(stderr, "%s: Cannot open output file `%s': ", progname, oname);
        perror("fopen");
        exit(EXIT_FAILURE);
    }
    
    start = parasail_time();
    if (qname == NULL) {
        read_and_pack_file(fname, T, n, progname);
    }
    else {
        read_and_pack_file(fname, T, t, progname);
        read_and_pack_file(qname, Q, q, progname);
        n = t+q;
        /* realloc T and copy Q into it */
        T = (unsigned char*)realloc(T, (n+1)*sizeof(unsigned char));
        if (T == NULL) {
            fprintf(stderr, "%s: Cannot reallocate memory.\n", progname);
            perror("realloc");
            exit(EXIT_FAILURE);
        }
        (void)memcpy(T+t, Q, q);
        free(Q);
    }
    T[n] = '\0';
    finish = parasail_time();
    fprintf(stdout, "%20s: %.4f seconds\n", "read and pack time", finish-start);

    /* Allocate memory for sequence ID array. */
    SID = (int *)malloc((size_t)n * sizeof(int));
    if(SID == NULL) {
        fprintf(stderr, "%s: Cannot allocate memory.\n", progname);
        perror("malloc");
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
    sid = 0;
    BEG.push_back(0);
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
    if (0 == sid) { /* no sentinal found */
        fprintf(stderr, "no sentinal(%c) found in input\n", sentinal);
        exit(EXIT_FAILURE);
    }
    fprintf(stdout, "%20s: %d\n", "number of sequences", sid);
    /* if we don't have a query file, clear the DB flags */
    if (qname == NULL) {
        DB.clear();
        sid_crossover = -1;
    }
    else {
        fprintf(stdout, "%20s: %d\n", "number of queries", sid - sid_crossover);
        fprintf(stdout, "%20s: %d\n", "number of db seqs", sid_crossover);
    }

    /* use the enhanced SA filter */
    if (use_filter) {
        /* Allocate memory for enhanced SA. */
        SA = (int *)malloc((size_t)(n+1) * sizeof(int)); /* +1 for LCP */
        LCP = (int *)malloc((size_t)(n+1) * sizeof(int)); /* +1 for lcp tree */
        BWT = (unsigned char *)malloc((size_t)(n+1) * sizeof(unsigned char));
        if((SA == NULL) || (LCP == NULL) || (BWT == NULL))
        {
            fprintf(stderr, "%s: Cannot allocate ESA memory.\n", progname);
            perror("malloc");
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
        fprintf(stdout,"%20s: %.4f seconds\n", "induced SA time", finish-start);

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
        count_possible = ((unsigned long)sid)*((unsigned long)sid-1)/2;
        fprintf(stdout, "%20s: %.4f seconds\n", "ESA time", finish-start);
        fprintf(stdout, "%20s: %lu\n", "possible pairs", count_possible);
        fprintf(stdout, "%20s: %lu\n", "generated pairs", count_generated);

        /* Deallocate memory. */
        free(SA);
        free(LCP);
        free(BWT);
    }
    else {
        /* don't use enhanced SA filter -- generate all pairs */
        start = parasail_time();
        if (qname == NULL) {
            /* no query file, so all against all comparison */
            for (int i=0; i<sid; ++i) {
                for (int j=i+1; j<sid; ++j) {
                    pairs.insert(make_pair(i,j));
                }
            }
        }
        else {
            /* query given, so only compare query against database */
            for (int i=sid_crossover; i<sid; ++i) {
                for (int j=0; j<sid_crossover; ++j) {
                    pairs.insert(make_pair(i,j));
                }
            }
        }
        finish = parasail_time();
        fprintf(stdout, "%20s: %.4f seconds\n", "enumerate time", finish-start);
    }
    fprintf(stdout, "%20s: %zu\n", "unique pairs", pairs.size());

    /* Deallocate memory. */
    free(SID);

#ifdef _OPENMP
    if (-1 == num_threads) {
        num_threads = omp_get_max_threads();
    }
    else if (num_threads >= 1) {
        omp_set_num_threads(num_threads);
    }
    else {
        fprintf(stderr, "invalid number of threads chosen (%d)\n", num_threads);
        exit(EXIT_FAILURE);
    }
    fprintf(stdout, "%20s: %d\n", "omp num threads", num_threads);
#endif

    /* OpenMP can't iterate over an STL set. Convert to STL vector. */
    start = parasail_time();
    vector<Pair> vpairs(pairs.begin(), pairs.end());
    vector<parasail_result_t*> results(vpairs.size(), NULL);
    finish = parasail_time();
    fprintf(stdout, "%20s: %.4f seconds\n", "openmp prep time", finish-start);

    /* create profiles, if necessary */
    vector<parasail_profile_t*> profiles;
    if (pfunction) {
        start = parasail_time();
        set<int> profile_indices_set;
        for (size_t index=0; index<vpairs.size(); ++index) {
            profile_indices_set.insert(vpairs[index].first);
        }
        vector<int> profile_indices(
                profile_indices_set.begin(),
                profile_indices_set.end());
        profiles.assign(sid, NULL);
        finish = parasail_time();
        fprintf(stdout, "%20s: %.4f seconds\n", "profile init", finish-start);
        start = parasail_time();
#pragma omp parallel
        {
#pragma omp for schedule(dynamic)
            for (size_t index=0; index<profile_indices.size(); ++index) {
                int i = profile_indices[index];
                int i_beg = BEG[i];
                int i_end = END[i];
                int i_len = i_end-i_beg;
                profiles[i] = pcreator((const char*)&T[i_beg], i_len, matrix);
            }
        }
        finish = parasail_time();
        fprintf(stdout, "%20s: %.4f seconds\n", "profile creation", finish-start);
    }

    /* align pairs */
    start = parasail_time();
    if (function) {
#pragma omp parallel
        {
#pragma omp for schedule(dynamic)
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
#pragma omp atomic
                work += local_work;
                results[index] = result;
            }
        }
    }
    else if (pfunction) {
#pragma omp parallel
        {
#pragma omp for schedule(dynamic)
            for (size_t index=0; index<vpairs.size(); ++index) {
                int i = vpairs[index].first;
                int j = vpairs[index].second;
                int j_beg = BEG[j];
                int j_end = END[j];
                int j_len = j_end-j_beg;
                parasail_profile_t *profile = profiles[i];
                if (NULL == profile) {
                    fprintf(stderr, "BAD PROFILE %d\n", i);
                    exit(EXIT_FAILURE);
                }
                unsigned long local_work = profile->s1Len * j_len;
                parasail_result_t *result = pfunction(
                        profile, (const char*)&T[j_beg], j_len,
                        gap_open, gap_extend);
#pragma omp atomic
                work += local_work;
                results[index] = result;
            }
        }
    }
    else {
        /* shouldn't get here */
        fprintf(stderr, "alignment function was not properly set (shouldn't happen)\n");
        exit(EXIT_FAILURE);
    }
    finish = parasail_time();
    fprintf(stdout, "%20s: %lu cells\n", "work", work);
    fprintf(stdout, "%20s: %.4f seconds\n", "alignment time", finish-start);
    fprintf(stdout, "%20s: %.4f \n", "gcups", double(work)/(finish-start)/1000000000);

    if (pfunction) {
        start = parasail_time();
#pragma omp parallel
        {
#pragma omp for schedule(dynamic)
            for (size_t index=0; index<profiles.size(); ++index) {
                if (NULL != profiles[index]) {
                    parasail_profile_free(profiles[index]);
                }
            }
            profiles.clear();
        }
        finish = parasail_time();
        fprintf(stdout, "%20s: %.4f seconds\n", "profile cleanup", finish-start);
    }

    /* Output results. */
    bool is_stats = (NULL != strstr(funcname, "stats"));
    for (size_t index=0; index<results.size(); ++index) {
        parasail_result_t *result = results[index];
        int i = vpairs[index].first;
        int j = vpairs[index].second;

        if (is_stats) {
            fprintf(fop, "%d,%d,%d,%d,%d,%d\n",
                    (NULL == qname) ? i : i - sid_crossover,
                    j,
                    result->score,
                    result->matches,
                    result->similar,
                    result->length);
        }
        else {
            fprintf(fop, "%d,%d,%d\n",
                    (NULL == qname) ? i : i - sid_crossover,
                    j,
                    result->score);
        }

        parasail_result_free(result);
    }
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

inline static void read_and_pack_file(
        const char *fname, 
        unsigned char * &T,
        long &n,
        const char *progname)
{
    FILE *fip = NULL;

    /* Open a file for reading. */
    if((fip = fopen(fname, "r")) == NULL) {
        fprintf(stderr, "%s: Cannot open input file `%s': ", progname, fname);
        perror("fopen");
        exit(EXIT_FAILURE);
    }

    /* Get the database file size. */
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
    if (T == NULL) {
        fprintf(stderr, "%s: Cannot allocate memory for file.\n", progname);
        perror("malloc");
        exit(EXIT_FAILURE);
    }
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
        long newlines = 0;
        for (i=0; i<n; ++i) { 
            if (T[i] == '>') {
                if (first) {
                    first = 0;
                }
                else {
                    T[w++] = '$';
                }
                /* skip rest of this line */
                while (T[i] != '\n' && T[i] != '\r') {
                    ++i;
                }
                newlines++;
                /* for the case of "\r\n" */
                if (T[i] == '\n') {
                    ++i;
                }
            }
            else if (isalpha(T[i])) {
                T[w++] = T[i];
            }
            else if (T[i] == '\n' || T[i] == '\r') {
                /* ignore newline */
                newlines++;
                /* for the case of "\r\n" */
                if (T[i] == '\r') {
                    ++i;
                    if (i >= n || T[i] != '\n') {
                        fprintf(stderr, "error: \\r without \\n "
                                "line %ld in input\n", newlines);
                        exit(EXIT_FAILURE);
                    }
                }
            }
            else if (isprint(T[i])) {
                fprintf(stderr, "error: non-alpha character "
                        "at pos %ld line %ld in input ('%c')\n",
                        i, newlines, T[i]);
                exit(EXIT_FAILURE);
            }
            else {
                fprintf(stderr, "error: non-printing character "
                        "at pos %ld line %ld in input ('%d')\n",
                        i, newlines, T[i]);
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
}

