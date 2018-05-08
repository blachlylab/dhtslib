/** htslib-1.7 bgzf.h as D module
 *
 *  Changes include:
 *      Removed if(n)defs
 *      Change numeric #defines to enum int
 *      Changed ^typedef struct {...} <name>$ to ^struct <name> {...}$
 *      extern const to __gshared
 *      made #define function macros into inline functions (tbx_itr* -> hts_itr*) 
 *      In D, const on either LHS or RHS of function declaration applies to the function, not return value, unless parents included:
 *      changed ^const <type> <fnname> to ^const(<type>) <fnname>
 */
module htslib.tbx;

import std.stdint : int32_t;

import htslib.hts;
import htslib.bgzf;

extern (C):

/// @file htslib/tbx.h
/// Tabix API functions.
/*
    Copyright (C) 2009, 2012-2015 Genome Research Ltd.
    Copyright (C) 2010, 2012 Broad Institute.
    Author: Heng Li <lh3@sanger.ac.uk>
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

enum int TBX_MAX_SHIFT = 31;

enum int TBX_GENERIC = 0;
enum int TBX_SAM     = 1;
enum int TBX_VCF     = 2;
enum int TBX_UCSC    = 0x10000;

struct tbx_conf_t {
    int32_t preset;
    int32_t sc, bc, ec; // seq col., beg col. and end col.
    int32_t meta_char, line_skip;
};

struct tbx_t {
    tbx_conf_t conf;
    hts_idx_t *idx;
    void *dict;
};

//extern const tbx_conf_t tbx_conf_gff, tbx_conf_bed, tbx_conf_psltbl, tbx_conf_sam, tbx_conf_vcf;
extern (C) extern __gshared immutable tbx_conf_t tbx_conf_gff, tbx_conf_bed, tbx_conf_psltbl, tbx_conf_sam, tbx_conf_vcf;

    //#define tbx_itr_destroy(iter) hts_itr_destroy(iter)
    alias tbx_itr_destroy = hts_itr_destroy;

    /* hts_itr_t *hts_itr_query(const hts_idx_t *idx, int tid, int beg, int end, hts_readrec_func *readrec); */
    //#define tbx_itr_queryi(tbx, tid, beg, end) hts_itr_query((tbx)->idx, (tid), (beg), (end), tbx_readrec)
    auto tbx_itr_queryi(const tbx_t *tbx, int tid, int beg, int end) { return hts_itr_query(tbx.idx, tid, beg, end, &tbx_readrec); }

    /* hts_itr_t *hts_itr_querys(const hts_idx_t *idx, const char *reg, hts_name2id_f getid, void *hdr, hts_itr_query_func *itr_query, hts_readrec_func *readrec); */
    //#define tbx_itr_querys(tbx, s) hts_itr_querys((tbx)->idx, (s), (hts_name2id_f)(tbx_name2id), (tbx), hts_itr_query, tbx_readrec)
    auto tbx_itr_querys(const tbx_t *tbx, const char *s) { return hts_itr_querys(tbx.idx, s, cast(hts_name2id_f)(&tbx_name2id), cast(void*)tbx, &hts_itr_query, &tbx_readrec); }

    /* int hts_itr_next(BGZF *fp, hts_itr_t *iter, void *r, void *data) HTS_RESULT_USED; */
    //#define tbx_itr_next(htsfp, tbx, itr, r) hts_itr_next(hts_get_bgzfp(htsfp), (itr), (r), (tbx))
    auto tbx_itr_next(htsFile *htsfp, tbx_t *tbx, hts_itr_t *itr, void *r) { return hts_itr_next(hts_get_bgzfp(htsfp), itr, r, tbx); }
    
    //#define tbx_bgzf_itr_next(bgzfp, tbx, itr, r) hts_itr_next((bgzfp), (itr), (r), (tbx))
    auto tbx_bgzf_itr_next(BGZF *bgzfp, tbx_t *tbx, hts_itr_t *itr, void *r) { return hts_itr_next(bgzfp, itr, r, tbx); }

    int tbx_name2id(tbx_t *tbx, const char *ss);

    /* Internal helper function used by tbx_itr_next() */
    BGZF *hts_get_bgzfp(htsFile *fp);
    int tbx_readrec(BGZF *fp, void *tbxv, void *sv, int *tid, int *beg, int *end);

    tbx_t *tbx_index(BGZF *fp, int min_shift, const tbx_conf_t *conf);
    int tbx_index_build(const char *fn, int min_shift, const tbx_conf_t *conf);
    int tbx_index_build2(const char *fn, const char *fnidx, int min_shift, const tbx_conf_t *conf);
    int tbx_index_build3(const char *fn, const char *fnidx, int min_shift, int n_threads, const tbx_conf_t *conf);
    tbx_t *tbx_index_load(const char *fn);
    tbx_t *tbx_index_load2(const char *fn, const char *fnidx);
    const(char **) tbx_seqnames(tbx_t *tbx, int *n);  // free the array but not the values
    void tbx_destroy(tbx_t *tbx);