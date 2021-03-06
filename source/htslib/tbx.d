/** htslib-1.9 tbx.h as D module
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
    Copyright (C) 2009, 2012-2015, 2019 Genome Research Ltd.
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

enum int TBX_MAX_SHIFT = 31;/// ???

enum int TBX_GENERIC = 0;   /// generic flat file
enum int TBX_SAM     = 1;   /// SAM
enum int TBX_VCF     = 2;   /// VCF
enum int TBX_UCSC    = 0x10000; /// ?UCSC flat file?

/// tabix config
struct tbx_conf_t {
    int32_t preset;     /// ?
    /// seq col., beg col. and end col.
    int32_t sc, bc, ec; // seq col., beg col. and end col.
    /// ?
    int32_t meta_char, line_skip;
}

/// tabix data
struct tbx_t {
    tbx_conf_t conf;    /// tabix config
    hts_idx_t *idx;     /// index data
    void *dict;         /// ?dictionary
}

//extern const tbx_conf_t tbx_conf_gff, tbx_conf_bed, tbx_conf_psltbl, tbx_conf_sam, tbx_conf_vcf;
/// prebaked TABIX config data for GFF3, BED, PSL table, SAM, VCF
extern (C) extern __gshared const tbx_conf_t tbx_conf_gff, tbx_conf_bed, tbx_conf_psltbl, tbx_conf_sam, tbx_conf_vcf;

    alias tbx_itr_destroy = hts_itr_destroy;

    /* hts_itr_t *hts_itr_query(const hts_idx_t *idx, int tid, int beg, int end, hts_readrec_func *readrec); */
    //#define tbx_itr_queryi(tbx, tid, beg, end) hts_itr_query((tbx)->idx, (tid), (beg), (end), tbx_readrec)
    /// tabix query by integer based tid(contig)/start/end
    pragma(inline, true)
    auto tbx_itr_queryi(const tbx_t *tbx, int tid, int beg, int end)
        { return hts_itr_query(tbx.idx, tid, beg, end, &tbx_readrec); }

    /* hts_itr_t *hts_itr_querys(const hts_idx_t *idx, const char *reg, hts_name2id_f getid, void *hdr, hts_itr_query_func *itr_query, hts_readrec_func *readrec); */
    //#define tbx_itr_querys(tbx, s) hts_itr_querys((tbx)->idx, (s), (hts_name2id_f)(tbx_name2id), (tbx), hts_itr_query, tbx_readrec)
    /// tabix query by string "chr:start-end"
    pragma(inline, true)
    auto tbx_itr_querys(const tbx_t *tbx, const char *s)
    {
        return hts_itr_querys(tbx.idx, s,
            cast(hts_name2id_f)(&tbx_name2id),
            cast(void*)tbx,
            &hts_itr_query,
            &tbx_readrec);
    }

    /* int hts_itr_next(BGZF *fp, hts_itr_t *iter, void *r, void *data) HTS_RESULT_USED; */
    //#define tbx_itr_next(htsfp, tbx, itr, r) hts_itr_next(hts_get_bgzfp(htsfp), (itr), (r), (tbx))
    /// advance tabix iterator
    pragma(inline, true)
    auto tbx_itr_next(htsFile *htsfp, tbx_t *tbx, hts_itr_t *itr, void *r)
        { return hts_itr_next(hts_get_bgzfp(htsfp), itr, r, tbx); }
    
    /* int hts_itr_next(BGZF *fp, hts_itr_t *iter, void *r, void *data) HTS_RESULT_USED; */
    //#define tbx_bgzf_itr_next(bgzfp, tbx, itr, r) hts_itr_next((bgzfp), (itr), (r), (tbx))
    /// advance tabix iterator
    pragma(inline, true)
    auto tbx_bgzf_itr_next(BGZF *bgzfp, tbx_t *tbx, hts_itr_t *itr, void *r)
        { return hts_itr_next(bgzfp, itr, r, tbx); }

    /// contig name to integer id
    int tbx_name2id(tbx_t *tbx, const char *ss);

    /** Internal helper function used by tbx_itr_next()  defined in hts.c -- do not use directly */
    BGZF *hts_get_bgzfp(htsFile *fp);
    /** Called by tabix iterator to read the next record */
    int tbx_readrec(BGZF *fp, void *tbxv, void *sv, int *tid, hts_pos_t *beg, hts_pos_t *end);

/// Build an index of the lines in a BGZF-compressed file
/** The index struct returned by a successful call should be freed
    via tbx_destroy() when it is no longer needed.
*/
    tbx_t *tbx_index(BGZF *fp, int min_shift, const(tbx_conf_t) *conf);
/*
 * All tbx_index_build* methods return: 0 (success), -1 (general failure) or -2 (compression not BGZF)
 */
    int tbx_index_build(const(char) *fn, int min_shift, const(tbx_conf_t) *conf);
    /// ditto
    int tbx_index_build2(const(char) *fn, const(char) *fnidx, int min_shift, const(tbx_conf_t) *conf);
    /// ditto
    int tbx_index_build3(const(char) *fn, const(char) *fnidx, int min_shift, int n_threads, const(tbx_conf_t) *conf);
    
/// Load or stream a .tbi or .csi index
/** @param fn     Name of the data file corresponding to the index

    Equivalent to tbx_index_load3(fn, NULL, HTS_IDX_SAVE_REMOTE);
*/
    tbx_t *tbx_index_load(const(char) *fn);

/// Load or stream a .tbi or .csi index
/** @param fn     Name of the data file corresponding to the index
    @param fnidx  Name of the indexed file
    @return The index, or NULL if an error occurred

    If @p fnidx is NULL, the index name will be derived from @p fn.

    Equivalent to tbx_index_load3(fn, fnidx, HTS_IDX_SAVE_REMOTE);
*/
    tbx_t *tbx_index_load2(const(char) *fn, const(char) *fnidx);

/// Load or stream a .tbi or .csi index
/** @param fn     Name of the data file corresponding to the index
    @param fnidx  Name of the indexed file
    @param flags  Flags to alter behaviour (see description)
    @return The index, or NULL if an error occurred

    If @p fnidx is NULL, the index name will be derived from @p fn.

    The @p flags parameter can be set to a combination of the following
    values:

        HTS_IDX_SAVE_REMOTE   Save a local copy of any remote indexes
        HTS_IDX_SILENT_FAIL   Fail silently if the index is not present

    The index struct returned by a successful call should be freed
    via tbx_destroy() when it is no longer needed.
*/
    tbx_t *tbx_index_load3(const(char) *fn, const(char) *fnidx, HTS_IDX_FLAG flags);

    /// return C-style array of sequence names (NB: free the array but not the values)
    const(char **) tbx_seqnames(tbx_t *tbx, int *n);  // free the array but not the values

    /// destroy/dealloc tabix data
    void tbx_destroy(tbx_t *tbx);
