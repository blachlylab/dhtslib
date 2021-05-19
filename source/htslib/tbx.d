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
module htslib.tbx;
import htslib.hts;
import htslib.bgzf: BGZF;

extern (C):

enum TBX_MAX_SHIFT = 31;

enum TBX_GENERIC = 0;
enum TBX_SAM = 1;
enum TBX_VCF = 2;
enum TBX_UCSC = 0x10000;

struct tbx_conf_t
{
    int preset;
    int sc;
    int bc;
    int ec; // seq col., beg col. and end col.
    int meta_char;
    int line_skip;
}

struct tbx_t
{
    tbx_conf_t conf;
    hts_idx_t* idx;
    void* dict;
}

extern __gshared const tbx_conf_t tbx_conf_gff;
extern __gshared const tbx_conf_t tbx_conf_bed;
extern __gshared const tbx_conf_t tbx_conf_psltbl;
extern __gshared const tbx_conf_t tbx_conf_sam;
extern __gshared const tbx_conf_t tbx_conf_vcf;

alias tbx_itr_destroy = hts_itr_destroy;

pragma(inline, true)
auto tbx_itr_queryi(const tbx_t *tbx, int tid, hts_pos_t beg, hts_pos_t end)
    { return hts_itr_query(tbx.idx, tid, beg, end, &tbx_readrec); }

pragma(inline, true)
auto tbx_itr_querys(const tbx_t *tbx, const char *s)
{
    return hts_itr_querys(tbx.idx, s,
        cast(hts_name2id_f)(&tbx_name2id),
        cast(void*)tbx,
        &hts_itr_query,
        &tbx_readrec);
}

pragma(inline, true)
auto tbx_itr_next(htsFile *htsfp, tbx_t *tbx, hts_itr_t *itr, void *r)
    { return hts_itr_next(hts_get_bgzfp(htsfp), itr, r, tbx); }

pragma(inline, true)
auto tbx_bgzf_itr_next(BGZF *bgzfp, tbx_t *tbx, hts_itr_t *itr, void *r)
    { return hts_itr_next(bgzfp, itr, r, tbx); }

// contig name to integer id
int tbx_name2id (tbx_t* tbx, const(char)* ss);

/* Internal helper function used by tbx_itr_next() defined in hts.c -- do not use directly*/
BGZF* hts_get_bgzfp (htsFile* fp);

int tbx_readrec (
    BGZF* fp,
    void* tbxv,
    void* sv,
    int* tid,
    hts_pos_t* beg,
    hts_pos_t* end);

/// Build an index of the lines in a BGZF-compressed file
/** The index struct returned by a successful call should be freed
    via tbx_destroy() when it is no longer needed.
*/
tbx_t* tbx_index (BGZF* fp, int min_shift, const(tbx_conf_t)* conf);
/*
 * All tbx_index_build* methods return: 0 (success), -1 (general failure) or -2 (compression not BGZF)
 */
int tbx_index_build (const(char)* fn, int min_shift, const(tbx_conf_t)* conf);

int tbx_index_build2 (
    const(char)* fn,
    const(char)* fnidx,
    int min_shift,
    const(tbx_conf_t)* conf);

int tbx_index_build3 (
    const(char)* fn,
    const(char)* fnidx,
    int min_shift,
    int n_threads,
    const(tbx_conf_t)* conf);

/// Load or stream a .tbi or .csi index
/** @param fn     Name of the data file corresponding to the index

    Equivalent to tbx_index_load3(fn, NULL, HTS_IDX_SAVE_REMOTE);
*/
tbx_t* tbx_index_load (const(char)* fn);

/// Load or stream a .tbi or .csi index
/** @param fn     Name of the data file corresponding to the index
    @param fnidx  Name of the indexed file
    @return The index, or NULL if an error occurred

    If @p fnidx is NULL, the index name will be derived from @p fn.

    Equivalent to tbx_index_load3(fn, fnidx, HTS_IDX_SAVE_REMOTE);
*/
tbx_t* tbx_index_load2 (const(char)* fn, const(char)* fnidx);

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
tbx_t* tbx_index_load3 (const(char)* fn, const(char)* fnidx, int flags);

/// return C-style array of sequence names (NB: free the array but not the values)
const(char*)* tbx_seqnames (tbx_t* tbx, int* n); // free the array but not the values

/// destroy/dealloc tabix data
void tbx_destroy (tbx_t* tbx);

