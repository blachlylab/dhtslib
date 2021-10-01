/// @file htslib/vcf.h
/// High-level VCF/BCF variant calling file operations.
/*
    Copyright (C) 2012, 2013 Broad Institute.
    Copyright (C) 2012-2020 Genome Research Ltd.

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

/*
    todo:
        - make the function names consistent
        - provide calls to abstract away structs as much as possible
 */

import core.stdc.config;

extern (C):

/* Included only for backwards compatibility with e.g. bcftools 1.10 */

/*****************
 * Header struct *
 *****************/

enum BCF_HL_FLT = 0; // header line
enum BCF_HL_INFO = 1;
enum BCF_HL_FMT = 2;
enum BCF_HL_CTG = 3;
enum BCF_HL_STR = 4; // structured header line TAG=<A=..,B=..>
enum BCF_HL_GEN = 5; // generic header line

enum BCF_HT_FLAG = 0; // header type
enum BCF_HT_INT = 1;
enum BCF_HT_REAL = 2;
enum BCF_HT_STR = 3;
enum BCF_HT_LONG = BCF_HT_INT | 0x100; // BCF_HT_INT, but for int64_t values; VCF only!

enum BCF_VL_FIXED = 0; // variable length
enum BCF_VL_VAR = 1;
enum BCF_VL_A = 2;
enum BCF_VL_G = 3;
enum BCF_VL_R = 4;

/* === Dictionary ===

   The header keeps three dictionaries. The first keeps IDs in the
   "FILTER/INFO/FORMAT" lines, the second keeps the sequence names and lengths
   in the "contig" lines and the last keeps the sample names. bcf_hdr_t::dict[]
   is the actual hash table, which is opaque to the end users. In the hash
   table, the key is the ID or sample name as a C string and the value is a
   bcf_idinfo_t struct. bcf_hdr_t::id[] points to key-value pairs in the hash
   table in the order that they appear in the VCF header. bcf_hdr_t::n[] is the
   size of the hash table or, equivalently, the length of the id[] arrays.
*/

enum BCF_DT_ID = 0; // dictionary type
enum BCF_DT_CTG = 1;
enum BCF_DT_SAMPLE = 2;

// Complete textual representation of a header line
struct bcf_hrec_t
{
    int type; // One of the BCF_HL_* type
    char* key; // The part before '=', i.e. FILTER/INFO/FORMAT/contig/fileformat etc.
    char* value; // Set only for generic lines, NULL for FILTER/INFO, etc.
    int nkeys; // Number of structured fields
    char** keys;
    char** vals; // The key=value pairs
}

struct bcf_idinfo_t
{
    ulong[3] info; // stores Number:20, var:4, Type:4, ColType:4 in info[0..2]
    // for BCF_HL_FLT,INFO,FMT and contig length in info[0] for BCF_HL_CTG
    bcf_hrec_t*[3] hrec;
    int id;
}

struct bcf_idpair_t
{
    const(char)* key;
    const(bcf_idinfo_t)* val;
}

// Note that bcf_hdr_t structs must always be created via bcf_hdr_init()
struct bcf_hdr_t
{
    int[3] n; // n:the size of the dictionary block in use, (allocated size, m, is below to preserve ABI)
    bcf_idpair_t*[3] id;
    void*[3] dict; // ID dictionary, contig dict and sample dict
    char** samples;
    bcf_hrec_t** hrec;
    int nhrec;
    int dirty;
    int ntransl;
    int*[2] transl; // for bcf_translate()
    int nsamples_ori; // for bcf_hdr_set_samples()
    ubyte* keep_samples;
    kstring_t mem;
    int[3] m; // m: allocated size of the dictionary block in use (see n above)
}

extern __gshared ubyte[] bcf_type_shift;

/**************
 * VCF record *
 **************/

enum BCF_BT_NULL = 0;
enum BCF_BT_INT8 = 1;
enum BCF_BT_INT16 = 2;
enum BCF_BT_INT32 = 3;
enum BCF_BT_INT64 = 4; // Unofficial, for internal use only.
enum BCF_BT_FLOAT = 5;
enum BCF_BT_CHAR = 7;

enum VCF_REF = 0;
enum VCF_SNP = 1;
enum VCF_MNP = 2;
enum VCF_INDEL = 4;
enum VCF_OTHER = 8;
enum VCF_BND = 16; // breakend
enum VCF_OVERLAP = 32; // overlapping deletion, ALT=*

struct bcf_variant_t
{
    int type;
    int n; // variant type and the number of bases affected, negative for deletions
}

struct bcf_fmt_t
{
    import std.bitmanip : bitfields;

    int id; // id: numeric tag id, the corresponding string is bcf_hdr_t::id[BCF_DT_ID][$id].key
    int n;
    int size;
    int type; // n: number of values per-sample; size: number of bytes per-sample; type: one of BCF_BT_* types
    ubyte* p; // same as vptr and vptr_* in bcf_info_t below
    uint p_len;

    mixin(bitfields!(
        uint, "p_off", 31,
        uint, "p_free", 1));
}

struct bcf_info_t
{
    import std.bitmanip : bitfields;

    int key; // key: numeric tag id, the corresponding string is bcf_hdr_t::id[BCF_DT_ID][$key].key
    int type; // type: one of BCF_BT_* types

    // integer value
    // float value
    union _Anonymous_0
    {
        long i;
        float f;
    }

    _Anonymous_0 v1; // only set if $len==1; for easier access
    ubyte* vptr; // pointer to data array in bcf1_t->shared.s, excluding the size+type and tag id bytes
    uint vptr_len;

    mixin(bitfields!(
        uint, "vptr_off", 31,
        uint, "vptr_free", 1)); // length of the vptr block or, when set, of the vptr_mod block, excluding offset
    // vptr offset, i.e., the size of the INFO key plus size+type bytes
    // indicates that vptr-vptr_off must be freed; set only when modified and the new
    //    data block is bigger than the original
    int len; // vector length, 1 for scalars
}

enum BCF1_DIRTY_ID = 1;
enum BCF1_DIRTY_ALS = 2;
enum BCF1_DIRTY_FLT = 4;
enum BCF1_DIRTY_INF = 8;

struct bcf_dec_t
{
    int m_fmt;
    int m_info;
    int m_id;
    int m_als;
    int m_allele;
    int m_flt; // allocated size (high-water mark); do not change
    int n_flt; // Number of FILTER fields
    int* flt; // FILTER keys in the dictionary
    char* id;
    char* als; // ID and REF+ALT block (\0-separated)
    char** allele; // allele[0] is the REF (allele[] pointers to the als block); all null terminated
    bcf_info_t* info; // INFO
    bcf_fmt_t* fmt; // FORMAT and individual sample
    bcf_variant_t* var; // $var and $var_type set only when set_variant_types called
    int n_var;
    int var_type;
    int shared_dirty; // if set, shared.s must be recreated on BCF output
    int indiv_dirty; // if set, indiv.s must be recreated on BCF output
}

enum BCF_ERR_CTG_UNDEF = 1;
enum BCF_ERR_TAG_UNDEF = 2;
enum BCF_ERR_NCOLS = 4;
enum BCF_ERR_LIMITS = 8;
enum BCF_ERR_CHAR = 16;
enum BCF_ERR_CTG_INVALID = 32;
enum BCF_ERR_TAG_INVALID = 64;

/*
    The bcf1_t structure corresponds to one VCF/BCF line. Reading from VCF file
    is slower because the string is first to be parsed, packed into BCF line
    (done in vcf_parse), then unpacked into internal bcf1_t structure. If it
    is known in advance that some of the fields will not be required (notably
    the sample columns), parsing of these can be skipped by setting max_unpack
    appropriately.
    Similarly, it is fast to output a BCF line because the columns (kept in
    shared.s, indiv.s, etc.) are written directly by bcf_write, whereas a VCF
    line must be formatted in vcf_format.
 */
struct bcf1_t
{
    import std.bitmanip : bitfields;

    hts_pos_t pos; // POS
    hts_pos_t rlen; // length of REF
    int rid; // CHROM
    float qual;

    mixin(bitfields!(
        uint, "n_info", 16,
        uint, "n_allele", 16,
        uint, "n_fmt", 8,
        uint, "n_sample", 24)); // QUAL

    kstring_t shared_;
    kstring_t indiv;
    bcf_dec_t d; // lazy evaluation: $d is not generated by bcf_read(), but by explicitly calling bcf_unpack()
    int max_unpack; // Set to BCF_UN_STR, BCF_UN_FLT, or BCF_UN_INFO to boost performance of vcf_parse when some of the fields won't be needed
    int unpacked; // remember what has been unpacked to allow calling bcf_unpack() repeatedly without redoing the work
    int[3] unpack_size; // the original block size of ID, REF+ALT and FILTER
    int errcode; // one of BCF_ERR_* codes
}

/*******
 * API *
 *******/

/***********************************************************************
 *  BCF and VCF I/O
 *
 *  A note about naming conventions: htslib internally represents VCF
 *  records as bcf1_t data structures, therefore most functions are
 *  prefixed with bcf_. There are a few exceptions where the functions must
 *  be aware of both BCF and VCF worlds, such as bcf_parse vs vcf_parse. In
 *  these cases, functions prefixed with bcf_ are more general and work
 *  with both BCF and VCF.
 *
 ***********************************************************************/

/** These macros are defined only for consistency with other parts of htslib */
alias bcf_init1 = bcf_init;
alias bcf_read1 = bcf_read;
alias vcf_read1 = vcf_read;
alias bcf_write1 = bcf_write;
alias vcf_write1 = vcf_write;
alias bcf_destroy1 = bcf_destroy;
alias bcf_empty1 = bcf_empty;
alias vcf_parse1 = vcf_parse;
alias bcf_clear1 = bcf_clear;
alias vcf_format1 = vcf_format;

/**
 *  bcf_hdr_init() - create an empty BCF header.
 *  @param mode    "r" or "w"
 *
 *  When opened for writing, the mandatory fileFormat and
 *  FILTER=PASS lines are added automatically.
 *
 * The bcf_hdr_t struct returned by a successful call should be freed
 * via bcf_hdr_destroy() when it is no longer needed.
 */
bcf_hdr_t* bcf_hdr_init(const(char)* mode);

/** Destroy a BCF header struct */
void bcf_hdr_destroy(bcf_hdr_t* h);

/** Allocate and initialize a bcf1_t object.
 *
 * The bcf1_t struct returned by a successful call should be freed
 * via bcf_destroy() when it is no longer needed.
 */
bcf1_t* bcf_init();

/** Deallocate a bcf1_t object */
void bcf_destroy(bcf1_t* v);

/**
 *  Same as bcf_destroy() but frees only the memory allocated by bcf1_t,
 *  not the bcf1_t object itself.
 */
void bcf_empty(bcf1_t* v);

/**
 *  Make the bcf1_t object ready for next read. Intended mostly for
 *  internal use, the user should rarely need to call this function
 *  directly.
 */
void bcf_clear(bcf1_t* v);

/** bcf_open and vcf_open mode: please see hts_open() in hts.h */
alias vcfFile = htsFile_;
alias bcf_open = hts_open;
alias vcf_open = hts_open;
alias bcf_close = hts_close;
alias vcf_close = hts_close;

/// Read a VCF or BCF header
/** @param  fp  The file to read the header from
    @return Pointer to a populated header structure on success;
            NULL on failure

    The bcf_hdr_t struct returned by a successful call should be freed
    via bcf_hdr_destroy() when it is no longer needed.
*/
bcf_hdr_t* bcf_hdr_read(htsFile* fp);

/**
 *  bcf_hdr_set_samples() - for more efficient VCF parsing when only one/few samples are needed
 *  @param samples  samples to include or exclude from file or as a comma-separated string.
 *              LIST|FILE   .. select samples in list/file
 *              ^LIST|FILE  .. exclude samples from list/file
 *              -           .. include all samples
 *              NULL        .. exclude all samples
 *  @param is_file  @p samples is a file (1) or a comma-separated list (0)
 *
 *  The bottleneck of VCF reading is parsing of genotype fields. If the
 *  reader knows in advance that only subset of samples is needed (possibly
 *  no samples at all), the performance of bcf_read() can be significantly
 *  improved by calling bcf_hdr_set_samples after bcf_hdr_read().
 *  The function bcf_read() will subset the VCF/BCF records automatically
 *  with the notable exception when reading records via bcf_itr_next().
 *  In this case, bcf_subset_format() must be called explicitly, because
 *  bcf_readrec() does not see the header.
 *
 *  Returns 0 on success, -1 on error or a positive integer if the list
 *  contains samples not present in the VCF header. In such a case, the
 *  return value is the index of the offending sample.
 */
int bcf_hdr_set_samples(bcf_hdr_t* hdr, const(char)* samples, int is_file);

int bcf_subset_format(const(bcf_hdr_t)* hdr, bcf1_t* rec);

/// Write a VCF or BCF header
/** @param  fp  Output file
    @param  h   The header to write
    @return 0 on success; -1 on failure
 */
int bcf_hdr_write(htsFile* fp, bcf_hdr_t* h);

/**
 * Parse VCF line contained in kstring and populate the bcf1_t struct
 * The line must not end with \n or \r characters.
 */
int vcf_parse(kstring_t* s, const(bcf_hdr_t)* h, bcf1_t* v);

/**
 * Complete the file opening mode, according to its extension.
 * @param mode      Preallocated mode string to be completed.
 * @param fn        File name to be opened.
 * @param format    Format string (vcf|bcf|vcf.gz)
 * @return          0 on success; -1 on failure
 */
int vcf_open_mode(char* mode, const(char)* fn, const(char)* format);

/** The opposite of vcf_parse. It should rarely be called directly, see vcf_write */
int vcf_format(const(bcf_hdr_t)* h, const(bcf1_t)* v, kstring_t* s);

/// Read next VCF or BCF record
/** @param fp  The file to read the record from
        @param h   The header for the vcf/bcf file
        @param v   The bcf1_t structure to populate
        @return 0 on success; -1 on end of file; < -1 on critical error

On errors which are not critical for reading, such as missing header
definitions in vcf files, zero will be returned but v->errcode will have been
set to one of BCF_ERR* codes and must be checked before calling bcf_write().
     */
int bcf_read(htsFile* fp, const(bcf_hdr_t)* h, bcf1_t* v);

/**
 *  bcf_unpack() - unpack/decode a BCF record (fills the bcf1_t::d field)
 *
 *  Note that bcf_unpack() must be called even when reading VCF. It is safe
 *  to call the function repeatedly, it will not unpack the same field
 *  twice.
 */
enum BCF_UN_STR = 1; // up to ALT inclusive
enum BCF_UN_FLT = 2; // up to FILTER
enum BCF_UN_INFO = 4; // up to INFO
enum BCF_UN_SHR = BCF_UN_STR | BCF_UN_FLT | BCF_UN_INFO; // all shared information
enum BCF_UN_FMT = 8; // unpack format and each sample
enum BCF_UN_IND = BCF_UN_FMT; // a synonym of BCF_UN_FMT
enum BCF_UN_ALL = BCF_UN_SHR | BCF_UN_FMT; // everything
int bcf_unpack(bcf1_t* b, int which);

/*
 *  bcf_dup() - create a copy of BCF record.
 *
 *  Note that bcf_unpack() must be called on the returned copy as if it was
 *  obtained from bcf_read(). Also note that bcf_dup() calls bcf_sync1(src)
 *  internally to reflect any changes made by bcf_update_* functions.
 *
 *  The bcf1_t struct returned by a successful call should be freed
 *  via bcf_destroy() when it is no longer needed.
 */
bcf1_t* bcf_dup(bcf1_t* src);

bcf1_t* bcf_copy(bcf1_t* dst, bcf1_t* src);

/// Write one VCF or BCF record. The type is determined at the open() call.
/** @param  fp  The file to write to
    @param  h   The header for the vcf/bcf file
    @param  v   The bcf1_t structure to write
    @return 0 on success; -1 on error
 */
int bcf_write(htsFile* fp, bcf_hdr_t* h, bcf1_t* v);

/**
 *  The following functions work only with VCFs and should rarely be called
 *  directly. Usually one wants to use their bcf_* alternatives, which work
 *  transparently with both VCFs and BCFs.
 */
/// Read a VCF format header
/** @param  fp  The file to read the header from
    @return Pointer to a populated header structure on success;
            NULL on failure

    Use bcf_hdr_read() instead.

    The bcf_hdr_t struct returned by a successful call should be freed
    via bcf_hdr_destroy() when it is no longer needed.
*/
bcf_hdr_t* vcf_hdr_read(htsFile* fp);

/// Write a VCF format header
/** @param  fp  Output file
    @param  h   The header to write
    @return 0 on success; -1 on failure

    Use bcf_hdr_write() instead
*/
int vcf_hdr_write(htsFile* fp, const(bcf_hdr_t)* h);

/// Read a record from a VCF file
/** @param fp  The file to read the record from
    @param h   The header for the vcf file
    @param v   The bcf1_t structure to populate
    @return 0 on success; -1 on end of file; < -1 on error

    Use bcf_read() instead
*/
int vcf_read(htsFile* fp, const(bcf_hdr_t)* h, bcf1_t* v);

/// Write a record to a VCF file
/** @param  fp  The file to write to
    @param h   The header for the vcf file
    @param v   The bcf1_t structure to write
    @return 0 on success; -1 on error

    Use bcf_write() instead
*/
int vcf_write(htsFile* fp, const(bcf_hdr_t)* h, bcf1_t* v);

/** Helper function for the bcf_itr_next() macro; internal use, ignore it */
int bcf_readrec(
    BGZF* fp,
    void* null_,
    void* v,
    int* tid,
    hts_pos_t* beg,
    hts_pos_t* end);

/// Write a line to a VCF file
/** @param line   Line to write
    @param fp     File to write it to
    @return 0 on success; -1 on failure

    @note No checks are done on the line being added, apart from
          ensuring that it ends with a newline.  This function
          should therefore be used with care.
*/
int vcf_write_line(htsFile* fp, kstring_t* line);

/**************************************************************************
 *  Header querying and manipulation routines
 **************************************************************************/

/** Create a new header using the supplied template
 *
 *  The bcf_hdr_t struct returned by a successful call should be freed
 *  via bcf_hdr_destroy() when it is no longer needed.
 *  @return NULL on failure, header otherwise
 */
bcf_hdr_t* bcf_hdr_dup(const(bcf_hdr_t)* hdr);

/**
 *  Copy header lines from src to dst if not already present in dst. See also bcf_translate().
 *  Returns 0 on success or sets a bit on error:
 *      1 .. conflicting definitions of tag length
 *      // todo
 */
int bcf_hdr_combine(bcf_hdr_t* dst, const(bcf_hdr_t)* src);

/**
 *  bcf_hdr_merge() - copy header lines from src to dst, see also bcf_translate()
 *  @param dst: the destination header to be merged into, NULL on the first pass
 *  @param src: the source header
 *  @return NULL on failure, header otherwise
 *
 *  Notes:
 *      - use as:
 *          bcf_hdr_t *dst = NULL;
 *          for (i=0; i<nsrc; i++) dst = bcf_hdr_merge(dst,src[i]);
 *
 *      - bcf_hdr_merge() replaces bcf_hdr_combine() which had a problem when
 *      combining multiple BCF headers. The current bcf_hdr_combine()
 *      does not have this problem, but became slow when used for many files.
 */
bcf_hdr_t* bcf_hdr_merge(bcf_hdr_t* dst, const(bcf_hdr_t)* src);

/**
 *  bcf_hdr_add_sample() - add a new sample.
 *  @param sample:  sample name to be added
 *
 *  Note:
 *      After all samples have been added, the internal header structure must be updated
 *      by calling bcf_hdr_sync(). This is normally done automatically by the first bcf_hdr_write()
 *      or bcf_write() call. Otherwise, the caller must force the update by calling bcf_hdr_sync()
 *      explicitly.
 */
int bcf_hdr_add_sample(bcf_hdr_t* hdr, const(char)* sample);

/** Read VCF header from a file and update the header */
int bcf_hdr_set(bcf_hdr_t* hdr, const(char)* fname);

/// Appends formatted header text to _str_.
/** If _is_bcf_ is zero, `IDX` fields are discarded.
 *  @return 0 if successful, or negative if an error occurred
 *  @since 1.4
 */
int bcf_hdr_format(const(bcf_hdr_t)* hdr, int is_bcf, kstring_t* str);

/** Returns formatted header (newly allocated string) and its length,
 *  excluding the terminating \0. If is_bcf parameter is unset, IDX
 *  fields are discarded.
 *  @deprecated Use bcf_hdr_format() instead as it can handle huge headers.
 */
char* bcf_hdr_fmt_text(const(bcf_hdr_t)* hdr, int is_bcf, int* len);

/** Append new VCF header line, returns 0 on success */
int bcf_hdr_append(bcf_hdr_t* h, const(char)* line);

int bcf_hdr_printf(bcf_hdr_t* h, const(char)* format, ...);

/** VCF version, e.g. VCFv4.2 */
const(char)* bcf_hdr_get_version(const(bcf_hdr_t)* hdr);

/// Set version in bcf header
/**
   @param hdr     BCF header struct
   @param version Version to set, e.g. "VCFv4.3"
   @return 0 on success; < 0 on error
 */
int bcf_hdr_set_version(bcf_hdr_t* hdr, const(char)* version_);

/**
 *  bcf_hdr_remove() - remove VCF header tag
 *  @param type:      one of BCF_HL_*
 *  @param key:       tag name or NULL to remove all tags of the given type
 */
void bcf_hdr_remove(bcf_hdr_t* h, int type, const(char)* key);

/**
 *  bcf_hdr_subset() - creates a new copy of the header removing unwanted samples
 *  @param n:        number of samples to keep
 *  @param samples:  names of the samples to keep
 *  @param imap:     mapping from index in @samples to the sample index in the original file
 *  @return NULL on failure, header otherwise
 *
 *  Sample names not present in h0 are ignored. The number of unmatched samples can be checked
 *  by comparing n and bcf_hdr_nsamples(out_hdr).
 *  This function can be used to reorder samples.
 *  See also bcf_subset() which subsets individual records.
 *  The bcf_hdr_t struct returned by a successful call should be freed
 *  via bcf_hdr_destroy() when it is no longer needed.
 */
bcf_hdr_t* bcf_hdr_subset(
    const(bcf_hdr_t)* h0,
    int n,
    char** samples,
    int* imap);

/** Creates a list of sequence names. It is up to the caller to free the list (but not the sequence names) */
const(char*)* bcf_hdr_seqnames(const(bcf_hdr_t)* h, int* nseqs);

/** Get number of samples */
extern (D) auto bcf_hdr_nsamples(T)(auto ref T hdr)
{
    return hdr.n[BCF_DT_SAMPLE];
}

/** The following functions are for internal use and should rarely be called directly */
int bcf_hdr_parse(bcf_hdr_t* hdr, char* htxt);

/// Synchronize internal header structures
/** @param h  Header
    @return 0 on success, -1 on failure

    This function updates the id, sample and contig arrays in the
    bcf_hdr_t structure so that they point to the same locations as
    the id, sample and contig dictionaries.
*/
int bcf_hdr_sync(bcf_hdr_t* h);

/**
 * bcf_hdr_parse_line() - parse a single line of VCF textual header
 * @param h     BCF header struct
 * @param line  One or more lines of header text
 * @param len   Filled out with length data parsed from 'line'.
 * @return bcf_hrec_t* on success;
 *         NULL on error or on end of header text.
 *         NB: to distinguish error from end-of-header, check *len:
 *           *len == 0 indicates @p line did not start with "##"
 *           *len == -1 indicates failure, likely due to out of memory
 *           *len > 0 indicates a malformed header line
 *
 * If *len > 0 on exit, it will contain the full length of the line
 * including any trailing newline (this includes cases where NULL was
 * returned due to a malformed line).  Callers can use this to skip to
 * the next header line.
 */
bcf_hrec_t* bcf_hdr_parse_line(
    const(bcf_hdr_t)* h,
    const(char)* line,
    int* len);
/// Convert a bcf header record to string form
/**
 * @param hrec    Header record
 * @param str     Destination kstring
 * @return 0 on success; < 0 on error
 */
int bcf_hrec_format(const(bcf_hrec_t)* hrec, kstring_t* str);

int bcf_hdr_add_hrec(bcf_hdr_t* hdr, bcf_hrec_t* hrec);

/**
 *  bcf_hdr_get_hrec() - get header line info
 *  @param type:  one of the BCF_HL_* types: FLT,INFO,FMT,CTG,STR,GEN
 *  @param key:   the header key for generic lines (e.g. "fileformat"), any field
 *                  for structured lines, typically "ID".
 *  @param value: the value which pairs with key. Can be be NULL for BCF_HL_GEN
 *  @param str_class: the class of BCF_HL_STR line (e.g. "ALT" or "SAMPLE"), otherwise NULL
 */
bcf_hrec_t* bcf_hdr_get_hrec(
    const(bcf_hdr_t)* hdr,
    int type,
    const(char)* key,
    const(char)* value,
    const(char)* str_class);

/// Duplicate a header record
/** @param hrec   Header record to copy
    @return A new header record on success; NULL on failure

    The bcf_hrec_t struct returned by a successful call should be freed
    via bcf_hrec_destroy() when it is no longer needed.
*/
bcf_hrec_t* bcf_hrec_dup(bcf_hrec_t* hrec);

/// Add a new header record key
/** @param hrec  Header record
    @param str   Key name
    @param len   Length of @p str
    @return 0 on success; -1 on failure
*/
int bcf_hrec_add_key(bcf_hrec_t* hrec, const(char)* str, size_t len);

/// Set a header record value
/** @param hrec      Header record
    @param i         Index of value
    @param str       Value to set
    @param len       Length of @p str
    @param is_quoted Value should be quoted
    @return 0 on success; -1 on failure
*/
int bcf_hrec_set_val(
    bcf_hrec_t* hrec,
    int i,
    const(char)* str,
    size_t len,
    int is_quoted);

int bcf_hrec_find_key(bcf_hrec_t* hrec, const(char)* key);

/// Add an IDX header record
/** @param hrec   Header record
    @param idx    IDX value to add
    @return 0 on success; -1 on failure
*/
int hrec_add_idx(bcf_hrec_t* hrec, int idx);

/// Free up a header record and associated structures
/** @param hrec  Header record
 */
void bcf_hrec_destroy(bcf_hrec_t* hrec);

/**************************************************************************
 *  Individual record querying and manipulation routines
 **************************************************************************/

/** See the description of bcf_hdr_subset() */
int bcf_subset(const(bcf_hdr_t)* h, bcf1_t* v, int n, int* imap);

/**
 *  bcf_translate() - translate tags ids to be consistent with different header. This function
 *                    is useful when lines from multiple VCF need to be combined.
 *  @dst_hdr:   the destination header, to be used in bcf_write(), see also bcf_hdr_combine()
 *  @src_hdr:   the source header, used in bcf_read()
 *  @src_line:  line obtained by bcf_read()
 */
int bcf_translate(
    const(bcf_hdr_t)* dst_hdr,
    bcf_hdr_t* src_hdr,
    bcf1_t* src_line);

/**
 *  bcf_get_variant_type[s]()  - returns one of VCF_REF, VCF_SNP, etc
 */
int bcf_get_variant_types(bcf1_t* rec);

int bcf_get_variant_type(bcf1_t* rec, int ith_allele);

int bcf_is_snp(bcf1_t* v);

/**
 *  bcf_update_filter() - sets the FILTER column
 *  @flt_ids:  The filter IDs to set, numeric IDs returned by bcf_hdr_id2int(hdr, BCF_DT_ID, "PASS")
 *  @n:        Number of filters. If n==0, all filters are removed
 */
int bcf_update_filter(const(bcf_hdr_t)* hdr, bcf1_t* line, int* flt_ids, int n);
/**
 *  bcf_add_filter() - adds to the FILTER column
 *  @flt_id:   filter ID to add, numeric ID returned by bcf_hdr_id2int(hdr, BCF_DT_ID, "PASS")
 *
 *  If flt_id is PASS, all existing filters are removed first. If other than PASS, existing PASS is removed.
 */
int bcf_add_filter(const(bcf_hdr_t)* hdr, bcf1_t* line, int flt_id);
/**
 *  bcf_remove_filter() - removes from the FILTER column
 *  @flt_id:   filter ID to remove, numeric ID returned by bcf_hdr_id2int(hdr, BCF_DT_ID, "PASS")
 *  @pass:     when set to 1 and no filters are present, set to PASS
 */
int bcf_remove_filter(
    const(bcf_hdr_t)* hdr,
    bcf1_t* line,
    int flt_id,
    int pass);
/**
 *  Returns 1 if present, 0 if absent, or -1 if filter does not exist. "PASS" and "." can be used interchangeably.
 */
int bcf_has_filter(const(bcf_hdr_t)* hdr, bcf1_t* line, char* filter);
/**
 *  bcf_update_alleles() and bcf_update_alleles_str() - update REF and ALT column
 *  @alleles:           Array of alleles
 *  @nals:              Number of alleles
 *  @alleles_string:    Comma-separated alleles, starting with the REF allele
 */
int bcf_update_alleles(
    const(bcf_hdr_t)* hdr,
    bcf1_t* line,
    const(char*)* alleles,
    int nals);

int bcf_update_alleles_str(
    const(bcf_hdr_t)* hdr,
    bcf1_t* line,
    const(char)* alleles_string);

/**
  *  bcf_update_id() - sets new ID string
  *  bcf_add_id() - adds to the ID string checking for duplicates
  */
int bcf_update_id(const(bcf_hdr_t)* hdr, bcf1_t* line, const(char)* id);

int bcf_add_id(const(bcf_hdr_t)* hdr, bcf1_t* line, const(char)* id);

/**
 *  bcf_update_info_*() - functions for updating INFO fields
 *  @param hdr:       the BCF header
 *  @param line:      VCF line to be edited
 *  @param key:       the INFO tag to be updated
 *  @param values:    pointer to the array of values. Pass NULL to remove the tag.
 *  @param n:         number of values in the array. When set to 0, the INFO tag is removed
 *  @return 0 on success or negative value on error.
 *
 *  The @p string in bcf_update_info_flag() is optional,
 *  @p n indicates whether the flag is set or removed.
 *
 *  Note that updating an END info tag will cause line->rlen to be
 *  updated as a side-effect (removing the tag will set it to the
 *  string length of the REF allele). If line->pos is being changed as
 *  well, it is important that this is done before calling
 *  bcf_update_info_int32() to update the END tag, otherwise rlen will be
 *  set incorrectly.  If the new END value is less than or equal to
 *  line->pos, a warning will be printed and line->rlen will be set to
 *  the length of the REF allele.
 */
extern (D) auto bcf_update_info_int32(T0, T1, T2, T3, T4)(auto ref T0 hdr, auto ref T1 line, auto ref T2 key, auto ref T3 values, auto ref T4 n)
{
    return bcf_update_info(hdr, line, key, values, n, BCF_HT_INT);
}

extern (D) auto bcf_update_info_float(T0, T1, T2, T3, T4)(auto ref T0 hdr, auto ref T1 line, auto ref T2 key, auto ref T3 values, auto ref T4 n)
{
    return bcf_update_info(hdr, line, key, values, n, BCF_HT_REAL);
}

extern (D) auto bcf_update_info_flag(T0, T1, T2, T3, T4)(auto ref T0 hdr, auto ref T1 line, auto ref T2 key, auto ref T3 string, auto ref T4 n)
{
    return bcf_update_info(hdr, line, key, string, n, BCF_HT_FLAG);
}

extern (D) auto bcf_update_info_string(T0, T1, T2, T3)(auto ref T0 hdr, auto ref T1 line, auto ref T2 key, auto ref T3 string)
{
    return bcf_update_info(hdr, line, key, string, 1, BCF_HT_STR);
}

int bcf_update_info(
    const(bcf_hdr_t)* hdr,
    bcf1_t* line,
    const(char)* key,
    const(void)* values,
    int n,
    int type);

/// Set or update 64-bit integer INFO values
/**
 *  @param hdr:       the BCF header
 *  @param line:      VCF line to be edited
 *  @param key:       the INFO tag to be updated
 *  @param values:    pointer to the array of values. Pass NULL to remove the tag.
 *  @param n:         number of values in the array. When set to 0, the INFO tag is removed
 *  @return 0 on success or negative value on error.
 *
 *  This function takes an int64_t values array as input.  The data
 *  actually stored will be shrunk to the minimum size that can
 *  accept all of the values.
 *
 *  INFO values outside of the range BCF_MIN_BT_INT32 to BCF_MAX_BT_INT32
 *  can only be written to VCF files.
 */
int bcf_update_info_int64(
    const(bcf_hdr_t)* hdr,
    bcf1_t* line,
    const(char)* key,
    const(long)* values,
    int n);

/*
 *  bcf_update_format_*() - functions for updating FORMAT fields
 *  @values:    pointer to the array of values, the same number of elements
 *              is expected for each sample. Missing values must be padded
 *              with bcf_*_missing or bcf_*_vector_end values.
 *  @n:         number of values in the array. If n==0, existing tag is removed.
 *
 *  The function bcf_update_format_string() is a higher-level (slower) variant of
 *  bcf_update_format_char(). The former accepts array of \0-terminated strings
 *  whereas the latter requires that the strings are collapsed into a single array
 *  of fixed-length strings. In case of strings with variable length, shorter strings
 *  can be \0-padded. Note that the collapsed strings passed to bcf_update_format_char()
 *  are not \0-terminated.
 *
 *  Returns 0 on success or negative value on error.
 */
extern (D) auto bcf_update_format_int32(T0, T1, T2, T3, T4)(auto ref T0 hdr, auto ref T1 line, auto ref T2 key, auto ref T3 values, auto ref T4 n)
{
    return bcf_update_format(hdr, line, key, values, n, BCF_HT_INT);
}

extern (D) auto bcf_update_format_float(T0, T1, T2, T3, T4)(auto ref T0 hdr, auto ref T1 line, auto ref T2 key, auto ref T3 values, auto ref T4 n)
{
    return bcf_update_format(hdr, line, key, values, n, BCF_HT_REAL);
}

extern (D) auto bcf_update_format_char(T0, T1, T2, T3, T4)(auto ref T0 hdr, auto ref T1 line, auto ref T2 key, auto ref T3 values, auto ref T4 n)
{
    return bcf_update_format(hdr, line, key, values, n, BCF_HT_STR);
}

extern (D) auto bcf_update_genotypes(T0, T1, T2, T3)(auto ref T0 hdr, auto ref T1 line, auto ref T2 gts, auto ref T3 n)
{
    return bcf_update_format(hdr, line, "GT", gts, n, BCF_HT_INT);
} // See bcf_gt_ macros below

int bcf_update_format_string(
    const(bcf_hdr_t)* hdr,
    bcf1_t* line,
    const(char)* key,
    const(char*)* values,
    int n);

int bcf_update_format(
    const(bcf_hdr_t)* hdr,
    bcf1_t* line,
    const(char)* key,
    const(void)* values,
    int n,
    int type);

// Macros for setting genotypes correctly, for use with bcf_update_genotypes only; idx corresponds
// to VCF's GT (1-based index to ALT or 0 for the reference allele) and val is the opposite, obtained
// from bcf_get_genotypes() below.
extern (D) auto bcf_gt_phased(T)(auto ref T idx)
{
    return (idx + 1) << 1 | 1;
}

extern (D) auto bcf_gt_unphased(T)(auto ref T idx)
{
    return (idx + 1) << 1;
}

enum bcf_gt_missing = 0;

extern (D) int bcf_gt_is_missing(T)(auto ref T val)
{
    return val >> 1 ? 0 : 1;
}

extern (D) auto bcf_gt_is_phased(T)(auto ref T idx)
{
    return idx & 1;
}

extern (D) auto bcf_gt_allele(T)(auto ref T val)
{
    return (val >> 1) - 1;
}

/** Conversion between alleles indexes to Number=G genotype index (assuming diploid, all 0-based) */
extern (D) auto bcf_alleles2gt(T0, T1)(auto ref T0 a, auto ref T1 b)
{
    return a > b ? (a * (a + 1) / 2 + b) : (b * (b + 1) / 2 + a);
}

void bcf_gt2alleles(int igt, int* a, int* b);

/**
 * bcf_get_fmt() - returns pointer to FORMAT's field data
 * @header: for access to BCF_DT_ID dictionary
 * @line:   VCF line obtained from vcf_parse1
 * @fmt:    one of GT,PL,...
 *
 * Returns bcf_fmt_t* if the call succeeded, or returns NULL when the field
 * is not available.
 */
bcf_fmt_t* bcf_get_fmt(const(bcf_hdr_t)* hdr, bcf1_t* line, const(char)* key);

bcf_info_t* bcf_get_info(const(bcf_hdr_t)* hdr, bcf1_t* line, const(char)* key);

/**
 * bcf_get_*_id() - returns pointer to FORMAT/INFO field data given the header index instead of the string ID
 * @line: VCF line obtained from vcf_parse1
 * @id:  The header index for the tag, obtained from bcf_hdr_id2int()
 *
 * Returns bcf_fmt_t* / bcf_info_t*. These functions do not check if the index is valid
 * as their goal is to avoid the header lookup.
 */
bcf_fmt_t* bcf_get_fmt_id(bcf1_t* line, const int id);

bcf_info_t* bcf_get_info_id(bcf1_t* line, const int id);

/**
 *  bcf_get_info_*() - get INFO values, integers or floats
 *  @param hdr:    BCF header
 *  @param line:   BCF record
 *  @param tag:    INFO tag to retrieve
 *  @param dst:    *dst is pointer to a memory location, can point to NULL
 *  @param ndst:   pointer to the size of allocated memory
 *  @return  >=0 on success
 *          -1 .. no such INFO tag defined in the header
 *          -2 .. clash between types defined in the header and encountered in the VCF record
 *          -3 .. tag is not present in the VCF record
 *          -4 .. the operation could not be completed (e.g. out of memory)
 *
 *  Returns negative value on error or the number of values (including
 *  missing values) put in *dst on success. bcf_get_info_string() returns
 *  on success the number of characters stored excluding the nul-
 *  terminating byte. bcf_get_info_flag() does not store anything in *dst
 *  but returns 1 if the flag is set or 0 if not.
 *
 *  *dst will be reallocated if it is not big enough (i.e. *ndst is too
 *  small) or NULL on entry.  The new size will be stored in *ndst.
 */
extern (D) auto bcf_get_info_int32(T0, T1, T2, T3, T4)(auto ref T0 hdr, auto ref T1 line, auto ref T2 tag, auto ref T3 dst, auto ref T4 ndst)
{
    return bcf_get_info_values(hdr, line, tag, cast(void**) dst, ndst, BCF_HT_INT);
}

extern (D) auto bcf_get_info_float(T0, T1, T2, T3, T4)(auto ref T0 hdr, auto ref T1 line, auto ref T2 tag, auto ref T3 dst, auto ref T4 ndst)
{
    return bcf_get_info_values(hdr, line, tag, cast(void**) dst, ndst, BCF_HT_REAL);
}

extern (D) auto bcf_get_info_string(T0, T1, T2, T3, T4)(auto ref T0 hdr, auto ref T1 line, auto ref T2 tag, auto ref T3 dst, auto ref T4 ndst)
{
    return bcf_get_info_values(hdr, line, tag, cast(void**) dst, ndst, BCF_HT_STR);
}

extern (D) auto bcf_get_info_flag(T0, T1, T2, T3, T4)(auto ref T0 hdr, auto ref T1 line, auto ref T2 tag, auto ref T3 dst, auto ref T4 ndst)
{
    return bcf_get_info_values(hdr, line, tag, cast(void**) dst, ndst, BCF_HT_FLAG);
}

int bcf_get_info_values(
    const(bcf_hdr_t)* hdr,
    bcf1_t* line,
    const(char)* tag,
    void** dst,
    int* ndst,
    int type);

/// Put integer INFO values into an int64_t array
/**
 *  @param hdr:    BCF header
 *  @param line:   BCF record
 *  @param tag:    INFO tag to retrieve
 *  @param dst:    *dst is pointer to a memory location, can point to NULL
 *  @param ndst:   pointer to the size of allocated memory
 *  @return  >=0 on success
 *          -1 .. no such INFO tag defined in the header
 *          -2 .. clash between types defined in the header and encountered in the VCF record
 *          -3 .. tag is not present in the VCF record
 *          -4 .. the operation could not be completed (e.g. out of memory)
 *
 *  Returns negative value on error or the number of values (including
 *  missing values) put in *dst on success.
 *
 *  *dst will be reallocated if it is not big enough (i.e. *ndst is too
 *  small) or NULL on entry.  The new size will be stored in *ndst.
 */
int bcf_get_info_int64(
    const(bcf_hdr_t)* hdr,
    bcf1_t* line,
    const(char)* tag,
    long** dst,
    int* ndst);

/**
 *  bcf_get_format_*() - same as bcf_get_info*() above
 *
 *  The function bcf_get_format_string() is a higher-level (slower) variant of bcf_get_format_char().
 *  see the description of bcf_update_format_string() and bcf_update_format_char() above.
 *  Unlike other bcf_get_format__*() functions, bcf_get_format_string() allocates two arrays:
 *  a single block of \0-terminated strings collapsed into a single array and an array of pointers
 *  to these strings. Both arrays must be cleaned by the user.
 *
 *  Returns negative value on error or the number of written values on success.
 *
 *  Use the returned number of written values for accessing valid entries of dst, as ndst is only a
 *  watermark that can be higher than the returned value, i.e. the end of dst can contain carry-over
 *  values from previous calls to bcf_get_format_*() on lines with more values per sample.
 *
 *  Example:
 *      int ndst = 0; char **dst = NULL;
 *      if ( bcf_get_format_string(hdr, line, "XX", &dst, &ndst) > 0 )
 *          for (i=0; i<bcf_hdr_nsamples(hdr); i++) printf("%s\n", dst[i]);
 *      free(dst[0]); free(dst);
 *
 *  Example:
 *      int i, j, ngt, nsmpl = bcf_hdr_nsamples(hdr);
 *      int32_t *gt_arr = NULL, ngt_arr = 0;
 *
 *      ngt = bcf_get_genotypes(hdr, line, &gt_arr, &ngt_arr);
 *      if ( ngt<=0 ) return; // GT not present
 *
 *      int max_ploidy = ngt/nsmpl;
 *      for (i=0; i<nsmpl; i++)
 *      {
 *        int32_t *ptr = gt_arr + i*max_ploidy;
 *        for (j=0; j<max_ploidy; j++)
 *        {
 *           // if true, the sample has smaller ploidy
 *           if ( ptr[j]==bcf_int32_vector_end ) break;
 *
 *           // missing allele
 *           if ( bcf_gt_is_missing(ptr[j]) ) continue;
 *
 *           // the VCF 0-based allele index
 *           int allele_index = bcf_gt_allele(ptr[j]);
 *
 *           // is phased?
 *           int is_phased = bcf_gt_is_phased(ptr[j]);
 *
 *           // .. do something ..
 *         }
 *      }
 *      free(gt_arr);
 *
 */
extern (D) auto bcf_get_format_int32(T0, T1, T2, T3, T4)(auto ref T0 hdr, auto ref T1 line, auto ref T2 tag, auto ref T3 dst, auto ref T4 ndst)
{
    return bcf_get_format_values(hdr, line, tag, cast(void**) dst, ndst, BCF_HT_INT);
}

extern (D) auto bcf_get_format_float(T0, T1, T2, T3, T4)(auto ref T0 hdr, auto ref T1 line, auto ref T2 tag, auto ref T3 dst, auto ref T4 ndst)
{
    return bcf_get_format_values(hdr, line, tag, cast(void**) dst, ndst, BCF_HT_REAL);
}

extern (D) auto bcf_get_format_char(T0, T1, T2, T3, T4)(auto ref T0 hdr, auto ref T1 line, auto ref T2 tag, auto ref T3 dst, auto ref T4 ndst)
{
    return bcf_get_format_values(hdr, line, tag, cast(void**) dst, ndst, BCF_HT_STR);
}

extern (D) auto bcf_get_genotypes(T0, T1, T2, T3)(auto ref T0 hdr, auto ref T1 line, auto ref T2 dst, auto ref T3 ndst)
{
    return bcf_get_format_values(hdr, line, "GT", cast(void**) dst, ndst, BCF_HT_INT);
}

int bcf_get_format_string(
    const(bcf_hdr_t)* hdr,
    bcf1_t* line,
    const(char)* tag,
    char*** dst,
    int* ndst);

int bcf_get_format_values(
    const(bcf_hdr_t)* hdr,
    bcf1_t* line,
    const(char)* tag,
    void** dst,
    int* ndst,
    int type);

/**************************************************************************
 *  Helper functions
 **************************************************************************/

/**
 *  bcf_hdr_id2int() - Translates string into numeric ID
 *  bcf_hdr_int2id() - Translates numeric ID into string
 *  @type:     one of BCF_DT_ID, BCF_DT_CTG, BCF_DT_SAMPLE
 *  @id:       tag name, such as: PL, DP, GT, etc.
 *
 *  Returns -1 if string is not in dictionary, otherwise numeric ID which identifies
 *  fields in BCF records.
 */
int bcf_hdr_id2int(const(bcf_hdr_t)* hdr, int type, const(char)* id);

extern (D) auto bcf_hdr_int2id(T0, T1, T2)(auto ref T0 hdr, auto ref T1 type, auto ref T2 int_id)
{
    return hdr.id[type][int_id].key;
}

/**
 *  bcf_hdr_name2id() - Translates sequence names (chromosomes) into numeric ID
 *  bcf_hdr_id2name() - Translates numeric ID to sequence name
 */
int bcf_hdr_name2id(const(bcf_hdr_t)* hdr, const(char)* id);
const(char)* bcf_hdr_id2name(const(bcf_hdr_t)* hdr, int rid);
const(char)* bcf_seqname(const(bcf_hdr_t)* hdr, const(bcf1_t)* rec);

/** Return CONTIG name, or "(unknown)"

    Like bcf_seqname(), but this function will never return NULL.  If
    the contig name cannot be found (either because @p hdr was not
    supplied or rec->rid was out of range) it returns the string
    "(unknown)".
*/
const(char)* bcf_seqname_safe(const(bcf_hdr_t)* hdr, const(bcf1_t)* rec);

/**
 *  bcf_hdr_id2*() - Macros for accessing bcf_idinfo_t
 *  @type:      one of BCF_HL_FLT, BCF_HL_INFO, BCF_HL_FMT
 *  @int_id:    return value of bcf_hdr_id2int, must be >=0
 *
 *  The returned values are:
 *     bcf_hdr_id2length   ..  whether the number of values is fixed or variable, one of BCF_VL_*
 *     bcf_hdr_id2number   ..  the number of values, 0xfffff for variable length fields
 *     bcf_hdr_id2type     ..  the field type, one of BCF_HT_*
 *     bcf_hdr_id2coltype  ..  the column type, one of BCF_HL_*
 *
 *  Notes: Prior to using the macros, the presence of the info should be
 *  tested with bcf_hdr_idinfo_exists().
 */
extern (D) auto bcf_hdr_id2length(T0, T1, T2)(auto ref T0 hdr, auto ref T1 type, auto ref T2 int_id)
{
    return hdr.id[BCF_DT_ID][int_id].val.info[type] >> 8 & 0xf;
}

extern (D) auto bcf_hdr_id2number(T0, T1, T2)(auto ref T0 hdr, auto ref T1 type, auto ref T2 int_id)
{
    return hdr.id[BCF_DT_ID][int_id].val.info[type] >> 12;
}

extern (D) auto bcf_hdr_id2type(T0, T1, T2)(auto ref T0 hdr, auto ref T1 type, auto ref T2 int_id)
{
    return cast(uint) hdr.id[BCF_DT_ID][int_id].val.info[type] >> 4 & 0xf;
}

extern (D) auto bcf_hdr_id2coltype(T0, T1, T2)(auto ref T0 hdr, auto ref T1 type, auto ref T2 int_id)
{
    return cast(uint) hdr.id[BCF_DT_ID][int_id].val.info[type] & 0xf;
}

extern (D) auto bcf_hdr_idinfo_exists(T0, T1, T2)(auto ref T0 hdr, auto ref T1 type, auto ref T2 int_id)
{
    return int_id >= 0 && bcf_hdr_id2coltype(hdr, type, int_id) != 0xf;
}

extern (D) auto bcf_hdr_id2hrec(T0, T1, T2, T3)(auto ref T0 hdr, auto ref T1 dict_type, auto ref T2 col_type, auto ref T3 int_id)
{
    return hdr.id[dict_type == BCF_DT_CTG ? BCF_DT_CTG : BCF_DT_ID][int_id].val.hrec[dict_type == BCF_DT_CTG ? 0 : col_type];
}

/// Convert BCF FORMAT data to string form
/**
 * @param s    kstring to write into
 * @param n    number of items in @p data
 * @param type type of items in @p data
 * @param data BCF format data
 * @return  0 on success
 *         -1 if out of memory
 */
int bcf_fmt_array(kstring_t* s, int n, int type, void* data);

ubyte* bcf_fmt_sized_array(kstring_t* s, ubyte* ptr);

/// Encode a variable-length char array in BCF format
/**
 * @param s    kstring to write into
 * @param l    length of input
 * @param a    input data to encode
 * @return 0 on success; < 0 on error
 */
int bcf_enc_vchar(kstring_t* s, int l, const(char)* a);

/// Encode a variable-length integer array in BCF format
/**
 * @param s      kstring to write into
 * @param n      total number of items in @p a (<= 0 to encode BCF_BT_NULL)
 * @param a      input data to encode
 * @param wsize  vector length (<= 0 is equivalent to @p n)
 * @return 0 on success; < 0 on error
 * @note @p n should be an exact multiple of @p wsize
 */
int bcf_enc_vint(kstring_t* s, int n, int* a, int wsize);

/// Encode a variable-length float array in BCF format
/**
 * @param s      kstring to write into
 * @param n      total number of items in @p a (<= 0 to encode BCF_BT_NULL)
 * @param a      input data to encode
 * @return 0 on success; < 0 on error
 */
int bcf_enc_vfloat(kstring_t* s, int n, float* a);

/**************************************************************************
 *  BCF index
 *
 *  Note that these functions work with BCFs only. See synced_bcf_reader.h
 *  which provides (amongst other things) an API to work transparently with
 *  both indexed BCFs and VCFs.
 **************************************************************************/

alias bcf_itr_destroy = hts_itr_destroy;

extern (D) auto bcf_itr_queryi(T0, T1, T2, T3)(auto ref T0 idx, auto ref T1 tid, auto ref T2 beg, auto ref T3 end)
{
    return hts_itr_query(idx, tid, beg, end, bcf_readrec);
}

extern (D) auto bcf_itr_querys(T0, T1, T2)(auto ref T0 idx, auto ref T1 hdr, auto ref T2 s)
{
    return hts_itr_querys(idx, s, cast(hts_name2id_f) bcf_hdr_name2id, hdr, hts_itr_query, bcf_readrec);
}

int bcf_itr_next(htsFile* htsfp, hts_itr_t* itr, void* r);
/// Load a BCF index
/** @param fn   BCF file name
    @return The index, or NULL if an error occurred.
     @note This only works for BCF files.  Consider synced_bcf_reader instead
which works for both BCF and VCF.
*/
extern (D) auto bcf_index_load(T)(auto ref T fn)
{
    return hts_idx_load(fn, HTS_FMT_CSI);
}

extern (D) auto bcf_index_seqnames(T0, T1, T2)(auto ref T0 idx, auto ref T1 hdr, auto ref T2 nptr)
{
    return hts_idx_seqnames(idx, nptr, cast(hts_id2name_f) bcf_hdr_id2name, hdr);
}

/// Load a BCF index from a given index file name
/**  @param fn     Input BAM/BCF/etc filename
     @param fnidx  The input index filename
     @return  The index, or NULL if an error occurred.
     @note This only works for BCF files.  Consider synced_bcf_reader instead
which works for both BCF and VCF.
*/
hts_idx_t* bcf_index_load2(const(char)* fn, const(char)* fnidx);

/// Load a BCF index from a given index file name
/**  @param fn     Input BAM/BCF/etc filename
     @param fnidx  The input index filename
     @param flags  Flags to alter behaviour (see description)
     @return  The index, or NULL if an error occurred.
     @note This only works for BCF files.  Consider synced_bcf_reader instead
which works for both BCF and VCF.

     The @p flags parameter can be set to a combination of the following
     values:

        HTS_IDX_SAVE_REMOTE   Save a local copy of any remote indexes
        HTS_IDX_SILENT_FAIL   Fail silently if the index is not present

     Equivalent to hts_idx_load3(fn, fnidx, HTS_FMT_CSI, flags);
*/
hts_idx_t* bcf_index_load3(const(char)* fn, const(char)* fnidx, int flags);

/**
 *  bcf_index_build() - Generate and save an index file
 *  @fn:         Input VCF(compressed)/BCF filename
 *  @min_shift:  log2(width of the smallest bin), e.g. a value of 14
 *  imposes a 16k base lower limit on the width of index bins.
 *  Positive to generate CSI, or 0 to generate TBI. However, a small
 *  value of min_shift would create a large index, which would lead to
 *  reduced performance when using the index. A recommended value is 14.
 *  For BCF files, only the CSI index can be generated.
 *
 *  Returns 0 if successful, or negative if an error occurred.
 *
 *  List of error codes:
 *      -1 .. indexing failed
 *      -2 .. opening @fn failed
 *      -3 .. format not indexable
 *      -4 .. failed to create and/or save the index
 */
int bcf_index_build(const(char)* fn, int min_shift);

/**
 *  bcf_index_build2() - Generate and save an index to a specific file
 *  @fn:         Input VCF/BCF filename
 *  @fnidx:      Output filename, or NULL to add .csi/.tbi to @fn
 *  @min_shift:  Positive to generate CSI, or 0 to generate TBI
 *
 *  Returns 0 if successful, or negative if an error occurred.
 *
 *  List of error codes:
 *      -1 .. indexing failed
 *      -2 .. opening @fn failed
 *      -3 .. format not indexable
 *      -4 .. failed to create and/or save the index
 */
int bcf_index_build2(const(char)* fn, const(char)* fnidx, int min_shift);

/**
 *  bcf_index_build3() - Generate and save an index to a specific file
 *  @fn:         Input VCF/BCF filename
 *  @fnidx:      Output filename, or NULL to add .csi/.tbi to @fn
 *  @min_shift:  Positive to generate CSI, or 0 to generate TBI
 *  @n_threads:  Number of VCF/BCF decoder threads
 *
 *  Returns 0 if successful, or negative if an error occurred.
 *
 *  List of error codes:
 *      -1 .. indexing failed
 *      -2 .. opening @fn failed
 *      -3 .. format not indexable
 *      -4 .. failed to create and/or save the index
 */
int bcf_index_build3(
    const(char)* fn,
    const(char)* fnidx,
    int min_shift,
    int n_threads);

/// Initialise fp->idx for the current format type, for VCF and BCF files.
/** @param fp        File handle for the data file being written.
    @param h         BCF header structured (needed for BAI and CSI).
    @param min_shift CSI bin size (CSI default is 14).
    @param fnidx     Filename to write index to.  This pointer must remain valid
                     until after bcf_idx_save is called.
    @return          0 on success, <0 on failure.
    @note This must be called after the header has been written, but before
          any other data.
*/
int bcf_idx_init(htsFile* fp, bcf_hdr_t* h, int min_shift, const(char)* fnidx);

/// Writes the index initialised with bcf_idx_init to disk.
/** @param fp        File handle for the data file being written.
    @return          0 on success, <0 on failure.
*/
int bcf_idx_save(htsFile* fp);

/*******************
 * Typed value I/O *
 *******************/

/*
    Note that in contrast with BCFv2.1 specification, HTSlib implementation
    allows missing values in vectors. For integer types, the values 0x80,
    0x8000, 0x80000000 are interpreted as missing values and 0x81, 0x8001,
    0x80000001 as end-of-vector indicators.  Similarly for floats, the value of
    0x7F800001 is interpreted as a missing value and 0x7F800002 as an
    end-of-vector indicator.
    Note that the end-of-vector byte is not part of the vector.

    This trial BCF version (v2.2) is compatible with the VCF specification and
    enables to handle correctly vectors with different ploidy in presence of
    missing values.
 */
enum bcf_int8_vector_end = -127; /* INT8_MIN  + 1 */
enum bcf_int16_vector_end = -32767; /* INT16_MIN + 1 */
enum bcf_int32_vector_end = -2147483647; /* INT32_MIN + 1 */
enum bcf_int64_vector_end = -9223372036854775807LL; /* INT64_MIN + 1 */
enum bcf_str_vector_end = 0;
enum bcf_int8_missing = -128; /* INT8_MIN  */
enum bcf_int16_missing = -32767 - 1; /* INT16_MIN */
enum bcf_int32_missing = -2147483647 - 1; /* INT32_MIN */
enum bcf_int64_missing = -9223372036854775807LL - 1LL; /* INT64_MIN */
enum bcf_str_missing = 0x07;

// Limits on BCF values stored in given types.  Max values are the same
// as for the underlying type.  Min values are slightly different as
// the last 8 values for each type were reserved by BCFv2.2.
enum BCF_MAX_BT_INT8 = 0x7f; /* INT8_MAX  */
enum BCF_MAX_BT_INT16 = 0x7fff; /* INT16_MAX */
enum BCF_MAX_BT_INT32 = 0x7fffffff; /* INT32_MAX */
enum BCF_MIN_BT_INT8 = -120; /* INT8_MIN  + 8 */
enum BCF_MIN_BT_INT16 = -32760; /* INT16_MIN + 8 */
enum BCF_MIN_BT_INT32 = -2147483640; /* INT32_MIN + 8 */

extern __gshared uint bcf_float_vector_end;
extern __gshared uint bcf_float_missing;
void bcf_float_set(float* ptr, uint value);

extern (D) auto bcf_float_set_vector_end(T)(auto ref T x)
{
    return bcf_float_set(&x, bcf_float_vector_end);
}

extern (D) auto bcf_float_set_missing(T)(auto ref T x)
{
    return bcf_float_set(&x, bcf_float_missing);
}

int bcf_float_is_missing(float f);
int bcf_float_is_vector_end(float f);

int bcf_format_gt(bcf_fmt_t* fmt, int isample, kstring_t* str);

int bcf_enc_size(kstring_t* s, int size, int type);

int bcf_enc_inttype(c_long x);

int bcf_enc_int1(kstring_t* s, int x);

/// Return the value of a single typed integer.
/** @param      p    Pointer to input data block.
    @param      type One of the BCF_BT_INT* type codes
    @param[out] q    Location to store an updated value for p
    @return The integer value, or zero if @p type is not valid.

If @p type is not one of BCF_BT_INT8, BCF_BT_INT16, BCF_BT_INT32 or
BCF_BT_INT64, zero will be returned and @p *q will not be updated.
Otherwise, the integer value will be returned and @p *q will be set
to the memory location immediately following the integer value.

Cautious callers can detect invalid type codes by checking that *q has
actually been updated.
*/

// Invalid type.
long bcf_dec_int1(const(ubyte)* p, int type, ubyte** q);

/// Return the value of a single typed integer from a byte stream.
/** @param      p    Pointer to input data block.
    @param[out] q    Location to store an updated value for p
    @return The integer value, or zero if the type code was not valid.

Reads a one-byte type code from @p p, and uses it to decode an integer
value from the following bytes in @p p.

If the type is not one of BCF_BT_INT8, BCF_BT_INT16 or BCF_BT_INT32, zero
will be returned and @p *q will unchanged.  Otherwise, the integer value will
be returned and @p *q will be set to the memory location immediately following
the integer value.

Cautious callers can detect invalid type codes by checking that *q has
actually been updated.
*/
long bcf_dec_typed_int1(const(ubyte)* p, ubyte** q);

int bcf_dec_size(const(ubyte)* p, ubyte** q, int* type);

