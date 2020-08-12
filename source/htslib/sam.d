// htslib-1.9 sam.h as D module
// Changes include:
// Removed if(n)defs
// 	including BAM_NO_ID which is probably archaic
// Changed #defines to const/immutable
// Removed all HTS_RESULT_USED (__attribute__ ((__warn_unused_result__)))
// Do not #include "hts_defs.h"
// Change numeric #defines to enum int
// Changed string #define to 8-bit (char) string literal "xxx"c
// typedef struct to alias
// removed typedefs
// modified bitfields in struct and aligned(1)
// removed redundant struct declarations when declaring struct pointers
// ref is a reserved keyword in D; changed 'ref' to '_ref'
// Function prototypes taking fixed size array (e.g. ..., const char tag[2], ) should include ref in the D prototype 
// const char *p -> const(char) *p
// Replace localal definition with import kstring
/*  
Aliased function pointer typedefs:

    typedef int (*bam_plp_auto_f)(void *data, bam1_t *b);
    alias bam_plp_auto_f = int *function(void *data, bam1_t *b);

Updated parameter function pointers:
    int (*func)(void *data, const bam1_t *b, bam_pileup_cd *cd));
    int function(void *data, const bam1_t *b, bam_pileup_cd *cd) func);

Functions returning const must be rewritten as const(type)*func_name
*/
module htslib.sam;

import std.bitmanip;

extern (C):
/// @file htslib/sam.h
/// High-level SAM/BAM/CRAM sequence file operations.
/*
    Copyright (C) 2008, 2009, 2013-2017 Genome Research Ltd.
    Copyright (C) 2010, 2012, 2013 Broad Institute.

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

import core.stdc.stdint;
import htslib.hts;
import htslib.bgzf: BGZF;
import htslib.kstring: __kstring_t, kstring_t;

/// Highest SAM format version supported by this library
auto SAM_FORMAT_VERSION = "1.6"c;

/**********************
 *** SAM/BAM header ***
 **********************/

/**! @typedef
 @abstract Structure for the alignment header.
 @field n_targets   number of reference sequences
 @field l_text      length of the plain text in the header
 @field target_len  lengths of the reference sequences
 @field target_name names of the reference sequences
 @field text        plain text
 @field sdict       header dictionary
 */

struct bam_hdr_t { // @suppress(dscanner.style.phobos_naming_convention)
    int  n_targets;         /// number of targets (contigs)
    int ignore_sam_err; /// ? Ignore SAM format errors?
    uint l_text;        /// ? text data length
    uint *target_len;   /// ? number of targets (contigs)
    byte *cigar_tab;    /// ? CIGAR table
    char **target_name; /// ? C-style array of target (contig) names
    char *text;         /// ?
    void *sdict;        /// ?
}

/****************************
 *** CIGAR related macros ***
 ****************************/

enum {
    BAM_CMATCH      = 0,    /// CIGAR M, match
    BAM_CINS        = 1,    /// CIGAR I, insertion
    BAM_CDEL        = 2,    /// CIGAR D, deletion
    BAM_CREF_SKIP   = 3,    /// CIGAR N, refskip
    BAM_CSOFT_CLIP  = 4,    /// CIGAR S, soft clip
    BAM_CHARD_CLIP  = 5,    /// CIGAR H, hard clip
    BAM_CPAD        = 6,    /// CIGAR P, padding
    BAM_CEQUAL      = 7,    /// CIGAR =, equal
    BAM_CDIFF       = 8,    /// CIGAR X, differs
    BAM_CBACK       = 9     /// CIGAR B, back(?)
}

auto BAM_CIGAR_STR       = "MIDNSHP=XB"c;   /// internal
enum int BAM_CIGAR_SHIFT = 4;       /// internal
enum int BAM_CIGAR_MASK  = 0xf;     /// internal
enum int BAM_CIGAR_TYPE  = 0x3C1A7; /// internal magic const for bam_cigar_type()

pragma(inline, true)
{
/// CIGAR opcode for CIGAR uint (4 LSB)
auto bam_cigar_op(uint c)    { return ( c & BAM_CIGAR_MASK); }
/// CIGAR operation length for CIGAR uint (28 MSB >> 4)
auto bam_cigar_oplen(uint c) { return ( c >> BAM_CIGAR_SHIFT ); }
// Note that BAM_CIGAR_STR is padded to length 16 bytes below so that
// the array look-up will not fall off the end.  '?' is chosen as the
// padding character so it's easy to spot if one is emitted, and will
// result in a parsing failure (in sam_parse1(), at least) if read.
/// CIGAR opcode character for CIGAR uint
auto bam_cigar_opchr(uint c)    { return (BAM_CIGAR_STR ~ "??????"[bam_cigar_op(c)]); }
/// Generate CIGAR uint from length and operator
auto bam_cigar_gen(uint l, uint o) { return (l << BAM_CIGAR_SHIFT | o); }

/** bam_cigar_type returns a bit flag with:
 *   bit 1 set if the cigar operation consumes the query
 *   bit 2 set if the cigar operation consumes the reference
 *
 * For reference, the unobfuscated truth table for this function is:
 * BAM_CIGAR_TYPE  QUERY  REFERENCE
 * --------------------------------
 * BAM_CMATCH      1      1
 * BAM_CINS        1      0
 * BAM_CDEL        0      1
 * BAM_CREF_SKIP   0      1
 * BAM_CSOFT_CLIP  1      0
 * BAM_CHARD_CLIP  0      0
 * BAM_CPAD        0      0
 * BAM_CEQUAL      1      1
 * BAM_CDIFF       1      1
 * BAM_CBACK       0      0
 * --------------------------------
 */
auto bam_cigar_type(uint o) { return (BAM_CIGAR_TYPE >> (o << 1) & 3); }    // bit 1: consume query; bit 2: consume reference
} // end pragma(inline)

/**! @abstract the read is paired in sequencing, no matter whether it is mapped in a pair */
enum int BAM_FPAIRED       = 1;
/**! @abstract the read is mapped in a proper pair */
enum int BAM_FPROPER_PAIR  = 2;
/**! @abstract the read itself is unmapped; conflictive with BAM_FPROPER_PAIR */
enum int BAM_FUNMAP        = 4;
/**! @abstract the mate is unmapped */
enum int BAM_FMUNMAP       = 8;
/**! @abstract the read is mapped to the reverse strand */
enum int BAM_FREVERSE      =16;
/**! @abstract the mate is mapped to the reverse strand */
enum int BAM_FMREVERSE     =32;
/**! @abstract this is read1 */
enum int BAM_FREAD1        =64;
/**! @abstract this is read2 */
enum int BAM_FREAD2        =128;
/**! @abstract not primary alignment */
enum int BAM_FSECONDARY    =256;
/**! @abstract QC failure */
enum int BAM_FQCFAIL       =512;
/**! @abstract optical or PCR duplicate */
enum int BAM_FDUP          =1024;
/**! @abstract supplementary alignment */
enum int BAM_FSUPPLEMENTARY=2048;

/*************************
 *** Alignment records ***
 *************************/

/**! @typedef
 @abstract Structure for core alignment information.
 @field  tid     chromosome ID, defined by bam_hdr_t
 @field  pos     0-based leftmost coordinate
 @field  bin     bin calculated by bam_reg2bin()
 @field  qual    mapping quality
 @field  l_qname length of the query name
 @field  flag    bitwise flag
 @field  l_extranul length of extra NULs between qname & cigar (for alignment)
 @field  n_cigar number of CIGAR operations
 @field  l_qseq  length of the query sequence (read)
 @field  mtid    chromosome ID of next read in template, defined by bam_hdr_t
 @field  mpos    0-based leftmost coordinate of next read in template
 */
struct bam1_core_t { // @suppress(dscanner.style.phobos_naming_convention)
    int32_t  tid;       /// chromosome ID, defined by bam_hdr_t
    int32_t  pos;       /// 0-based leftmost coordinate
    uint16_t bin;       /// bin calculated by bam_reg2bin()
    uint8_t  qual;      /// mapping quality
    uint8_t  l_qname;   /// length of the query name
    uint16_t flag;      /// bitwise flag
    uint8_t  unused1;   /// padding
    uint8_t  l_extranul;/// length of extra NULs between qname & cigar (for alignment)
    uint32_t n_cigar;   /// number of CIGAR operations
    int32_t  l_qseq;    /// length of the query sequence (read)
    int32_t  mtid;      /// chromosome ID of next read in template, defined by bam_hdr_t
    int32_t  mpos;      /// 0-based leftmost coordinate of next read in template (mate)
    int32_t  isize;     /// ? template length
}

/**! @typedef
 @abstract Structure for one alignment.
 @field  core       core information about the alignment
 @field  l_data     current length of bam1_t::data
 @field  m_data     maximum length of bam1_t::data
 @field  data       all variable-length data, concatenated; structure: qname-cigar-seq-qual-aux

 @discussion Notes:

 1. qname is terminated by one to four NULs, so that the following
 cigar data is 32-bit aligned; core.l_qname includes these trailing NULs,
 while core.l_extranul counts the excess NULs (so 0 <= l_extranul <= 3).
 2. l_qseq is calculated from the total length of an alignment block
 on reading or from CIGAR.
 3. cigar data is encoded 4 bytes per CIGAR operation.
 4. seq is nybble-encoded according to bam_nt16_table.
 */
struct bam1_t { // @suppress(dscanner.style.phobos_naming_convention)
    bam1_core_t core;   /// core information about the alignment
    int         l_data; /// current length of bam1_t::data
    uint32_t    m_data; /// maximum length of bam1_t::data
    uint8_t     *data;  /// all variable-length data, concatenated; structure: qname-cigar-seq-qual-aux
    uint64_t    id;     /// ???
}

pragma(inline, true) {
/**! @function
 @abstract  Get whether the query is on the reverse strand
 @param  b  pointer to an alignment
 @return    boolean true if query is on the reverse strand
 */
bool bam_is_rev(bam1_t *b) { return ( ((*b).core.flag & BAM_FREVERSE) != 0 ); }
/**! @function
 @abstract  Get whether the query's mate is on the reverse strand
 @param  b  pointer to an alignment
 @return    boolean true if query's mate on the reverse strand
 */
bool bam_is_mrev(bam1_t *b) { return( ((*b).core.flag & BAM_FMREVERSE) != 0); }
/**! @function
 @abstract  Get the name of the query
 @param  b  pointer to an alignment
 @return    pointer to the name string, null terminated
 */
auto bam_get_qname(bam1_t *b) { return (cast(char*)(*b).data); }
/**! @function
 @abstract  Get the CIGAR array
 @param  b  pointer to an alignment
 @return    pointer to the CIGAR array

 @discussion In the CIGAR array, each element is a 32-bit integer. The
 lower 4 bits gives a CIGAR operation and the higher 28 bits keep the
 length of a CIGAR.
 */
auto bam_get_cigar(bam1_t *b) { return cast(uint *)((*b).data + (*b).core.l_qname); }
/**! @function
 @abstract  Get query sequence
 @param  b  pointer to an alignment
 @return    pointer to sequence

 @discussion Each base is encoded in 4 bits: 1 for A, 2 for C, 4 for G,
 8 for T and 15 for N. Two bases are packed in one byte with the base
 at the higher 4 bits having smaller coordinate on the read. It is
 recommended to use bam_seqi() macro to get the base.
 */
auto bam_get_seq(bam1_t *b) { return ((*b).data + ((*b).core.n_cigar<<2) + (*b).core.l_qname); }
/**! @function
 @abstract  Get query quality
 @param  b  pointer to an alignment
 @return    pointer to quality string
 */
////#define bam_get_qual(b)  ((b)->data + ((b)->core.n_cigar<<2) + (b)->core.l_qname + (((b)->core.l_qseq + 1)>>1))
auto bam_get_qual(bam1_t *b) { return (*b).data + ((*b).core.n_cigar<<2) + 
                                    (*b).core.l_qname + 
                                        (((*b).core.l_qseq + 1)>>1); }

/**! @function
 @abstract  Get auxiliary data
 @param  b  pointer to an alignment
 @return    pointer to the concatenated auxiliary data
 */
////#define bam_get_aux(b)   ((b)->data + ((b)->core.n_cigar<<2) + (b)->core.l_qname + (((b)->core.l_qseq + 1)>>1) + (b)->core.l_qseq)
auto bam_get_aux(bam1_t *b) { return ((*b).data + ((*b).core.n_cigar<<2) + 
                                    (*b).core.l_qname + 
                                        (((*b).core.l_qseq + 1)>>1) + 
                                            (*b).core.l_qseq); }
/**! @function
 @abstract  Get length of auxiliary data
 @param  b  pointer to an alignment
 @return    length of the concatenated auxiliary data
 */
////#define bam_get_l_aux(b) ((b)->l_data - ((b)->core.n_cigar<<2) - (b)->core.l_qname - (b)->core.l_qseq - (((b)->core.l_qseq + 1)>>1))
auto bam_get_l_aux(bam1_t *b) { return ((*b).l_data - 
                                ((*b).core.n_cigar<<2) - 
                                    (*b).core.l_qname - 
                                        (*b).core.l_qseq - 
                                            (((*b).core.l_qseq + 1)>>1)); }
/**! @function
 @abstract  Get a base on read
 @param  s  Query sequence returned by bam_get_seq()
 @param  i  The i-th position, 0-based
 @return    4-bit integer representing the base.
 */
////#define bam_seqi(s, i) ((s)[(i)>>1] >> ((~(i)&1)<<2) & 0xf)
auto bam_seqi(ubyte *s, uint i) { return ((s)[(i)>>1] >> ((~(i)&1)<<2) & 0xf);   }
//auto bam_seqi(char *s, uint i) { return ((s)[(i)>>1] >> ((~(i)&1)<<2) & 0xf); }
} // end pragma(inline)

/**************************
 *** Exported functions ***
 **************************/

    /***************
     *** BAM I/O ***
     ***************/

    /// init a BAM header
    bam_hdr_t *bam_hdr_init();
    /// read BAM header from fp
    bam_hdr_t *bam_hdr_read(BGZF *fp);
    /// write BAM header to fp
    int bam_hdr_write(BGZF *fp, const(bam_hdr_t) *h);
    /// destroy a BAM header
    void bam_hdr_destroy(bam_hdr_t *h);
    /// Get contig tid 
    int bam_name2id(bam_hdr_t *h, const(char) *_ref);
    /// duplicate a BAM header
    bam_hdr_t* bam_hdr_dup(const(bam_hdr_t) *h0);

    /// init a BAM record
    bam1_t *bam_init1();
    /// destroy a BAM record
    void bam_destroy1(bam1_t *b);
    /// read BAM record from fp
    int bam_read1(BGZF *fp, bam1_t *b);
    /// write BAM record to fp
    int bam_write1(BGZF *fp, const(bam1_t) *b);
    /// copy a BAM record from dst to src
    bam1_t *bam_copy1(bam1_t *bdst, const(bam1_t) *bsrc);
    /// duplicate a BAM record
    bam1_t *bam_dup1(const(bam1_t) *bsrc);

    /// How much query is consumed by CIGAR op(s)?
    int bam_cigar2qlen(int n_cigar, const(uint32_t) *cigar);
    /// How much reference is consumed by CIGAR op(s)?
    int bam_cigar2rlen(int n_cigar, const(uint32_t) *cigar);

    /**!
      @abstract Calculate the rightmost base position of an alignment on the
      reference genome.

      @param  b  pointer to an alignment
      @return    the coordinate of the first base after the alignment, 0-based

      @discussion For a mapped read, this is just b->core.pos + bam_cigar2rlen.
      For an unmapped read (either according to its flags or if it has no cigar
      string), we return b->core.pos + 1 by convention.
    */
    int32_t bam_endpos(const(bam1_t) *b);

    /// returns negative value on error
    int   bam_str2flag(const(char) *str);    /** returns negative value on error */
    /// The string must be freed by the user
    char *bam_flag2str(int flag);   /** The string must be freed by the user */

    /*************************
     *** BAM/CRAM indexing ***
     *************************/

    // These BAM iterator functions work only on BAM files.  To work with either
    // BAM or CRAM files use the sam_index_load() & sam_itr_*() functions.
    deprecated("These BAM iterator functions work only on BAM files. "
        ~ "To work with either BAM or CRAM files use the sam_index_load() & sam_itr_*() functions.")
    {
    ////#define bam_itr_destroy(iter) hts_itr_destroy(iter)
    alias bam_itr_destroy = hts_itr_destroy;
    ////#define bam_itr_queryi(idx, tid, beg, end) sam_itr_queryi(idx, tid, beg, end)
    alias bam_itr_queryi = sam_itr_queryi;
    ////#define bam_itr_querys(idx, hdr, region) sam_itr_querys(idx, hdr, region)
    alias bam_itr_querys = sam_itr_querys;
    ////#define bam_itr_next(htsfp, itr, r) hts_itr_next((htsfp)->fp.bgzf, (itr), (r), 0)
    pragma(inline, true)
    auto bam_itr_next(htsFile *htsfp, hts_itr_t *itr, void *r) { return hts_itr_next(htsfp.fp.bgzf, itr, r, null); }
    }

// Load/build .csi or .bai BAM index file.  Does not work with CRAM.
// It is recommended to use the sam_index_* functions below instead.
deprecated("Load/build .csi or .bai BAM index file.  Does not work with CRAM. "
    ~ "It is recommended to use the sam_index_* functions below instead.")
{
////#define bam_index_load(fn) hts_idx_load((fn), HTS_FMT_BAI)
pragma(inline, true) auto bam_index_load(const(char) *fn) { return hts_idx_load(fn, HTS_FMT_BAI); }
////#define bam_index_build(fn, min_shift) (sam_index_build((fn), (min_shift)))
alias bam_index_build = sam_index_build;
}

/// Load a BAM (.csi or .bai) or CRAM (.crai) index file
/** @param fp  File handle of the data file whose index is being opened
    @param fn  BAM/CRAM/etc filename to search alongside for the index file
    @return  The index, or NULL if an error occurred.
*/
hts_idx_t *sam_index_load(htsFile *fp, const(char) *fn);

/// Load a specific BAM (.csi or .bai) or CRAM (.crai) index file
/** @param fp     File handle of the data file whose index is being opened
    @param fn     BAM/CRAM/etc data file filename
    @param fnidx  Index filename, or NULL to search alongside @a fn
    @return  The index, or NULL if an error occurred.
*/
hts_idx_t *sam_index_load2(htsFile *fp, const(char) *fn, const(char) *fnidx);

/// Generate and save an index file
/** @param fn        Input BAM/etc filename, to which .csi/etc will be added
    @param min_shift Positive to generate CSI, or 0 to generate BAI
    @return  0 if successful, or negative if an error occurred (usually -1; or
             -2: opening fn failed; -3: format not indexable; -4:
             failed to create and/or save the index)
*/
int sam_index_build(const(char) *fn, int min_shift);

/// Generate and save an index to a specific file
/** @param fn        Input BAM/CRAM/etc filename
    @param fnidx     Output filename, or NULL to add .bai/.csi/etc to @a fn
    @param min_shift Positive to generate CSI, or 0 to generate BAI
    @return  0 if successful, or negative if an error occurred (see
             sam_index_build for error codes)
*/
int sam_index_build2(const(char) *fn, const(char) *fnidx, int min_shift);
/// ditto
int sam_index_build3(const(char) *fn, const(char) *fnidx, int min_shift, int nthreads);

    ////#define sam_itr_destroy(iter) hts_itr_destroy(iter)
    alias sam_itr_destroy = hts_itr_destroy;
    /// SAM/BAM/CRAM iterator query by integer tid/start/end
    hts_itr_t *sam_itr_queryi(const(hts_idx_t) *idx, int tid, int beg, int end);
    /// SAM/BAM/CRAM iterator query by string ("chr:start-end")
    hts_itr_t *sam_itr_querys(const(hts_idx_t) *idx, bam_hdr_t *hdr, const(char) *region);
    /// SAM/BAM/CRAM iterator query by region list
    hts_itr_multi_t *sam_itr_regions(const(hts_idx_t) *idx, bam_hdr_t *hdr, hts_reglist_t *reglist, uint regcount);

    ////#define sam_itr_next(htsfp, itr, r) hts_itr_next((htsfp)->fp.bgzf, (itr), (r), (htsfp))
    pragma(inline, true)
    auto sam_itr_next(htsFile *htsfp, hts_itr_t *itr, void *r) { return hts_itr_next(htsfp.fp.bgzf, itr, r, htsfp); }
    ////#define sam_itr_multi_next(htsfp, itr, r) hts_itr_multi_next((htsfp), (itr), (r))
    //auto sam_itr_multi_next(htsFile *htsfp, hts_itr_t *itr, void *r) { return hts_itr_multi_next(htsfp, itr, r); }
    alias sam_itr_multi_next = hts_itr_multi_next;

    /***************
     *** SAM I/O ***
     ***************/

    //#define sam_open(fn, mode) (hts_open((fn), (mode)))
    alias sam_open = hts_open;
    //#define sam_open_format(fn, mode, fmt) (hts_open_format((fn), (mode), (fmt)))
    alias sam_open_format = hts_open_format;
    //#define sam_close(fp) hts_close(fp)
    alias sam_close = hts_close;

    /// open SAM/BAM/CRAM
    int sam_open_mode(char *mode, const(char) *fn, const(char) *format);

    /// A version of sam_open_mode that can handle ,key=value options.
    /// The format string is allocated and returned, to be freed by the caller.
    /// Prefix should be "r" or "w",
    char *sam_open_mode_opts(const(char) *fn,
                             const(char) *mode,
                             const(char) *format);

    alias samFile = htsFile;
    /// parse SAM/BAM/CRAM header from text
    bam_hdr_t *sam_hdr_parse(int l_text, const(char) *text);
    /// read SAM/BAM/CRAM header frpm fp
    bam_hdr_t *sam_hdr_read(samFile *fp);
    /// write SAM/BAM/CRAM header to fp
    int sam_hdr_write(samFile *fp, const(bam_hdr_t) *h);
    /// edit SAM/BAM/CRAM header @HD line by key/value (see SAM specs section 1.3)
    int sam_hdr_change_HD(bam_hdr_t *h, const(char) *key, const(char) *val);

    /// parse text SAM line
    int sam_parse1(kstring_t *s, bam_hdr_t *h, bam1_t *b);
    /// emit text SAM line
    int sam_format1(const(bam_hdr_t) *h, const(bam1_t) *b, kstring_t *str);

    /**!
     *  @return >= 0 on successfully reading a new record, -1 on end of stream, < -1 on error
     **/
    int sam_read1(samFile *fp, bam_hdr_t *h, bam1_t *b);
    /// Return >= 0 on successfully writing a new record; <= -1 on error
    int sam_write1(samFile *fp, const(bam_hdr_t) *h, const(bam1_t) *b);

    /*************************************
     *** Manipulating auxiliary fields ***
     *************************************/

/// Return a pointer to an aux record
/** @param b   Pointer to the bam record
    @param tag Desired aux tag
    @return Pointer to the tag data, or NULL if tag is not present or on error
    If the tag is not present, this function returns NULL and sets errno to
    ENOENT.  If the bam record's aux data is corrupt (either a tag has an
    invalid type, or the last record is incomplete) then errno is set to
    EINVAL and NULL is returned.
 */
uint8_t *bam_aux_get(const(bam1_t) *b, const ref char[2] tag);

/// Get an integer aux value
/** @param s Pointer to the tag data, as returned by bam_aux_get()
    @return The value, or 0 if the tag was not an integer type
    If the tag is not an integer type, errno is set to EINVAL.  This function
    will not return the value of floating-point tags.
*/
int64_t bam_aux2i(const(uint8_t) *s);

/// Get an integer aux value
/** @param s Pointer to the tag data, as returned by bam_aux_get()
    @return The value, or 0 if the tag was not an integer type
    If the tag is not an numeric type, errno is set to EINVAL.  The value of
    integer flags will be returned cast to a double.
*/
double bam_aux2f(const(uint8_t) *s);

/// Get a character aux value
/** @param s Pointer to the tag data, as returned by bam_aux_get().
    @return The value, or 0 if the tag was not a character ('A') type
    If the tag is not a character type, errno is set to EINVAL.
*/
char bam_aux2A(const(uint8_t) *s);

/// Get a string aux value
/** @param s Pointer to the tag data, as returned by bam_aux_get().
    @return Pointer to the string, or NULL if the tag was not a string type
    If the tag is not a string type ('Z' or 'H'), errno is set to EINVAL.
*/
char *bam_aux2Z(const(uint8_t) *s);

/// Get the length of an array-type ('B') tag
/** @param s Pointer to the tag data, as returned by bam_aux_get().
    @return The length of the array, or 0 if the tag is not an array type.
    If the tag is not an array type, errno is set to EINVAL.
 */
uint32_t bam_auxB_len(const(uint8_t) *s);

/// Get an integer value from an array-type tag
/** @param s   Pointer to the tag data, as returned by bam_aux_get().
    @param idx 0-based Index into the array
    @return The idx'th value, or 0 on error.
    If the array is not an integer type, errno is set to EINVAL.  If idx
    is greater than or equal to  the value returned by bam_auxB_len(s),
    errno is set to ERANGE.  In both cases, 0 will be returned.
 */
int64_t bam_auxB2i(const(uint8_t) *s, uint32_t idx);

/// Get a floating-point value from an array-type tag
/** @param s   Pointer to the tag data, as returned by bam_aux_get().
    @param idx 0-based Index into the array
    @return The idx'th value, or 0.0 on error.
    If the array is not a numeric type, errno is set to EINVAL.  This can
    only actually happen if the input record has an invalid type field.  If
    idx is greater than or equal to  the value returned by bam_auxB_len(s),
    errno is set to ERANGE.  In both cases, 0.0 will be returned.
 */
double bam_auxB2f(const(uint8_t) *s, uint32_t idx);

/// Append tag data to a bam record
/* @param b    The bam record to append to.
   @param tag  Tag identifier
   @param type Tag data type
   @param len  Length of the data in bytes
   @param data The data to append
   @return 0 on success; -1 on failure.
If there is not enough space to store the additional tag, errno is set to
ENOMEM.  If the type is invalid, errno may be set to EINVAL.  errno is
also set to EINVAL if the bam record's aux data is corrupt.
*/
int bam_aux_append(bam1_t *b, const ref char[2] tag, char type, int len, const(uint8_t) *data);

/// Delete tag data from a bam record
/* @param b The bam record to update
   @param s Pointer to the tag to delete, as returned by bam_aux_get().
   @return 0 on success; -1 on failure
   If the bam record's aux data is corrupt, errno is set to EINVAL and this
   function returns -1;
*/
int bam_aux_del(bam1_t *b, uint8_t *s);

/// Update or add a string-type tag
/* @param b    The bam record to update
   @param tag  Tag identifier
   @param len  The length of the new string
   @param data The new string
   @return 0 on success, -1 on failure
   This function will not change the ordering of tags in the bam record.
   New tags will be appended to any existing aux records.

   On failure, errno may be set to one of the following values:

   EINVAL: The bam record's aux data is corrupt or an existing tag with the
   given ID is not of type 'Z'.

   ENOMEM: The bam data needs to be expanded and either the attempt to
   reallocate the data buffer failed or the resulting buffer would be
   longer than the maximum size allowed in a bam record (2Gbytes).
*/
int bam_aux_update_str(bam1_t *b, const ref char[2] tag, int len, const(char) *data);

/// Update or add an integer tag
/* @param b    The bam record to update
   @param tag  Tag identifier
   @param val  The new value
   @return 0 on success, -1 on failure
   This function will not change the ordering of tags in the bam record.
   New tags will be appended to any existing aux records.

   On failure, errno may be set to one of the following values:

   EINVAL: The bam record's aux data is corrupt or an existing tag with the
   given ID is not of an integer type (c, C, s, S, i or I).

   EOVERFLOW (or ERANGE on systems that do not have EOVERFLOW): val is
   outside the range that can be stored in an integer bam tag (-2147483647
   to 4294967295).

   ENOMEM: The bam data needs to be expanded and either the attempt to
   reallocate the data buffer failed or the resulting buffer would be
   longer than the maximum size allowed in a bam record (2Gbytes).
*/
int bam_aux_update_int(bam1_t *b, const ref char[2] tag, int64_t val);

/// Update or add a floating-point tag
/* @param b    The bam record to update
   @param tag  Tag identifier
   @param val  The new value
   @return 0 on success, -1 on failure
   This function will not change the ordering of tags in the bam record.
   New tags will be appended to any existing aux records.

   On failure, errno may be set to one of the following values:

   EINVAL: The bam record's aux data is corrupt or an existing tag with the
   given ID is not of a float type.

   ENOMEM: The bam data needs to be expanded and either the attempt to
   reallocate the data buffer failed or the resulting buffer would be
   longer than the maximum size allowed in a bam record (2Gbytes).
*/
int bam_aux_update_float(bam1_t *b, const ref char[2] tag, float val);

/// Update or add an array tag
/* @param b     The bam record to update
   @param tag   Tag identifier
   @param type  Data type (one of c, C, s, S, i, I or f)
   @param items Number of items
   @param data  Pointer to data
   @return 0 on success, -1 on failure
   The type parameter indicates the how the data is interpreted:

   Letter code | Data type | Item Size (bytes)
   ----------- | --------- | -----------------
   c           | int8_t    | 1
   C           | uint8_t   | 1
   s           | int16_t   | 2
   S           | uint16_t  | 2
   i           | int32_t   | 4
   I           | uint32_t  | 4
   f           | float     | 4

   This function will not change the ordering of tags in the bam record.
   New tags will be appended to any existing aux records.  The bam record
   will grow or shrink in order to accomodate the new data.

   The data parameter must not point to any data in the bam record itself or
   undefined behaviour may result.

   On failure, errno may be set to one of the following values:

   EINVAL: The bam record's aux data is corrupt, an existing tag with the
   given ID is not of an array type or the type parameter is not one of
   the values listed above.

   ENOMEM: The bam data needs to be expanded and either the attempt to
   reallocate the data buffer failed or the resulting buffer would be
   longer than the maximum size allowed in a bam record (2Gbytes).
*/
int bam_aux_update_array(bam1_t *b, const ref char[2] tag,
                         uint8_t type, uint32_t items, void *data);

/**************************
 *** Pileup and Mpileup ***
 **************************/

/*! @typedef
 @abstract Generic pileup 'client data'.

 @discussion The pileup iterator allows setting a constructor and
 destructor function, which will be called every time a sequence is
 fetched and discarded.  This permits caching of per-sequence data in
 a tidy manner during the pileup process.  This union is the cached
 data to be manipulated by the "client" (the caller of pileup).
*/
union bam_pileup_cd {
    void    *p; /// data ptr
    int64_t i;  /// ?
    double  f;  /// ?
}

/**! @typedef
 @abstract Structure for one alignment covering the pileup position.
 @field  b          pointer to the alignment
 @field  qpos       position of the read base at the pileup site, 0-based
 @field  indel      indel length; 0 for no indel, positive for ins and negative for del
 @field  level      the level of the read in the "viewer" mode
 @field  is_del     1 iff the base on the padded read is a deletion
 @field  is_head    ???
 @field  is_tail    ???
 @field  is_refskip ???
 @field  aux        ???

 @discussion See also bam_plbuf_push() and bam_lplbuf_push(). The
 difference between the two functions is that the former does not
 set bam_pileup1_t::level, while the later does. Level helps the
 implementation of alignment viewers, but calculating this has some
 overhead.
 */
struct bam_pileup1_t {
    bam1_t  *b;     /// bam record
    int32_t qpos;   /// query position
    /// ???
    int indel, level;
    pragma(msg, "bam_pileup1_t: bitfield order assumed starting with LSB");
    //uint32_t is_del:1, is_head:1, is_tail:1, is_refskip:1, aux:28;
    mixin(bitfields!(
        bool, "is_del",  1,
        bool, "is_head", 1,
        bool, "is_tail", 1,
        bool, "is_refskip", 1,
        uint, "aux", 28 ));
    bam_pileup_cd cd; /// generic per-struct data, owned by caller.
}

///typedef int (*bam_plp_auto_f)(void *data, bam1_t *b);
alias bam_plp_auto_f = int *function(void *data, bam1_t *b);

/// opaque pileup data defined in sam.c
struct __bam_plp_t;
///typedef struct __bam_plp_t *bam_plp_t;
alias bam_plp_t = __bam_plp_t*;

/// opaque mpileup data defined in sam.c
struct __bam_mplp_t;
///typedef struct __bam_mplp_t *bam_mplp_t;
alias bam_mplp_t = __bam_mplp_t*;

    /**
     *  bam_plp_init() - sets an iterator over multiple
     *  @func:      see mplp_func in bam_plcmd.c in samtools for an example. Expected return
     *              status: 0 on success, -1 on end, < -1 on non-recoverable errors
     *  @data:      user data to pass to @func
     */
    /// NB: maxcnt records default is 8000
    bam_plp_t bam_plp_init(bam_plp_auto_f func, void *data);
    /// destroy pileup iterator
    void bam_plp_destroy(bam_plp_t iter);
    /// add bam record to pileup iterator
    int bam_plp_push(bam_plp_t iter, const(bam1_t) *b);
    /// Prepares next pileup position in bam records collected by bam_plp_auto -> user func -> bam_plp_push. Returns
    /// pointer to the piled records if next position is ready or NULL if there is not enough records in the
    /// buffer yet (the current position is still the maximum position across all buffered reads).
    const(bam_pileup1_t)*bam_plp_next(bam_plp_t iter, int *_tid, int *_pos, int *_n_plp);
    /// ???
    const(bam_pileup1_t)*bam_plp_auto(bam_plp_t iter, int *_tid, int *_pos, int *_n_plp);
    /// set maximum pileup records returned (init default is 8000)
    void bam_plp_set_maxcnt(bam_plp_t iter, int maxcnt);
    /// reset pileup
    void bam_plp_reset(bam_plp_t iter);

    /**
     *  bam_plp_constructor() - sets a callback to initialise any per-pileup1_t fields.
     *  @plp:       The bam_plp_t initialised using bam_plp_init.
     *  @func:      The callback function itself.  When called, it is given the
     *              data argument (specified in bam_plp_init), the bam structure and
     *              a pointer to a locally allocated bam_pileup_cd union.  This union
     *              will also be present in each bam_pileup1_t created.
     */
    void bam_plp_constructor(bam_plp_t plp,
                             int function(void *data, const(bam1_t) *b, bam_pileup_cd *cd) func);
    /// per-pileup1_t field destructor (may be needed if used bam_plp_constructor())
    void bam_plp_destructor(bam_plp_t plp,
                            int function(void *data, const(bam1_t) *b, bam_pileup_cd *cd) func);

    /// initialize new mpileup iterator
    bam_mplp_t bam_mplp_init(int n, bam_plp_auto_f func, void **data);
    /**
     *  bam_mplp_init_overlaps() - if called, mpileup will detect overlapping
     *  read pairs and for each base pair set the base quality of the
     *  lower-quality base to zero, thus effectively discarding it from
     *  calling. If the two bases are identical, the quality of the other base
     *  is increased to the sum of their qualities (capped at 200), otherwise
     *  it is multiplied by 0.8.
     */
    void bam_mplp_init_overlaps(bam_mplp_t iter);
    /// destroy mpileup iterator
    void bam_mplp_destroy(bam_mplp_t iter);
    /// set maximum mpileup records returned per subpileup (init default of each pileup in the mpileup is 8000)
    void bam_mplp_set_maxcnt(bam_mplp_t iter, int maxcnt);
    /// ???
    int bam_mplp_auto(bam_mplp_t iter, int *_tid, int *_pos, int *n_plp, const bam_pileup1_t **plp);
    /// reset mpileup
    void bam_mplp_reset(bam_mplp_t iter);
    /// see bam_plp_constructor
    void bam_mplp_constructor(bam_mplp_t iter,
                              int function(void *data, const(bam1_t) *b, bam_pileup_cd *cd) func);
    /// see bam_plp_destructor
    void bam_mplp_destructor(bam_mplp_t iter,
                             int function(void *data, const(bam1_t) *b, bam_pileup_cd *cd) func);


/***********************************
 * BAQ calculation and realignment *
 ***********************************/

int sam_cap_mapq(bam1_t *b, const(char) *_ref, int ref_len, int thres);

/// Calculate BAQ scores
/** @param b   BAM record
    @param ref     Reference sequence
    @param ref_len Reference sequence length
    @param flag    Flags, see description
    @return 0 on success \n
           -1 if the read was unmapped, zero length, had no quality values, did not have at least one M, X or = CIGAR operator, or included a reference skip. \n
           -3 if BAQ alignment has already been done and does not need to be applied, or has already been applied. \n
           -4 if alignment failed (most likely due to running out of memory)

This function calculates base alignment quality (BAQ) values using the method
described in "Improving SNP discovery by base alignment quality", Heng Li,
Bioinformatics, Volume 27, Issue 8 (https://doi.org/10.1093/bioinformatics/btr076).

The following @param flag bits can be used:

Bit 0: Adjust the quality values using the BAQ values

 If set, the data in the BQ:Z tag is used to adjust the quality values, and
 the BQ:Z tag is renamed to ZQ:Z.

 If clear, and a ZQ:Z tag is present, the quality values are reverted using
 the data in the tag, and the tag is renamed to BQ:Z.

Bit 1: Use "extended" BAQ.

 Changes the BAQ calculation to increase sensitivity at the expense of
 reduced specificity.

Bit 2: Recalculate BAQ, even if a BQ tag is present.

 Force BAQ to be recalculated.  Note that a ZQ:Z tag will always disable
 recalculation.

@bug
If the input read has both BQ:Z and ZQ:Z tags, the ZQ:Z one will be removed.
Depending on what previous processing happened, this may or may not be the
correct thing to do.  It would be wise to avoid this situation if possible.
*/

int sam_prob_realn(bam1_t *b, const(char) *_ref, int ref_len, int flag);