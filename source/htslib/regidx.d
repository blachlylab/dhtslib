/** htslib-1.9 regidx.h as D module
 *
 *  Changes include:
 *$(LIST
 *      Removed if(n)defs
 *      Changed ^typedef struct {...} <name>$ to ^struct <name> {...}$
 *      Commented out the REGITR_{START,END,PAYLOAD,OVERLAP} #define macros
 *      aliased typedef'd function pointers
 *)
 */
module htslib.regidx;

import std.stdint : uint32_t;

extern (C):
/// @file htslib/regidx.h
/// Region indexing.
/*
    Copyright (C) 2014 Genome Research Ltd.

    Author: Petr Danecek <pd3@sanger.ac.uk>

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
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
    THE SOFTWARE.
*/

/*
    Regions indexing with an optional payload. Inspired by samtools/bedidx.c.
    This code is intended as future replacement of bcf_sr_regions_t.

    Example of usage:

        // Init the parser and print regions. In this example the payload is a
        // pointer to a string. For the description of parse_custom and
        // free_custom functions, see regidx_parse_f and regidx_free_f below,
        // and for working example see test/test-regidx.c.
        regidx_t *idx = regidx_init(in_fname,parse_custom,free_custom,sizeof(char*),NULL);

        // Query overlap with chr:from-to
        regitr_t itr;
        if ( regidx_overlap(idx, chr,from,to, &itr) ) printf("There is an overlap!\n");

        while ( REGITR_OVERLAP(itr,from,to) )
        {
            printf("[%d,%d] overlaps with [%d,%d], payload=%s\n", from,to,
                REGITR_START(itr), REGITR_END(itr), REGITR_PAYLOAD(itr,char*));
            itr.i++;
        }

        regidx_destroy(regs);
*/

/// struct defined in regidx.c, we will leave as opaque
extern struct _regidx_t;    // @suppress(dscanner.style.phobos_naming_convention)
alias regidx_t = _regidx_t;
/// region (start, end)
struct reg_t // @suppress(dscanner.style.phobos_naming_convention)
{
    /// (start, end) -- unclear if 0-based/open
    uint32_t start, end;
}
/// region iterator
struct regitr_t // @suppress(dscanner.style.phobos_naming_convention)
{
    /// ???
    int i, n;
    reg_t *reg;     /// ???
    void *payload;  /// ???
}

/+
#define REGITR_START(itr) (itr).reg[(itr).i].start
#define REGITR_END(itr)   (itr).reg[(itr).i].end
#define REGITR_PAYLOAD(itr,type_t) ((type_t*)(itr).payload)[(itr).i]
#define REGITR_OVERLAP(itr,from,to) (itr.i < itr.n && REGITR_START(itr)<=to && REGITR_END(itr)>=from )
+/
pragma(inline, true)
{
    /// Get the start or end coordinate of the region iterator's current region
    uint32_t regitr_start(regitr_t *itr){ return itr.reg[itr.i].start; }
    /// ditto
    uint32_t regitr_end(regitr_t *itr)  { return itr.reg[itr.i].end;   }

    /// Get the payload of the region iterator's current region
    auto regitr_payload(T)(regitr_t *itr)   { return cast(T*)itr.payload[itr.i]; }  // looks super unsafe

    /// Does the given (from, to) overlap the region?    
    auto regitr_overlap(T)(regitr_t *itr, T from, T to)
    {
        return itr.i < itr.n && regitr_start(itr) <= to && regitr_end(itr) >= from;
    }
}

/*
 *  regidx_parse_f - Function to parse one input line, such as regidx_parse_bed
 *  or regidx_parse_tab below. The function is expected to set `chr_from` and
 *  `chr_to` to point to first and last character of chromosome name and set
 *  coordinates `reg->start` and `reg->end` (0-based, inclusive). If
 *  regidx_init() was called with non-zero payload_size, the `payload` points
 *  to a memory location of the payload_size and `usr` is data passed to
 *  regidx_init(). Any memory allocated by the function will be freed by
 *  regidx_free_f on regidx_destroy().
 *
 *  Return value: 0 on success, -1 to skip a record, -2 on fatal error.
 */
//typedef int  (*regidx_parse_f)(const char *line, char **chr_beg, char **chr_end, reg_t *reg, void *payload, void *usr);
alias regidx_parse_f =
    int function(const(char) *line, char **chr_beg, char **chr_end, reg_t *reg, void *payload, void *usr);
//typedef void (*regidx_free_f)(void *payload);
alias regidx_free_f = void function(void *payload);

/// regidx_parse_f for BED, CHROM,FROM,TO (0-based,right-open)
int regidx_parse_bed(const(char)*, char**, char**, reg_t*, void*, void*);   // CHROM,FROM,TO (0-based,right-open)
/// regidx_Parse_f for "TAB", CHROM,POS (1-based, inclusive)
int regidx_parse_tab(const(char)*, char**, char**, reg_t*, void*, void*);   // CHROM,POS (1-based, inclusive)

/**
 *  regidx_init() - creates new index
 *  @param fname:  input file name or NULL if regions will be added one-by-one via regidx_insert()
 *  @param parsef: regidx_parse_bed, regidx_parse_tab or see description of regidx_parse_f. If NULL,
 *                 the format will be autodected, currently either regidx_parse_tab (the default) or
 *                 regidx_parse_bed (file must be named 'bed' or 'bed.gz') will be used. Note that
 *                 the exact autodetection algorithm will change.
 *  @param freef:  NULL or see description of regidx_parse_f
 *  @param payload_size: 0 with regidx_parse_bed, regidx_parse_tab or see regidx_parse_f
 *  @param usr:    optional user data passed to regidx_parse_f
 *
 *  Returns index on success or NULL on error.
 */
regidx_t *regidx_init(const(char) *fname, regidx_parse_f parsef, regidx_free_f freef, size_t payload_size, void *usr);

/**
 *  regidx_destroy() - free memory allocated by regidx_init
 */
void regidx_destroy(regidx_t *idx);

/**
 *  regidx_overlap() - check overlap of the location chr:from-to with regions
 *  @param start,end:   0-based start, end coordinate (inclusive)
 *  @param itr:         pointer to iterator, can be NULL if not needed
 *
 *  Returns 0 if there is no overlap or 1 if overlap is found. The overlapping
 *  regions can be iterated as shown in the example above.
 */
int regidx_overlap(regidx_t *idx, const(char) *chr, uint32_t start, uint32_t end, regitr_t *itr);

/**
 *  regidx_insert() - add a new region.
 *
 *  After last region has been added, call regidx_insert(idx,NULL) to
 *  build the index.
 *
 *  Returns 0 on success or -1 on error.
 */
int regidx_insert(regidx_t *idx, char *line);

/**
 *  regidx_seq_names() - return list of all sequence names
 */
char **regidx_seq_names(regidx_t *idx, int *n);

/**
 *  regidx_seq_nregs() - number of regions
 *  regidx_nregs()  - total number of regions
 */
int regidx_seq_nregs(regidx_t *idx, const(char) *seq);
/// ditto
int regidx_nregs(regidx_t *idx);
