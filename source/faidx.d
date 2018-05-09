module dhtslib.faidx;

import std.string;
import core.stdc.stdlib : malloc, free;

import dhtslib.htslib.faidx;

/// Build index for a FASTA or bgzip-compressed FASTA file.
/**  @param  fn  FASTA file name
     @param  fnfai Name of .fai file to build.
     @param  fngzi Name of .gzi file to build (if fn is bgzip-compressed).
     @return     0 on success; or -1 on failure

If fnfai is NULL, ".fai" will be appended to fn to make the FAI file name.
If fngzi is NULL, ".gzi" will be appended to fn for the GZI file.  The GZI
file will only be built if fn is bgzip-compressed.
*/
bool buildFastaIndex(string fn, string fnfai = "", string fngzi = "")
{
    // int fai_build3(const char *fn, const char *fnfai, const char *fngzi);
    int ret = fai_build3( toStringz(fn),
                            fnfai ? toStringz(fnfai) : null,
                            fngzi ? toStringz(fngzi) : null);
    
    if (ret == 0) return true;
    if (ret == -1) return false;
}

struct IndexedFastaFile {

    faidx_t *faix;

    this(string fn, bool create=false)
    {
        if (create) {
            this.faidx = fai_load3( toStringz(fn), null, null, FAI_CREATE);
        }
        else {
            this.faidx = fai_load( toStringz(fn) , null, null, FAI_CREATE);
        }
    }
    ~this()
    {
        fai_destroy(this.faidx);
    }

    /** fetchSequence (overloaded)
     *
     *  Region in the format "chr2:20,000-30,000"
     */
    string fetchSequence(string region)
    {
        char *fetchedSeq;
        int fetchedLen;
        // char *fai_fetch(const faidx_t *fai, const char *reg, int *len);
        // @param  len  Length of the region; -2 if seq not present, -1 general error
        fetchedSeq = fai_fetch(this.faidx, toStringz(region), &fetchedLen);

        string seq = fromStringz(fetchedSeq);
        free(fetchedSeq);

        if (fetchedLen == -1) throw new Exception("fai_fetch: unknown error");
        else if (fetchedLen == -2) throw new Exception("fai_fetch: sequence not found");

        return seq;
    }

    /** fetchSequence (overloaded)
     *
     * htslib API for my reference:
     * @param  fai  Pointer to the faidx_t struct
     * @param  c_name Region name
     * @param  p_beg_i  Beginning position number (zero-based)
     * @param  p_end_i  End position number (zero-based)
     * @param  len  Length of the region; -2 if c_name not present, -1 general error
     * @return      Pointer to the sequence; null on failure
     *  The returned sequence is allocated by `malloc()` family and should be destroyed
     *  by end users by calling `free()` on it.
     *
     */
    string fetchSequence(string contig, int start, int end)
    {
        // char *faidx_fetch_seq(const faidx_t *fai, const char *c_name, int p_beg_i, int p_end_i, int *len);
        char *fetchedSeq;
        int fetchedLen;
        fetchedSeq = fai_fetch(this.faidx, toStringz(contig), start, end, &fetchedLen);
        
        string seq = fromStringz(fetchedSeq);
        free(fetchedSeq);

        if (fetchedLen == -1) throw new Exception("fai_fetch: unknown error");
        else if (fetchedLen == -2) throw new Exception("fai_fetch: sequence not found");

        return seq;
    }

    bool hasSeq(string seqname)
    {
        // int faidx_has_seq(const faidx_t *fai, const char *seq);
        return cast(bool) faidx_has_seq(this.faidx, toStringz(seqname) );
    }

    @property nSeq()
    {
        return faidx_nseq(this.faidx);
    }

    /// Return the name of the i'th sequence
    string seqName(int i)
    {
        // const(char) *faidx_iseq(const faidx_t *fai, int i);
        // TODO determine if this property is zero indexed or one indexed
        if (i > this.nSeq) throw new Exception("seqName: sequece number not present");

        return fromStringz( faidx_iseq(this.faidx, i) );
    }

    /// Return sequence length, -1 if not present
    int seqLen(string seqname)
    {
        // TODO should I check for -1 and throw exception or pass to caller?
        // int faidx_seq_len(const faidx_t *fai, const char *seq);
        int l = faidx_seq_len(this.faidx, toStringz(seqname) );
        if ( l == -1 ) throw new Exception("seqLen: sequence name not found");
        return l;
    }
}