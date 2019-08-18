module dhtslib.faidx;

import std.string;
import core.stdc.stdlib : malloc, free;

import dhtslib.htslib.faidx;

/** Coodinate Systems, since htslib sequence fetching methods surprisingly use zero-based, closed */
enum CoordSystem
{
    zbho = 0,   /// zero-based, half-open
    zbc,        /// zero-based, closed (behavior of faidx_fetch_seq)
    obho,       /// one-based, half-open
    obc         /// one-based, closed
}

/** Build index for a FASTA or bgzip-compressed FASTA file.

    @param  fn  FASTA file name
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
    const int ret = fai_build3( toStringz(fn),
                            fnfai ? toStringz(fnfai) : null,
                            fngzi ? toStringz(fngzi) : null);
    
    if (ret == 0) return true;
    if (ret == -1) return false;
    assert(0);
}

/** FASTA file with .fai or .gzi index

Reads existing FASTA file, optionally creating FASTA index if one does not exist.

Convenient member fns to get no. of sequences, get sequence names and lengths,
test for membership, and rapidly fetch sequence at offset.
*/
struct IndexedFastaFile {

    private faidx_t *faidx;
    private int refct;      // Postblit refcounting in case the object is passed around

    this(this)
    {
        refct++;
    }
    /// construct from filename, optionally creating index if it does not exist
    /// throws Exception (TODO: remove) if file DNE, or if index DNE unless create->true
    this(string fn, bool create=false)
    {
        if (create) {
            this.faidx = fai_load3( toStringz(fn), null, null, fai_load_options.FAI_CREATE);
            if (this.faidx is null) throw new Exception("Unable to load or create the FASTA index.");
        }
        else {
            this.faidx = fai_load3( toStringz(fn) , null, null, 0);
            if (this.faidx is null) throw new Exception("Unable to load the FASTA index.");
        }
        refct = 1;
    }
    ~this()
    {
        if (--refct == 0)
            fai_destroy(this.faidx);
    }

    /// Enable BGZF cacheing (size: bytes)
    void setCacheSize(int size)
    {
        import dhtslib.htslib.bgzf : BGZF, bgzf_set_cache_size;
        bgzf_set_cache_size(this.faidx.bgzf, size);
    }

    /// Enable multithreaded decompression of this FA/FQ
    /// Reading fn body of bgzf_mt, this actually ADDS threads (rather than setting)
    /// but we'll retain name for consistency with setCacheSize
    /// NB: IN A REAL-WORLD TEST (swiftover) CALLING setThreads(1) doubled runtime(???)
    void setThreads(int nthreads)
    {
        import dhtslib.htslib.bgzf : BGZF, bgzf_mt;
        // third parameter, n_sub_blks is not used in htslib 1.9; .h suggests value 64-256
        bgzf_mt(this.faidx.bgzf, nthreads, 64);
    }

    /// Fetch sequence in region by assoc array-style lookup:
    /// `string sequence = fafile["chr2:20123-30456"]`
    auto opIndex(string region)
    {
        return fetchSequence(region);
    }

    /// Fetch sequence in region by assoc array-style lookup:
    /// `string sequence = fafile.fetchSequence("chr2:20123-30456")`
    string fetchSequence(string region)
    {
        char *fetchedSeq;
        int fetchedLen;
        // char *fai_fetch(const faidx_t *fai, const char *reg, int *len);
        // @param  len  Length of the region; -2 if seq not present, -1 general error
        fetchedSeq = fai_fetch(this.faidx, toStringz(region), &fetchedLen);

        string seq = fromStringz(fetchedSeq).idup;
        free(fetchedSeq);

        if (fetchedLen == -1) throw new Exception("fai_fetch: unknown error");
        else if (fetchedLen == -2) throw new Exception("fai_fetch: sequence not found");

        return seq;
    }

    /// Fetch sequence in region by multidimensional slicing:
    /// `string sequence = fafile["chr2", 20123 .. 30456]`
    ///
    /// Sadly, $ to represent max length is not supported
    auto opIndex(string contig, int[2] pos)
    {
        return fetchSequence(contig, pos[0], pos[1]);
    }
    /// ditto
    int[2] opSlice(size_t dim)(int start, int end) if (dim == 1)
    in { assert(start >= 0); }
    do
    {
        return [start, end];
    }

    /// Fetch sequence in region by multidimensional slicing:
    /// `string sequence = fafile.fetchSequence("chr2", 20123, 30456)`
    string fetchSequence(CoordSystem cs = CoordSystem.zbho)(string contig, int start, int end)
    {
        char *fetchedSeq;
        int fetchedLen;

        static if (cs == CoordSystem.zbho) {
            end--;
        }
        else static if (cs == CoordSystem.zbc) {
            // ok
        }
        else static if (cs == CoordSystem.obho) {
            start--;
            end--;
        }
        else static if (cs == CoordSystem.obc) {
            start--;
        }
        /* htslib API for my reference:
         *
         * char *faidx_fetch_seq(const faidx_t *fai, const char *c_name, int p_beg_i, int p_end_i, int *len);
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
        fetchedSeq = faidx_fetch_seq(this.faidx, toStringz(contig), start, end, &fetchedLen);
        
        string seq = fromStringz(fetchedSeq).idup;
        free(fetchedSeq);

        if (fetchedLen == -1) throw new Exception("fai_fetch: unknown error");
        else if (fetchedLen == -2) throw new Exception("fai_fetch: sequence not found");

        assert(seq.length == fetchedLen);
        return seq;
    }

    /// Test whether the FASTA file/index contains string seqname
    bool hasSeq(string seqname)
    {
        // int faidx_has_seq(const faidx_t *fai, const char *seq);
        return cast(bool) faidx_has_seq(this.faidx, toStringz(seqname) );
    }

    /// Return the number of sequences in the FASTA file/index
    @property auto nSeq()
    {
        return faidx_nseq(this.faidx);
    }

    /// Return the name of the i'th sequence
    string seqName(int i)
    {
        // const(char) *faidx_iseq(const faidx_t *fai, int i);
        // TODO determine if this property is zero indexed or one indexed
        if (i > this.nSeq) throw new Exception("seqName: sequece number not present");

        return fromStringz( faidx_iseq(this.faidx, i) ).idup;
    }

    /// Return sequence length, -1 if not present
    int seqLen(string seqname)
    {
        // TODO should I check for -1 and throw exception or pass to caller?
        // int faidx_seq_len(const faidx_t *fai, const char *seq);
        const int l = faidx_seq_len(this.faidx, toStringz(seqname) );
        if ( l == -1 ) throw new Exception("seqLen: sequence name not found");
        return l;
    }
}
