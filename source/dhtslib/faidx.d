/**
This module provides a wrapper, `IndexedFastaFile` over a FASTA.
If an index exists, it will be used for rapid random access.
If an index does not exist, one will be built.

The wrapper provides the ability to list sequence names (i.e., chromosomes/contigs)
in the FASTA, efficiently retrieve sequences (by contig, start, end)

Sequence caching and multithreaded BGZF decompression are supported.
*/

module dhtslib.faidx;

import std.string;
import core.stdc.stdlib : malloc, free;

import dhtslib.coordinates;
import htslib.faidx;

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

    invariant(){
        assert(refct >= 0);
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
        fai_set_cache_size(this.faidx, size);
    }

    /// Enable multithreaded decompression of this FA/FQ
    /// Reading fn body of bgzf_mt, this actually ADDS threads (rather than setting)
    /// but we'll retain name for consistency with setCacheSize
    /// NB: IN A REAL-WORLD TEST (swiftover) CALLING setThreads(1) doubled runtime(???)
    deprecated("disabled until faidx_t again made non-opaque")
    void setThreads(int nthreads)
    {
        import htslib.bgzf : BGZF, bgzf_mt;
        // third parameter, n_sub_blks is not used in htslib 1.9; .h suggests value 64-256
        //bgzf_mt(this.faidx.bgzf, nthreads, 64);
    }

    /// Fetch sequence in region by assoc array-style lookup:
    /// `string sequence = fafile["chr2:20123-30456"]`
    auto opIndex(CoordSystem cs)(ChromCoordinates!cs region)
    {
        return fetchSequence(region);
    }

    /// Fetch sequence in region by multidimensional slicing:
    /// `string sequence = fafile["chr2", 20123 .. 30456]`
    ///
    /// Sadly, $ to represent max length is not supported
    auto opIndex(CoordSystem cs)(string contig, Coordinates!cs coords)
    {
        return fetchSequence(contig, coords);
    }
    
    /// ditto
    auto opSlice(size_t dim, Basis bs)(Coordinate!bs start, Coordinate!bs end) if (dim == 1)
    in { assert(start >= 0); assert(start <= end); }
    do
    {
        auto coords = Coordinates!(getCoordinateSystem!(bs, End.open))(start, end);
        return coords;
    }

    /// Fetch sequence in region by assoc array-style lookup:
    /// `string sequence = fafile.fetchSequence("chr2:20123-30456")`

    string fetchSequence(CoordSystem cs)(ChromCoordinates!cs region)
    {
        return fetchSequence(region.chrom, region.coords);
    }

    /// Fetch sequencing in a region by function call with contig, start, end
    /// `string sequence = fafile.fetchSequence("chr2", 20123, 30456)`
    string fetchSequence(CoordSystem cs)(string contig, Coordinates!cs coords)
    {
        char *fetchedSeq;
        long fetchedLen;

        // Convert given coordinates in system cs to zero-based, closed (faidx_fetch_seq)
        auto newcoords = coords.to!(CoordSystem.zbc);
        /* htslib API for my reference:
         *
         * char *faidx_fetch_seq64(const faidx_t *fai, const char *c_name, hts_pos_t p_beg_i, hts_pos_t p_end_i, hts_pos_t *len);
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
        fetchedSeq = faidx_fetch_seq64(this.faidx, toStringz(contig), newcoords.start, newcoords.end, &fetchedLen);
        
        if (fetchedLen == -1) throw new Exception("fai_fetch: unknown error");
        else if (fetchedLen == -2) throw new Exception("fai_fetch: sequence not found");

        string seq = fromStringz(fetchedSeq).idup;
        free(fetchedSeq);

        assert(seq.length == fetchedLen);
        return seq;
    }

    /// Test whether the FASTA file/index contains string seqname
    bool hasSeq(const(char)[] seqname)
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
    /// NOTE: there is no 64 bit equivalent of this function (yet) in htslib-1.10
    int seqLen(const(char)[] seqname)
    {
        // TODO should I check for -1 and throw exception or pass to caller?
        // int faidx_seq_len(const faidx_t *fai, const char *seq);
        const int l = faidx_seq_len(this.faidx, toStringz(seqname) );
        if ( l == -1 ) throw new Exception("seqLen: sequence name not found");
        return l;
    }
}

debug(dhtslib_unittest) unittest
{
    import dhtslib.sam;
    import htslib.hts_log;
    import std.path : buildPath, dirName;

    hts_set_log_level(htsLogLevel.HTS_LOG_INFO);
    hts_log_info(__FUNCTION__, "Testing IndexedFastaFile");
    hts_log_info(__FUNCTION__, "Loading test file");
    auto fai = IndexedFastaFile(buildPath(dirName(dirName(dirName(__FILE__))),"htslib","test","ce.fa"));

    fai.setCacheSize(4000000);

    assert(fai.fetchSequence("CHROMOSOME_I",ZBHO(0, 29)) == "GCCTAAGCCTAAGCCTAAGCCTAAGCCTA");
    assert(fai.fetchSequence("CHROMOSOME_I",OBHO(1, 30)) == "GCCTAAGCCTAAGCCTAAGCCTAAGCCTA");
    assert(fai.fetchSequence("CHROMOSOME_I",ZBC(0, 28)) == "GCCTAAGCCTAAGCCTAAGCCTAAGCCTA");
    assert(fai.fetchSequence("CHROMOSOME_I",OBC(1, 29)) == "GCCTAAGCCTAAGCCTAAGCCTAAGCCTA");

    assert(fai.fetchSequence("CHROMOSOME_I",ZBHO(1, 29)) == "CCTAAGCCTAAGCCTAAGCCTAAGCCTA");
    assert(fai.fetchSequence("CHROMOSOME_I",OBHO(2, 30)) == "CCTAAGCCTAAGCCTAAGCCTAAGCCTA");
    assert(fai.fetchSequence("CHROMOSOME_I",ZBC(1, 28)) == "CCTAAGCCTAAGCCTAAGCCTAAGCCTA");
    assert(fai.fetchSequence("CHROMOSOME_I",OBC(2, 29)) == "CCTAAGCCTAAGCCTAAGCCTAAGCCTA");

    assert(fai["CHROMOSOME_I",ZB(0) .. ZB(29)] == "GCCTAAGCCTAAGCCTAAGCCTAAGCCTA");
    import std.stdio;
    writeln(fai["CHROMOSOME_I",OB(1) .. OB(30)]);
    assert(fai["CHROMOSOME_I",OB(1) .. OB(30)] == "GCCTAAGCCTAAGCCTAAGCCTAAGCCTA");
    
    assert(fai.fetchSequence(ChromZBHO("CHROMOSOME_I:0-29")) == "GCCTAAGCCTAAGCCTAAGCCTAAGCCTA");
    assert(fai.fetchSequence(ChromOBHO("CHROMOSOME_I:1-30")) == "GCCTAAGCCTAAGCCTAAGCCTAAGCCTA");
    assert(fai.fetchSequence(ChromZBC("CHROMOSOME_I:0-28")) == "GCCTAAGCCTAAGCCTAAGCCTAAGCCTA");
    assert(fai.fetchSequence(ChromOBC("CHROMOSOME_I:1-29")) == "GCCTAAGCCTAAGCCTAAGCCTAAGCCTA");

    assert(fai.seqLen("CHROMOSOME_I") == 1009800);
    assert(fai.nSeq == 7);

    assert(fai.hasSeq("CHROMOSOME_I"));
    assert(!fai.hasSeq("CHROMOSOME_"));

    assert(fai.seqName(0) == "CHROMOSOME_I");
}
