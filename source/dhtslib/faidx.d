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
import std.typecons : Tuple;
import core.stdc.stdlib : malloc, free;

import dhtslib.coordinates;
import dhtslib.memory;
import dhtslib.util;
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

    private Faidx faidx;
    
    /// construct from filename, optionally creating index if it does not exist
    /// throws Exception (TODO: remove) if file DNE, or if index DNE unless create->true
    this(string fn, bool create=false)
    {
        if (create) {
            this.faidx = Faidx(fai_load3( toStringz(fn), null, null, fai_load_options.FAI_CREATE));
            if (this.faidx is null) throw new Exception("Unable to load or create the FASTA index.");
        }
        else {
            this.faidx = Faidx(fai_load3( toStringz(fn) , null, null, 0));
            if (this.faidx is null) throw new Exception("Unable to load the FASTA index.");
        }
    }

    /// Explicit postblit to avoid 
    /// https://github.com/blachlylab/dhtslib/issues/122
    this(this)
    {
        this.faidx = faidx;
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

    private struct OffsetType
    {
        ptrdiff_t offset;
        alias offset this;

        // supports e.g. $ - x
        OffsetType opBinary(string s, T)(T val)
        {
            mixin("return OffsetType(offset " ~ s ~ " val);");
        }

        invariant
        {
            assert(this.offset <= 0, "Offset from end should be zero or negative");
        }
    }
    /** Array-end `$` indexing hack courtesy of Steve Schveighoffer
        https://forum.dlang.org/post/rl7a56$nad$1@digitalmars.com

        Requires in addition to opDollar returning a bespoke non-integral type
        a series of overloads for opIndex and opSlice taking this type
    */
    OffsetType opDollar(size_t dim)() if(dim == 1)
    {
        return OffsetType.init;
    }
    /// opSlice as Coordinate and an offset
    /// i.e [ZB(2) .. $]
    auto opSlice(size_t dim, Basis bs)(Coordinate!bs start, OffsetType off) if(dim == 1)
    {
        return Tuple!(Coordinate!bs, OffsetType)(start, off);
    }
    /// opSlice as two offset
    /// i.e [$-2 .. $]
    auto opSlice(size_t dim)(OffsetType start, OffsetType end) if(dim == 1)
    {
        return Tuple!(OffsetType, OffsetType)(start, end);
    }
    /// opIndex coordinate and Offset
    /// i.e fai["chrom1", ZB(1) .. $]
    auto opIndex(Basis bs)(string ctg, Tuple!(Coordinate!bs, OffsetType) coords)
    {
        auto end = this.seqLen(ctg) + coords[1];
        auto endCoord = ZB(end);
        auto newEndCoord = endCoord.to!bs;
        auto c = Interval!(getCoordinateSystem!(bs,End.open))(coords[0], newEndCoord);
        return fetchSequence(ctg, c);
    }

    /// opIndex two Offsets
    /// i.e fai["chrom1", $-2 .. $]
    auto opIndex(string ctg, Tuple!(OffsetType, OffsetType) coords)
    {
        auto start = this.seqLen(ctg) + coords[0];
        auto end = this.seqLen(ctg) + coords[1];
        auto c = ZBHO(start, end);
        return fetchSequence(ctg, c);
    }
    /// opIndex one offset
    /// i.e fai["chrom1", $-1]
    auto opIndex(string ctg, OffsetType endoff)
    {
        auto end = this.seqLen(ctg) + endoff.offset;
        auto coords = Interval!(CoordSystem.zbho)(end, end + 1);
        return fetchSequence(ctg, coords);
    }

    /// Fetch sequence in region by assoc array-style lookup:
    /// Uses htslib string region parsing
    /// `string sequence = fafile["chr2:20123-30456"]`
    auto opIndex(string region)
    {
        auto coords = getIntervalFromString(region);
        return fetchSequence(coords.contig, coords.interval);
    }

    /// Fetch sequence in region by multidimensional slicing:
    /// `string sequence = fafile["chr2", 20123 .. 30456]`
    ///
    /// Sadly, $ to represent max length is not supported
    auto opIndex(CoordSystem cs)(string contig, Interval!cs coords)
    {
        return fetchSequence(contig, coords);
    }
    
    /// ditto
    auto opSlice(size_t dim, Basis bs)(Coordinate!bs start, Coordinate!bs end) if (dim == 1)
    in { assert(start >= 0); assert(start <= end); }
    do
    {
        auto coords = Interval!(getCoordinateSystem!(bs, End.open))(start, end);
        return coords;
    }

    /// Fetch sequencing in a region by function call with contig, start, end
    /// `string sequence = fafile.fetchSequence("chr2", 20123, 30456)`
    string fetchSequence(CoordSystem cs)(string contig, Interval!cs coords)
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

    assert(fai.seqLen("CHROMOSOME_I") == 1009800);
    assert(fai.nSeq == 7);

    assert(fai.hasSeq("CHROMOSOME_I"));
    assert(!fai.hasSeq("CHROMOSOME_"));

    assert(fai.seqName(0) == "CHROMOSOME_I");
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

    assert(fai.fetchSequence("CHROMOSOME_II",ZBHO(0, 29)) == "CCTAAGCCTAAGCCTAAGCCTAAGCCTAA");
    assert(fai.fetchSequence("CHROMOSOME_II",OBHO(1, 30)) == "CCTAAGCCTAAGCCTAAGCCTAAGCCTAA");
    assert(fai.fetchSequence("CHROMOSOME_II",ZBC(0, 28)) == "CCTAAGCCTAAGCCTAAGCCTAAGCCTAA");
    assert(fai.fetchSequence("CHROMOSOME_II",OBC(1, 29)) == "CCTAAGCCTAAGCCTAAGCCTAAGCCTAA");
    assert(fai["CHROMOSOME_II:1-29"] == "CCTAAGCCTAAGCCTAAGCCTAAGCCTAA");
    // "GTCAACACAGACCGTTAATTTTGGGAAGTTGAGAAATTCGCTAGTTTCTG"
    import std.stdio;
    assert(fai["CHROMOSOME_II",$-50 .. $] == "GTCAACACAGACCGTTAATTTTGGGAAGTTGAGAAATTCGCTAGTTTCTG");
    
    assert(fai["CHROMOSOME_II",OB(4951) .. $] == "GTCAACACAGACCGTTAATTTTGGGAAGTTGAGAAATTCGCTAGTTTCTG");
    assert(fai["CHROMOSOME_II",ZB(4950) .. $] == "GTCAACACAGACCGTTAATTTTGGGAAGTTGAGAAATTCGCTAGTTTCTG");
    assert(fai["CHROMOSOME_II",$-1] == "G");
}
