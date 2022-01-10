module dhtslib.sam.reader;

import core.stdc.stdio : SEEK_SET;
import std.format;
import std.stdio : writeln, writefln, stderr, File;
import std.string : fromStringz, toStringz;
import std.typecons : Tuple;
import std.range.interfaces : InputRange, inputRangeObject;
import std.range.primitives : isInputRange, ElementType;
import std.algorithm : map;
import std.array : array;
import std.utf : toUTFz;

import htslib.hts;
import htslib.hfile : hdopen, hclose, hFILE, off_t, htell, hseek;
import htslib.bgzf : bgzf_tell, bgzf_seek;

import htslib.hts_log;
import htslib.kstring;
import htslib.sam;
import dhtslib.sam.record;
import dhtslib.coordinates;
import dhtslib.sam.header;
import dhtslib.sam : cmpInterval, cmpRegList;
import dhtslib.memory;
import dhtslib.util;
import dhtslib.file;

alias SAMFile = SAMReader;
/**
Encapsulates a SAM/BAM file.

Implements InputRange interface using htslib calls.
If indexed, Random-access query via multidimensional slicing.
*/
struct SAMReader
{
    /// HtslibFile
    private HtslibFile f;

    /// header struct
    SAMHeader header;

    /** Create a representation of SAM/BAM/CRAM file from given filename or File

    Params:
        fn =            string filename (complete path passed to htslib; may support S3:// https:// etc.)
        extra_threads = extra threads for compression/decompression
                        -1 use all available cores (default)
                        0  use no extra threads
                        >1 add indicated number of threads (to a default of 1)
    */
    this(T)(T f, int extra_threads = -1)
    if (is(T == string) || is(T == File))
    {
        import std.parallelism : totalCPUs;
        
        // open file
        this.f = HtslibFile(f);

        if (extra_threads == -1)
        {
            if ( totalCPUs > 1)
            {
                hts_log_info(__FUNCTION__,
                        format("%d CPU cores detected; enabling multithreading", totalCPUs));
                // hts_set_threads adds N _EXTRA_ threads, so totalCPUs - 1 seemed reasonable,
                // but overcomitting by 1 thread (i.e., passing totalCPUs) buys an extra 3% on my 2-core 2013 Mac
                this.f.setExtraThreads(totalCPUs);
            }
        }
        else if (extra_threads > 0)
        {
            if ((extra_threads + 1) > totalCPUs)
                hts_log_warning(__FUNCTION__, "More threads requested than CPU cores detected");
            this.f.setExtraThreads(extra_threads);
        }
        else if (extra_threads == 0)
        {
            hts_log_debug(__FUNCTION__, "Zero extra threads requested");
        }
        else
        {
            hts_log_warning(__FUNCTION__, "Invalid negative number of extra threads requested");
        }

        // read header
        this.f.loadHeader;
        this.header = SAMHeader(this.f.bamHdr);
    }

    /// number of reference sequences; from bam_hdr_t
    deprecated("use these properties from SAMHeader")
    @property int nTargets() const
    {
        return this.header.nTargets;
    }
    alias n_targets = nTargets;

    /// length of specific reference sequence by number (`tid`)
    deprecated("use these properties from SAMHeader")
    uint targetLength(int target) const
    in (target >= 0)
    in (target < this.nTargets)
    {
        return this.header.targetLength(target);
    }
    alias targetLen = targetLength;
    alias target_len = targetLength;

    /// lengths of the reference sequences
    deprecated("use these properties from SAMHeader")
    @property uint[] targetLengths() const
    {
        return this.header.targetLengths;
    }
    alias targetLens = targetLengths;
    alias target_lens = targetLengths;

    /// names of the reference sequences
    deprecated("use these properties from SAMHeader")
    @property string[] targetNames() const
    {
        return this.header.targetNames;
    }
    alias target_names = targetNames;

    /// reference contig name to integer id
    deprecated("use these properties from SAMHeader")
    int targetId(string name)
    {
        return this.header.targetId(name);
    }
    alias target_id = targetId;


    /// fetch is provided as a PySAM compatible synonym for `query`
    alias fetch = query;

    /** Query a region and return matching alignments as InputRange
    *
    *   Query on (chr, start, end) may take several forms:
    *
    *   1. `query(region)` with a string-based "region" form (e.g. chr1:1000-2000)
            - Variant: pass an array of query region strings: `query([reg1, reg, ...])`
    *   2. `query(chr, start, end)` with a combination of function parameters for
    *       contig, start, and end (where contig may be either a string or the numeric
    *       `tid` from BAM header; it would be uncommon to use this directly)
    *
    *   NOTE THAT THERE IS AN OFF-BY-ONE DIFFERENCE IN THE TWO METHODS ABOVE!
    *   Region string based coordinates assume the first base of the reference
    *   is 1 (e.g., chrX:1-100 yields the first 100 bases), whereas with the
    *   integer function parameter versions, the coordinates are zero-based, half-open
    *   (e.g., <chrX, 0, 100> yields the first 100 bases).
    *
    *   We also support array indexing on object of type SAMReader directly
    *   in one of two above styles:
    *       1. `bamfile[region-string]`
    *       2. `bamfile[contig, start .. end]` with contig like no. 2 above
    *
    *   The D convention `$` operator marking length of array is supported.
    *
    *   Finally, the region string is parsed by underlying htslib's `hts_parse_region`
    *   and has special semantics available:
    *
    *   region          | Outputs
    *   --------------- | -------------
    *   REF             | All reads with RNAME REF
    *   REF:            | All reads with RNAME REF
    *   REF:START       | Reads with RNAME REF overlapping START to end of REF
    *   REF:-END        | Reads with RNAME REF overlapping start of REF to END
    *   REF:START-END   | Reads with RNAME REF overlapping START to END
    *   .               | All reads from the start of the file
    *   *               | Unmapped reads at the end of the file (RNAME '*' in SAM)
    *
    *
    *   Examples:
    *   ```
    *   bamfile = SAMReader("whatever.bam");
    *   auto reads1 = bamfile.query("chr1:1-500");
    *   auto reads2 = bamfile.query("chr2", 0, 500);
    *   auto reads3 = bamfile["chr3", 0 .. 500];
    *
    *   auto reads4 = bamfile["chrX", $-500 .. $];  // last 500 nt
    *
    *   auto reads5 = bamfile.query("chrY");    // entirety of chrY
    *
    *   // When colon present in reference name (e.g. HLA additions in GRCh38)
    *   // wrap the ref name in { } (this is an htslib convention; see hts_parse_region)
    *   auto reads6 = bamfile.query("{HLA-DRB1*12:17}:1-100");
    *   ```
    */ 
    auto query(CoordSystem cs)(string chrom, Interval!cs coords)
    in (!this.header.isNull)
    {
        auto tid = this.header.targetId(chrom);
        return query(tid, coords);
    }

    /// ditto
    auto query(CoordSystem cs)(int tid, Interval!cs coords)
    in (!this.header.isNull)
    {
        /// convert to zero-based half-open
        auto newcoords = coords.to!(CoordSystem.zbho);

        /// load index
        if(this.f.idx == null)
            this.f.loadHtsIndex;
        if (this.f.idx == null)
        {
            if(this.f.tbx == null)
                this.f.loadTabixIndex;
            if(this.f.tbx == null){
                hts_log_error(__FUNCTION__, "BAI/TABIX index not found");
                throw new Exception("Cannot query");
            }
        }
        auto newF =  this.f.dup;
        newF.resetToFirstRecord;
        return newF.query!Bam1(tid, newcoords.start, newcoords.end)
                .map!(x=>SAMRecord(bam_dup1(x), header));
    }


    /// ditto
    auto query(string[] regions)
    in (!this.header.isNull)
    {
        /// load index
        if(this.f.idx == null)
            this.f.loadHtsIndex;
        if (this.f.idx == null)
        {
            hts_log_error(__FUNCTION__, "BAI index not found");
            throw new Exception("Cannot multi region query");
        }
        return this.f.query(regions).map!(x=>SAMRecord(x, header));
    }

    /// opIndex with a list of string regions
    /// bam[["chr1:1-3","chr2:4-50"]]
    auto opIndex(string[] regions)
    {
        return query(regions);
    }

    /// opIndex with a string contig and an Interval
    /// bam["chr1", ZBHO(1,3)]
    auto opIndex(CoordSystem cs)(string contig, Interval!cs coords)
    {
        return query(contig, coords);
    }

    /// opIndex with an int tid and an Interval
    /// bam[0, ZBHO(1,3)]
    auto opIndex(CoordSystem cs)(int tid, Interval!cs coords)
    {
        return query(tid, coords);
    }

    /// opIndex with a string contig and a Coordinate
    /// bam["chr1", ZB(1)]
    auto opIndex(Basis bs)(string contig, Coordinate!bs pos)
    {
        auto coords = Interval!(getCoordinateSystem!(bs,End.open))(pos, pos + 1);
        return query(contig, coords);
    }

    /// opIndex with an int tid and a Coordinate
    /// bam[0, ZB(1)]
    auto opIndex(Basis bs)(int tid, Coordinate!bs pos)
    {
        auto coords = Interval!(getCoordinateSystem!(bs,End.open))(pos, pos + 1);
        return query(tid, coords);
    }

    /// opIndex with a string contig and two Coordinates
    /// bam["chr1", ZB(1), ZB(3)]
    deprecated("use multidimensional slicing with second parameter as range ([\"chr1\", 1 .. 2])")
    auto opIndex(Basis bs)(string tid, Coordinate!bs pos1, Coordinate!bs pos2)
    {
        auto coords = Interval!(getCoordinateSystem!(bs,End.open))(pos1, pos2);
        return query(tid, coords);
    }

    /// opIndex with an int tid and two Coordinates
    /// bam[0, ZB(1), ZB(3)]
    deprecated("use multidimensional slicing with second parameter as range ([20, 1 .. 2])")
    auto opIndex(Basis bs)(int tid, Coordinate!bs pos1, Coordinate!bs pos2)
    {
        auto coords = Interval!(getCoordinateSystem!(bs,End.open))(pos1, pos2);
        return query(tid, coords);
    }

    /// opSlice with two Coordinates
    /// [ZB(1) .. ZB(3)]
    auto opSlice(size_t dim, Basis bs)(Coordinate!bs start, Coordinate!bs end) if(dim == 1)
    {
        assert(end > start);
        return Interval!(getCoordinateSystem!(bs, End.open))(start, end);
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
    /// opIndex with a string contig and an Offset
    /// bam["chr1",$-2]
    auto opIndex(string ctg, OffsetType endoff)
    {
        auto tid = this.header.targetId(ctg);
        auto end = this.header.targetLength(tid) + endoff.offset;
        // TODO review: is targetLength the last nt, or targetLength - 1 the last nt?
        auto coords = Interval!(CoordSystem.zbho)(end, end + 1);
        return query(tid, coords);
    }
    /// opIndex with an int tid and an Offset
    /// bam[0,$-2]
    auto opIndex(int tid, OffsetType endoff)
    {
        auto end = this.header.targetLength(tid) + endoff.offset;
        // TODO review: is targetLength the last nt, or targetLength - 1 the last nt?
        auto coords = Interval!(CoordSystem.zbho)(end, end + 1);
        return query(tid, coords);
    }

    /// opSlice as Coordinate and an offset
    /// i.e [ZB(2) .. $]
    auto opSlice(size_t dim, Basis bs)(Coordinate!bs start, OffsetType off) if(dim == 1)
    {
        return Tuple!(Coordinate!bs, OffsetType)(start, off);
    }

    /// opIndex with a string contig and a Coordinate and Offset
    /// bam["chr1",ZB(1) .. $]
    auto opIndex(Basis bs)(string ctg, Tuple!(Coordinate!bs, OffsetType) coords)
    {
        auto tid = this.header.targetId(ctg);
        auto end = this.header.targetLength(tid) + coords[1];
        auto endCoord = ZB(end);
        auto newEndCoord = endCoord.to!bs;
        auto c = Interval!(getCoordinateSystem!(bs,End.open))(coords[0], newEndCoord);
        return query(tid, c);
    }

    /// opIndex with an int tid and a Coordinate and Offset
    /// bam[0,ZB(1) .. $]
    auto opIndex(Basis bs)(int tid, Tuple!(Coordinate!bs, OffsetType) coords)
    {
        auto end = this.header.targetLength(tid) + coords[1];
        auto endCoord = ZB(end);
        auto newEndCoord = endCoord.to!bs;
        auto c = Interval!(getCoordinateSystem!(bs,End.open))(coords[0], newEndCoord);
        return query(tid, c);
    }

    /// opSlice as two offset
    /// i.e [$-2 .. $]
    auto opSlice(size_t dim)(OffsetType start, OffsetType end) if(dim == 1)
    {
        return Tuple!(OffsetType, OffsetType)(start, end);
    }

    /// opIndex two Offsets
    /// i.e fai["chrom1", $-2 .. $]
    auto opIndex(string ctg, Tuple!(OffsetType, OffsetType) coords)
    {
        auto tid = this.header.targetId(ctg);
        auto start = this.header.targetLength(tid) + coords[0];
        auto end = this.header.targetLength(tid) + coords[1];
        auto c = ZBHO(start, end);
        return query(tid, c);
    }

    /// opIndex with an int tid and a Coordinate and Offset
    /// bam[0,ZB(1) .. $]
    auto opIndex(int tid, Tuple!(OffsetType, OffsetType) coords)
    {
        auto start = this.header.targetLength(tid) + coords[0];
        auto end = this.header.targetLength(tid) + coords[1];
        auto c = ZBHO(start, end);
        return query(tid, c);
    }


    /// Return an InputRange representing all records in the SAM/BAM/CRAM
    auto allRecords()
    {
        auto newf = this.f.dup;
        newf.resetToFirstRecord;
        return newf.byRecord!Bam1.map!(x=>SAMRecord(bam_dup1(x), header));
    }

    deprecated("Avoid snake_case names")
    alias all_records = allRecords;

}

///
debug(dhtslib_unittest) unittest
{
    import dhtslib.sam;
    import htslib.hts_log : hts_log_info;
    import std.path : buildPath, dirName;
    import std.string : fromStringz;
    import std.array : array; 

    hts_set_log_level(htsLogLevel.HTS_LOG_WARNING);
    hts_log_info(__FUNCTION__, "Testing SAMFile & SAMRecord");
    hts_log_info(__FUNCTION__, "Loading test file");
    auto sam = SAMFile(buildPath(dirName(dirName(dirName(dirName(__FILE__)))),"htslib","test","auxf#values.sam"), 0);
    auto sam2 = SAMWriter("/tmp/test.bam", sam.header);
    auto readrange = sam.allRecords;
    hts_log_info(__FUNCTION__, "Getting read 1");
    assert(readrange.empty == false);
    auto read = readrange.front();
    
    // writeln(read.sequence);
    assert(read.sequence=="GCTAGCTCAG");
    assert(sam.allRecords.array.length == 2);
    sam2.write(read);
    destroy(sam2);


    // testing with multiple specified threads
    sam = SAMFile("/tmp/test.bam", 2);
    readrange = sam.allRecords;
    assert(readrange.empty == false);
    read = readrange.front();
    assert(read.sequence=="GCTAGCTCAG");
    assert(sam.allRecords.array.length == 1);

    // testing with no additional threads
    sam = SAMFile("/tmp/test.bam", 0);
    readrange = sam.allRecords;
    assert(readrange.empty == false);
    read = readrange.front();
    assert(read.sequence=="GCTAGCTCAG");
    assert(sam.allRecords.array.length == 1);

    // testing SAMReader targets/tid functions
    assert(sam.header.nTargets == 1);
    assert(sam.header.targetId("Sheila") == 0);
    assert(sam.header.targetLength(0) == 20);
    assert(sam.header.targetLengths == [20]);
    assert(sam.header.targetNames == ["Sheila"]);

}

///
debug(dhtslib_unittest) unittest
{
    import std.stdio;
    import dhtslib.sam;
    import dhtslib.sam.md : MDItr;
    import std.algorithm : map;
    import std.array : array;
    import std.path : buildPath,dirName;
    import std.range : drop;
    hts_set_log_level(htsLogLevel.HTS_LOG_WARNING);
    hts_log_info(__FUNCTION__, "Testing SAMFile query");
    hts_log_info(__FUNCTION__, "Loading test file");

    auto bam = SAMFile(buildPath(dirName(dirName(dirName(dirName(__FILE__)))),"htslib","test","range.bam"), 0);
    assert(bam.allRecords.array.length == 112);
    // assert(bam["CHROMOSOME_I"].array.length == 18);
    // assert(bam["CHROMOSOME_II"].array.length == 34);
    // assert(bam["CHROMOSOME_III"].array.length == 41);
    // assert(bam["CHROMOSOME_IV"].array.length == 19);
    // assert(bam["CHROMOSOME_V"].array.length == 0);
    assert(bam.query("CHROMOSOME_I", ZBHO(900, 2000)) .array.length == 14);
    assert(bam["CHROMOSOME_I",ZB(900) .. ZB(2000)].array.length == 14);
    assert(bam[0, ZB(900) .. ZB(2000)].array.length == 14);

    assert(bam["CHROMOSOME_I",ZB(940)].array.length == 2);
    assert(bam[0, ZB(940)].array.length == 2);


    assert(bam["CHROMOSOME_I",ZB(900) .. $].array.length == 18);
    assert(bam[0, ZB(900) .. $].array.length == 18);
    assert(bam["CHROMOSOME_I",$].array.length == 0);
    assert(bam[0, $].array.length == 0);
    assert(bam[["CHROMOSOME_I:900-2000","CHROMOSOME_II:900-2000"]].array.length == 33);

    assert(bam.query("CHROMOSOME_I", OBHO(901, 2000)) .array.length == 14);
    assert(bam["CHROMOSOME_I",OB(901) .. OB(2001)].array.length == 14);
    assert(bam[0, OB(901) .. OB(2001)].array.length == 14);

    assert(bam["CHROMOSOME_I",OB(941)].array.length == 2);
    assert(bam[0, OB(941)].array.length == 2);


    assert(bam["CHROMOSOME_I",OB(901) .. $].array.length == 18);
    assert(bam[0, OB(901) .. $].array.length == 18);
    assert(bam["CHROMOSOME_I",$].array.length == 0);
    assert(bam[0, $].array.length == 0);
    assert(bam[["CHROMOSOME_I:900-2000","CHROMOSOME_II:900-2000"]].array.length == 33);

    assert(bam["CHROMOSOME_II",$-1918 .. $].array.length == 0);
    assert(bam["CHROMOSOME_II", ZB(3082) .. $].array.length == 0);
    assert(bam["CHROMOSOME_II",$-1919 .. $].array.length == 1);
    assert(bam["CHROMOSOME_II", ZB(3081) .. $].array.length == 1);
    assert(bam["CHROMOSOME_II",$-2018 .. $].array.length == 2);

    auto range = bam[["CHROMOSOME_I:900-2000","CHROMOSOME_II:900-2000"]];
    auto range1 = range.save;
    range = range.drop(5);
    auto range2 = range.save;
    range = range.drop(5);
    auto range3 = range.save;
    range = range.drop(10);
    auto range4 = range.save;
    assert(range1.array.length == 33);
    assert(range2.array.length == 28);
    assert(range3.array.length == 23);
    assert(range4.array.length == 13);
}
debug(dhtslib_unittest) unittest
{
    import std.stdio;
    import dhtslib.sam;
    import std.array : array;
    import std.path : buildPath,dirName;
    import std.range : drop;
    hts_set_log_level(htsLogLevel.HTS_LOG_WARNING);
    hts_log_info(__FUNCTION__, "Testing AllRecordsRange save()");
    hts_log_info(__FUNCTION__, "Loading test file");

    auto bam = SAMFile(buildPath(dirName(dirName(dirName(dirName(__FILE__)))),"htslib","test","range.bam"), 0);
    assert(bam.allRecords.array.length == 112);
    assert(bam.allRecords.array.length == 112);

    auto range = bam.allRecords;

    auto range1 = range.save();
    range = range.drop(5);

    auto range2 = range.save();
    range = range.drop(5);

    auto range3 = range.save();
    range = range.drop(5);

    auto range4 = range.save();
    range = range.drop(50);

    auto range5 = range.save();
    range.popFront;

    auto range6 = range.save();

    assert(range1.array.length == 112);
    assert(range2.array.length == 107);
    assert(range3.array.length == 102);
    assert(range4.array.length == 97);
    assert(range5.array.length == 47);
    assert(range6.array.length == 46);
}
///
debug(dhtslib_unittest) unittest
{
    import std.stdio;
    import dhtslib.sam;
    import std.array : array;
    import std.path : buildPath,dirName;
    import std.range : drop;
    hts_set_log_level(htsLogLevel.HTS_LOG_WARNING);
    hts_log_info(__FUNCTION__, "Testing RecordsRange save()");
    hts_log_info(__FUNCTION__, "Loading test file");

    auto bam = SAMReader(buildPath(dirName(dirName(dirName(dirName(__FILE__)))),"htslib","test","range.bam"), 4);
    
    auto range = bam.query("CHROMOSOME_I", ZBHO(900, 2000));
    assert(bam.query("CHROMOSOME_I", ZBHO(900, 2000)).array.length == 14);
    
    auto range1 =  range.save;
    range = range.drop(1);
    
    auto range2 =  range.save;
    range = range.drop(2);
    
    auto range3 =  range.save;
    range = range.drop(3);
    
    auto range4 =  range.save;
    range = range.drop(5);
    
    auto range5 =  range.save;
    range.popFront;
    
    auto range6 = range.save;

    assert(range1.array.length == 14);
    assert(range2.array.length == 13);
    assert(range3.array.length == 11);
    assert(range4.array.length == 8);
    assert(range5.array.length == 3);
    assert(range6.array.length == 2);
}