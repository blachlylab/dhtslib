module dhtslib.sam.reader;

import core.stdc.stdio : SEEK_SET;
import std.format;
import std.stdio : writeln, writefln, stderr, File;
import std.string : fromStringz, toStringz;
import std.typecons : Tuple;
import std.range.interfaces : InputRange, inputRangeObject;
import std.range.primitives : isInputRange, ElementType;

import htslib.hts : htsFile, hts_open, hts_close, hts_hopen, hts_set_threads, htsExactFormat;
import htslib.hts : hts_idx_t, hts_itr_t, hts_itr_multi_t, hts_reglist_t, hts_pair32_t;
import htslib.hfile : hdopen, hclose, hFILE, off_t, htell, hseek;
import htslib.bgzf : bgzf_tell, bgzf_seek;

import htslib.hts_log;
import htslib.kstring;
import htslib.sam;
import dhtslib.sam.record;
import dhtslib.coordinates;
import dhtslib.sam.header;
import dhtslib.sam : cmpInterval, cmpRegList;

alias SAMFile = SAMReader;
/**
Encapsulates a SAM/BAM file.

Implements InputRange interface using htslib calls.
If indexed, Random-access query via multidimensional slicing.
*/
struct SAMReader
{
    /// filename; as usable from D
    string filename;

    /// filename \0-terminated C string; reference needed to avoid GC reaping result of toStringz when ctor goes out of scope
    private const(char)* fn;

    /// htsFile
    private htsFile* fp;

    /// hFILE if required
    private hFILE* f;

    /// header struct
    SAMHeader header;

    private off_t header_offset;

    /// SAM/BAM/CRAM index 
    private hts_idx_t* idx;

    private kstring_t line;

    /// disallow copying
    @disable this(this);

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
        static if (is(T == string))
        {
            this.filename = f;
            this.fn = toStringz(f);
            this.fp = hts_open(this.fn, cast(immutable(char)*) "r");
        }
        else static if (is(T == File))
        {
            this.filename = f.name();
            this.fn = toStringz(f.name);
            this.f = hdopen(f.fileno, cast(immutable(char)*) "r");
            this.fp = hts_hopen(this.f, this.fn, cast(immutable(char)*) "r");
        }
        else assert(0);

        if (extra_threads == -1)
        {
            if ( totalCPUs > 1)
            {
                hts_log_info(__FUNCTION__,
                        format("%d CPU cores detected; enabling multithreading", totalCPUs));
                // hts_set_threads adds N _EXTRA_ threads, so totalCPUs - 1 seemed reasonable,
                // but overcomitting by 1 thread (i.e., passing totalCPUs) buys an extra 3% on my 2-core 2013 Mac
                hts_set_threads(this.fp, totalCPUs);
            }
        }
        else if (extra_threads > 0)
        {
            if ((extra_threads + 1) > totalCPUs)
                hts_log_warning(__FUNCTION__, "More threads requested than CPU cores detected");
            hts_set_threads(this.fp, extra_threads);
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
        this.header = SAMHeader(sam_hdr_read(this.fp));
        this.idx = sam_index_load(this.fp, this.fn);
        if (this.idx == null)
        {
            hts_log_info(__FUNCTION__, "SAM index not found");
            // TODO: attempt to build
            // TODO: edit range to return empty immediately if no idx
        }
        hts_log_debug(__FUNCTION__, format("SAM index: %s", this.idx));

        // collect header offset just in case it is a sam file
        // we would like to iterate
        if(this.fp.is_bgzf) this.header_offset = bgzf_tell(this.fp.fp.bgzf);
        else this.header_offset = htell(this.fp.fp.hfile);
    }

    ~this()
    {
        debug (dhtslib_debug)
        {
            writeln("SAMFile dtor");
        }
        
        // SAMHeader now takes care of this
        // bam_hdr_destroy(this.header);

        if((this.fp !is null) && (this.f is null))
        {
            const auto ret = hts_close(fp);
            if (ret < 0)
                writefln("There was an error closing %s", fromStringz(this.fn));
        }
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
    auto query(CoordSystem cs)(string chrom, Coordinates!cs coords)
    in (!this.header.isNull)
    {
        auto tid = this.header.targetId(chrom);
        return query(tid, coords);
    }

    /// ditto
    auto query(CoordSystem cs)(int tid, Coordinates!cs coords)
    in (!this.header.isNull)
    {
        /// convert to zero-based half-open
        auto newcoords = coords.to!(CoordSystem.zbho);
        auto itr = sam_itr_queryi(this.idx, tid, newcoords.start, newcoords.end);
        return RecordRange(this.fp, this.header, itr);
    }

    /// ditto
    auto query(CoordSystem cs)(ChromCoordinates!cs region)
    in (!this.header.isNull)
    {
        return query(region.chrom, region.coords);
    }

    /// ditto
    auto query(string[] regions)
    in (!this.header.isNull)
    {
        return RecordRangeMulti(this.fp, this.idx, this.header, &(this), regions);
    }

    /// ditto
    auto opIndex(CoordSystem cs)(ChromCoordinates!cs region)
    {
        return query!cs(region);
    }

    /// ditto
    auto opIndex(string[] regions)
    {
        return query(regions);
    }

    /// ditto
    auto opIndex(Basis bs)(string tid, Coordinate!bs[] pos)
    {
        assert(pos.length == 2);
        auto coords = Coordinates!(getCoordinateSystem!(bs,End.open))(pos[0], pos[1]);
        return query(tid, coords);
    }

    /// ditto
    auto opIndex(Basis bs)(int tid, Coordinate!bs[] pos)
    {
        assert(pos.length == 2);
        auto coords = Coordinates!(getCoordinateSystem!(bs,End.open))(pos[0], pos[1]);
        return query(tid, coords);
    }

    /// ditto
    auto opIndex(Basis bs)(string tid, Coordinate!bs pos)
    {
        auto coords = Coordinates!(getCoordinateSystem!(bs,End.open))(pos, pos + 1);
        return query(tid, coords);
    }

    /// ditto
    auto opIndex(Basis bs)(int tid, Coordinate!bs pos)
    {
        auto coords = Coordinates!(getCoordinateSystem!(bs,End.open))(pos, pos + 1);
        return query(tid, coords);
    }

    /// ditto
    deprecated("use multidimensional slicing with second parameter as range ([\"chr1\", 1 .. 2])")
    auto opIndex(Basis bs)(string tid, Coordinate!bs pos1, Coordinate!bs pos2)
    {
        auto coords = Coordinates!(getCoordinateSystem!(bs,End.open))(pos1, pos2);
        return query(tid, coords);
    }

    /// ditto
    deprecated("use multidimensional slicing with second parameter as range ([20, 1 .. 2])")
    auto opIndex(Basis bs)(int tid, Coordinate!bs pos1, Coordinate!bs pos2)
    {
        auto coords = Coordinates!(getCoordinateSystem!(bs,End.open))(pos1, pos2);
        return query(tid, coords);
    }

    /// ditto
    auto opSlice(size_t dim, Basis bs)(Coordinate!bs start, Coordinate!bs end) if(dim == 1)
    {
        assert(end > start);
        return [start, end];
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
    /// ditto
    auto opIndex(string ctg, OffsetType endoff)
    {
        auto tid = this.header.targetId(ctg);
        auto end = this.header.targetLength(tid) + endoff.offset;
        // TODO review: is targetLength the last nt, or targetLength - 1 the last nt?
        auto coords = Coordinates!(CoordSystem.zbho)(end, end + 1);
        return query(tid, coords);
    }
    /// ditto
    auto opIndex(int tid, OffsetType endoff)
    {
        auto end = this.header.targetLength(tid) + endoff.offset;
        // TODO review: is targetLength the last nt, or targetLength - 1 the last nt?
        auto coords = Coordinates!(CoordSystem.zbho)(end, end + 1);
        return query(tid, coords);
    }
    /// ditto
    auto opSlice(size_t dim, Basis bs)(Coordinate!bs start, OffsetType off) if(dim == 1)
    {
        return Tuple!(Coordinate!bs, OffsetType)(start, off);
    }
    /// ditto
    auto opIndex(Basis bs)(string ctg, Tuple!(Coordinate!bs, OffsetType) coords)
    {
        auto tid = this.header.targetId(ctg);
        auto end = this.header.targetLength(tid) + coords[1];
        auto endCoord = ZeroBased(end);
        auto newEndCoord = endCoord.to!bs;
        auto c = Coordinates!(getCoordinateSystem!(bs,End.open))(coords[0], newEndCoord);
        return query(tid, c);
    }
    /// ditto
    auto opIndex(Basis bs)(int tid, Tuple!(Coordinate!bs, OffsetType) coords)
    {
        auto end = this.header.targetLength(tid) + coords[1];
        auto endCoord = ZeroBased(end);
        auto newEndCoord = endCoord.to!bs;
        auto c = Coordinates!(getCoordinateSystem!(bs,End.open))(coords[0], newEndCoord);
        return query(tid, c);
    }


    /// Return an InputRange representing all records in the SAM/BAM/CRAM
    InputRange!SAMRecord allRecords()
    {
        if (this.fp.format.format == htsExactFormat.sam)
            return SamRecordsRange(this.fp, this.header, this.header_offset).inputRangeObject;
        else if(this.idx == null){
            return BamRecordsRange(this.fp, this.header, this.header_offset).inputRangeObject;
        }else{
            assert(this.fp.format.format == htsExactFormat.bam);
            //auto range = AllRecordsRange(this.fp, this.header);
            import htslib.hts : HTS_IDX_START;
            auto itr = sam_itr_queryi(this.idx, HTS_IDX_START, 0, 0);
            return RecordRange(this.fp, this.header, itr).inputRangeObject;
        }
    }

    deprecated("Avoid snake_case names")
    alias all_records = allRecords;

    /// Iterate through all records in the BAM or BGZF compressed SAM with no index
    struct BamRecordsRange
    {
        private htsFile*    fp;     // belongs to parent; shared
        private off_t header_offset;
        private SAMHeader header; // belongs to parent; shared
        private bam1_t*     b;
        private bool initialized;   // Needed to support both foreach and immediate .front()
        private int success;        // sam_read1 return code

        this(htsFile* fp, SAMHeader header, off_t header_offset)
        {
            this.fp = fp;
            this.header_offset = header_offset;
            this.header = header;
            this.b = bam_init1();
            
            // Need to seek past header offset
            bgzf_seek(fp.fp.bgzf, this.header_offset, SEEK_SET);
        }

        ~this()
        {
            //debug(dhtslib_debug) hts_log_debug(__FUNCTION__, "dtor");
            //TODO ?: free(this.b);
        }

        /// InputRange interface
        @property bool empty() // @suppress(dscanner.suspicious.incorrect_infinite_range)
        {
            if (!this.initialized) {
                this.popFront();
                this.initialized = true;
            }
            if (success >= 0)
                return false;
            else if (success == -1)
                return true;
            else
            {
                hts_log_error(__FUNCTION__, "*** ERROR sam_read1 < -1");
                return true;
            }
        }
        /// ditto
        void popFront()
        {
            success = sam_read1(this.fp, this.header.h, this.b);
            //bam_destroy1(this.b);
        }
        /// ditto
        SAMRecord front()
        {
            assert(this.initialized, "front called before empty");
            return SAMRecord(bam_dup1(this.b), this.header);
        }

        SAMRecord moveFront()
        {
            auto ret = SAMRecord(bam_dup1(this.b), this.header);
            popFront;
            return ret;
        }

    }

    /// Iterate through all records in the SAM
    struct SamRecordsRange
    {
        private htsFile*    fp;     // belongs to parent; shared
        private SAMHeader  header; // belongs to parent; shared
        private off_t header_offset;
        private bam1_t*     b;
        private bool initialized;   // Needed to support both foreach and immediate .front()
        private int success;        // sam_read1 return code

        this(htsFile* fp, SAMHeader header, off_t header_offset)
        {
            this.fp = fp;
            this.header_offset = header_offset;
            assert(this.fp.format.format == htsExactFormat.sam);
            this.header = header;
            this.b = bam_init1();

            // Need to seek past header offset
            hseek(fp.fp.hfile, this.header_offset, SEEK_SET);
        }

        ~this()
        {
            //debug(dhtslib_debug) hts_log_debug(__FUNCTION__, "dtor");
            //TODO ?: free(this.b);
        }

        /// InputRange interface
        @property bool empty() // @suppress(dscanner.suspicious.incorrect_infinite_range)
        {
            if (!this.initialized) {
                this.popFront();
                this.initialized = true;
            }
            if (success >= 0)
                return false;
            else if (success == -1)
                return true;
            else
            {
                hts_log_error(__FUNCTION__, "*** ERROR sam_read1 < -1");
                return true;
            }
        }
        /// ditto
        void popFront()
        {
            success = sam_read1(this.fp, this.header.h, this.b);
            //bam_destroy1(this.b);
        }
        /// ditto
        SAMRecord front()
        {
            assert(this.initialized, "front called before empty");
            return SAMRecord(bam_dup1(this.b), this.header);
        }

        SAMRecord moveFront()
        {
            auto ret = SAMRecord(bam_dup1(this.b), this.header);
            popFront;
            return ret;
        }

    }

    /// Iterate over records falling within a queried region (TODO: itr_multi_query)
    /// TODO destroy the itr with dtor
    struct RecordRange
    {
        private htsFile* fp;
        private SAMHeader h;
        private hts_itr_t* itr;
        private bam1_t* b;

        private int r;  // sam_itr_next: >= 0 on success; -1 when there is no more data; < -1 on error

        /// Constructor relieves caller of calling bam_init1 and simplifies first-record flow 
        this(htsFile* fp, SAMHeader header, hts_itr_t* itr)
        {
            this.fp = fp;
            this.h = header;
            this.itr = itr;
            this.b = bam_init1();
            
            //assert(itr !is null, "itr was null");
 
            if (this.itr !is null)
                popFront();
        }

        /// InputRange interface
        @property bool empty()
        {
            // TODO, itr.finished shouldn't be used
            if (this.itr is null) return true;
            return (r < 0 || itr.finished) ? true : false;
        }
        /// ditto
        void popFront()
        {
            this.r = sam_itr_next(this.fp, this.itr, this.b);
        }
        /// ditto
        SAMRecord front()
        {
            return SAMRecord(bam_dup1(b), this.h);
        }

        SAMRecord moveFront()
        {
            auto ret = SAMRecord(bam_dup1(this.b), this.h);
            popFront;
            return ret;
        }
    }

    static assert(isInputRange!RecordRange);
    static assert(is(ElementType!RecordRange == SAMRecord));

    /// Iterate over records falling within queried regions using a RegionList
    struct RecordRangeMulti
    {
        private htsFile* fp;
        private SAMHeader h;
        private hts_itr_multi_t* itr;
        private bam1_t* b;
        private hts_reglist_t[] rlist;
        private int r;

        ///
        this(htsFile* fp, hts_idx_t* idx, SAMHeader header, SAMFile* sam, string[] regions)
        {
            rlist = RegionList(sam, regions).getRegList();
            this.fp = fp;
            this.h = header;
            b = bam_init1();
            itr = sam_itr_regions(idx, this.h.h, rlist.ptr, cast(uint) rlist.length);
            //debug(dhtslib_debug) { writeln("sam_itr null? ",(cast(int)itr)==0); }
            hts_log_debug(__FUNCTION__, format("SAM itr null?: %s", cast(int) itr == 0));
            popFront();
        }

        /// InputRange interface
        @property bool empty()
        {
            return (r <= 0 && itr.finished) ? true : false;
        }
        /// ditto
        void popFront()
        {
            r = sam_itr_multi_next(this.fp, this.itr, this.b);
        }
        /// ditto
        SAMRecord front()
        {
            return SAMRecord(bam_dup1(b), this.h);
        }
    }

    /// List of regions based on sam/bam
    struct RegionList
    {
        import std.algorithm.iteration : splitter;
        import std.algorithm.sorting : sort;
        import std.range : drop, array;
        import std.conv : to;

        private hts_reglist_t[string] rlist;
        private SAMFile* sam;

        ///
        this(SAMFile* sam, string[] queries)
        {
            this.sam = sam;
            foreach (q; queries)
            {
                addRegion(q);
            }
        }

        /// Add a region in standard format (chr1:2000-3000) to the RegionList
        void addRegion(string reg)
        {
            //chr:1-2
            //split into chr and 1-2
            auto split = reg.splitter(":");
            auto chr = split.front;
            //split 1-2 into 1 and 2
            split = split.drop(1).front.splitter("-");
            if (sam.header.targetId(chr) < 0)
            {
                hts_log_error(__FUNCTION__, "tid not present in sam/bam");
            }
            addRegion(sam.header.targetId(chr), split.front.to!int, split.drop(1).front.to!int);
        }

        /// Add a region by {target/contig/chr id, start coord, end coord} to the RegionList
        void addRegion(int tid, int beg, int end)
        {
            if (tid > this.sam.header.nTargets || tid < 0)
                hts_log_error(__FUNCTION__, "tid not present in sam/bam");

            auto val = (this.sam.header.targetNames[tid] in this.rlist);
            hts_pair32_t p;

            if (beg < 0)
                hts_log_error(__FUNCTION__, "first coordinate < 0");

            if (beg >= this.sam.header.targetLength(tid))
                hts_log_error(__FUNCTION__, "first coordinate larger than tid length");

            if (end < 0)
                hts_log_error(__FUNCTION__, "second coordinate < 0");

            if (end >= this.sam.header.targetLength(tid))
                hts_log_error(__FUNCTION__, "second coordinate larger than tid length");

            p.beg = beg;
            p.end = end;
            hts_pair32_t[] plist;
            if (val is null)
            {
                hts_reglist_t r;

                //set tid
                r.tid = tid;

                //create intervals
                plist = plist ~ p;
                r.intervals = plist.ptr;
                r.count = cast(uint) plist.length;
                r.min_beg = p.beg;
                r.max_end = p.end;
                this.rlist[this.sam.header.targetNames[tid]] = r;
            }
            else
            {
                plist = (val.intervals[0 .. val.count] ~ p).sort!(cmpInterval).array;
                val.intervals = plist.ptr;
                val.count = cast(uint) plist.length;
                val.min_beg = plist[0].beg;
                val.max_end = plist[$ - 1].end;
            }
        }

        hts_reglist_t[] getRegList()
        {
            return rlist.byValue.array.sort!(cmpRegList).array;
        }
    }
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
    
    writeln(read.sequence);
    assert(read.sequence=="GCTAGCTCAG");
    assert(sam.allRecords.array.length == 2);
    sam2.write(read);
    sam2.close;


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
    assert(bam.query(Obc("CHROMOSOME_I:900-2000")).array.length == 14);
    assert(bam.query("CHROMOSOME_I", Zbho(900, 2000)) .array.length == 14);
    assert(bam["CHROMOSOME_I",ZeroBased(900) .. ZeroBased(2000)].array.length == 14);
    assert(bam["CHROMOSOME_I",[ZeroBased(900), ZeroBased(2000)]].array.length == 14);
    assert(bam[0, [ZeroBased(900), ZeroBased(2000)]].array.length == 14);

    assert(bam["CHROMOSOME_I",ZeroBased(940)].array.length == 2);
    assert(bam[0, ZeroBased(940)].array.length == 2);


    assert(bam["CHROMOSOME_I",ZeroBased(900) .. $].array.length == 18);
    assert(bam[0, ZeroBased(900) .. $].array.length == 18);
    assert(bam["CHROMOSOME_I",$].array.length == 0);
    assert(bam[0, $].array.length == 0);
    assert(bam[["CHROMOSOME_I:900-2000","CHROMOSOME_II:900-2000"]].array.length == 33);
}