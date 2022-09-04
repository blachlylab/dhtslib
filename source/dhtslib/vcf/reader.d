module dhtslib.vcf.reader;

import std.string: fromStringz, toStringz;
import std.utf: toUTFz;
import std.traits : ReturnType;
import std.parallelism : totalCPUs;
import std.format : format;

import dhtslib.vcf;
import dhtslib.tabix;
import dhtslib.coordinates;
import dhtslib.memory;
import dhtslib.util;
import dhtslib.file;

import htslib.vcf;
import htslib.hts_log;
import htslib.kstring;

alias BCFReader = VCFReader;

/** Basic support for reading VCF, BCF files */
struct VCFReader
{
    // htslib data structures
    HtslibFile     file;    /// rc htsFile wrapper
    VCFHeader   header;    /// rc header wrapper
    UnpackLevel MAX_UNPACK;     /// see htslib.vcf
    HtsIdx idx;
    Tbx tbx;

    @disable this();
    /// read existing VCF file
    /// MAX_UNPACK: setting alternate value could speed reading
    this(F)(F f, int threads = -1, UnpackLevel MAX_UNPACK = UnpackLevel.None)
    {
        this.file = HtslibFile(f);
        
        if (threads == -1)
        {
            hts_log_info(__FUNCTION__,
                        format("%d CPU cores detected; enabling multithreading", totalCPUs));
        } else if (threads > totalCPUs)
            hts_log_warning(__FUNCTION__, "More threads requested than CPU cores detected");
        else if (threads == 0)
            hts_log_debug(__FUNCTION__, "Zero threads requested");
        else
            hts_log_warning(__FUNCTION__, "Invalid negative number of extra threads requested");

        if (threads > 0 || threads == -1) {
            this.file.setThreads(threads);
        }

        this.file.loadHeader;
        this.header = VCFHeader(this.file.bcfHdr);
    }

    /// Explicit postblit to avoid 
    /// https://github.com/blachlylab/dhtslib/issues/122
    this(this)
    {
        this.file = file;
        this.header = header;
        this.MAX_UNPACK = MAX_UNPACK;
    }

    invariant
    {
        assert(this.header.hdr != null);
    }

    VCFHeader getHeader()
    {
        return this.header;
    }

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
    {
        auto tid = this.header.targetId(chrom);
        return query(tid, coords);
    }

    /// ditto
    auto query(CoordSystem cs)(int tid, Interval!cs coords)
    {
        /// convert to zero-based half-open
        auto newcoords = coords.to!(CoordSystem.zbho);

        /// load index
        if(this.idx == null)
            this.idx = this.file.loadHtsIndex;
        if (this.idx == null)
        {
            if(this.tbx == null)
                this.tbx = this.file.loadTabixIndex;
            if (this.idx == null && this.tbx == null) {
                this.idx = this.file.loadHtsIndex;
                hts_log_error(__FUNCTION__, "Couldn't find tbx index or hts index");
            }
        }
        auto fcopy = this.file.dup;
        fcopy.resetToFirstRecord;
        return VCFReaderItr(this.header, fcopy.query!Bcf1(tid, newcoords.start.pos, newcoords.end.pos), this.MAX_UNPACK);
    }


    /// ditto
    auto query(string[] regions)
    {
        /// load index
        if(this.idx == null)
            this.idx = this.file.loadHtsIndex;
        if (this.idx == null)
        {
            hts_log_error(__FUNCTION__, "SAM index not found");
            throw new Exception("SAM index not found");
        }
        auto fcopy = this.file.dup;
        fcopy.resetToFirstRecord;
        return fcopy.query(regions);
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

    /// Return an InputRange representing all records in the SAM/BAM/CRAM
    auto allRecords()
    {
        this.file.resetToFirstRecord;
        return VCFReaderItr(this.header, this.file.byRecord!Bcf1(), this.MAX_UNPACK);
    }

    struct VCFReaderItr {
        VCFHeader hdr;
        HtslibIterator!Bcf1 itr;
        UnpackLevel MAX_UNPACK;

        this(VCFHeader hdr, HtslibIterator!Bcf1 itr, UnpackLevel MAX_UNPACK) {
            this.hdr = hdr;
            this.itr = itr;
            this.MAX_UNPACK = MAX_UNPACK;
        }

        auto front() {
            return VCFRecord(this.hdr, this.itr.front, MAX_UNPACK);
        }

        void popFront() {
            this.itr.popFront;
        }

        auto save() {
            return VCFReaderItr(this.hdr, this.itr.save, this.MAX_UNPACK);
        }

        bool empty() {
            return this.itr.empty;
        }
    }
}

debug(dhtslib_unittest) unittest
{
    import std.stdio;
    import htslib.hts_log;
    import std.algorithm : map, count;
    import std.array : array;
    import std.path : buildPath, dirName;
    import std.math : approxEqual;
    hts_set_log_level(htsLogLevel.HTS_LOG_INFO);
    hts_log_info(__FUNCTION__, "Testing VCFReader (Tabix)");
    hts_log_info(__FUNCTION__, "Loading test file");

    auto vcf = VCFReader(buildPath(dirName(dirName(dirName(dirName(__FILE__)))),"htslib","test","tabix","vcf_file.vcf"));
    auto recs = vcf.allRecords;
    assert(recs.count == 15);
    vcf = VCFReader(buildPath(dirName(dirName(dirName(dirName(__FILE__)))),"htslib","test","tabix","vcf_file.vcf"));
    recs = vcf.allRecords;
    VCFRecord rec = recs.front;
    assert(rec.chrom == "1");
    assert(rec.pos == 3000149);
    assert(rec.refAllele == "C");
    assert(rec.altAllelesAsArray == ["T"]);
    assert(rec.allelesAsArray == ["C","T"]);
    assert(approxEqual(rec.qual,59.2));
    assert(rec.filter == "PASS");
    
    recs.popFront;
    rec = recs.front;

    assert(rec.chrom == "1");
    assert(rec.pos == 3000150);
    assert(rec.refAllele == "C");
    assert(rec.altAllelesAsArray == ["T"]);
    assert(rec.allelesAsArray == ["C","T"]);
    assert(approxEqual(rec.qual,59.2));
    assert(rec.filter == "PASS");
    
    vcf = VCFReader(buildPath(dirName(dirName(dirName(dirName(__FILE__)))),"htslib","test","tabix","vcf_file.vcf"));
    recs = vcf.allRecords;
    auto range1 = recs.save;
    recs.popFront;
    auto range2 = recs.save;
    recs.popFront;
    recs.popFront;
    auto range3 = recs.save;
    recs.popFront;
    recs.popFront;
    recs.popFront;
    auto range4 = recs.save;
    assert(range1.array.length == 15);
    assert(range2.array.length == 14);
    assert(range3.array.length == 12);
    assert(range4.array.length == 9);
}

debug(dhtslib_unittest) unittest
{
    import dhtslib.util;
    import std.stdio;
    import htslib.hts_log;
    import std.algorithm : map, count;
    import std.array : array;
    import std.path : buildPath, dirName;
    import std.math : approxEqual;
    hts_set_log_level(htsLogLevel.HTS_LOG_INFO);
    hts_log_info(__FUNCTION__, "Testing VCFReader");
    hts_log_info(__FUNCTION__, "Loading test file");
    auto fn = buildPath(dirName(dirName(dirName(dirName(__FILE__)))),"htslib","test","tabix","vcf_file.vcf.gz");
    auto tbx = TabixIndexedFile(fn);
    auto reg = getIntervalFromString("1:3000151-3062916");
    auto vcf = VCFReader(fn);
    auto recs = vcf.query(reg.contig, reg.interval);
    assert(recs.count == 3);
    recs = vcf.query(reg.contig, reg.interval);
    VCFRecord rec = recs.front;
    assert(rec.chrom == "1");
    assert(rec.pos == 3000150);
    assert(rec.refAllele == "C");
    assert(rec.altAllelesAsArray == ["T"]);
    assert(rec.allelesAsArray == ["C","T"]);
    assert(approxEqual(rec.qual, 59.2));
    assert(rec.filter == "PASS");
    assert(vcf.header.getSamples == ["A", "B"]);

    recs.popFront;
    rec = recs.front;

    assert(rec.chrom == "1");
    assert(rec.pos == 3062914);
    assert(rec.refAllele == "GTTT");
    assert(rec.altAllelesAsArray == ["G"]);
    assert(rec.allelesAsArray == ["GTTT","G"]);
    assert(approxEqual(rec.qual,12.9));
    assert(rec.filter == "q10");
    // assert(rec. == ["C","T"]);

    recs = vcf.query(reg.contig, reg.interval);
    auto range1 = recs.save;
    recs.popFront;
    auto range2 = recs.save;
    recs.popFront;
    auto range3 = recs.save;
    recs.popFront;
    auto range4 = recs.save;
    assert(range1.array.length == 3);
    assert(range2.array.length == 2);
    assert(range3.array.length == 1);
    assert(range4.array.length == 0);
}

debug(dhtslib_unittest) unittest
{
    import dhtslib.util;
    import std.stdio;
    import htslib.hts_log;
    import std.algorithm : map, count;
    import std.array : array;
    import std.path : buildPath, dirName;
    import std.math : approxEqual;
    hts_set_log_level(htsLogLevel.HTS_LOG_INFO);
    hts_log_info(__FUNCTION__, "Testing bcf1_t Unpacking");
    hts_log_info(__FUNCTION__, "Loading test file");
    auto fn = buildPath(dirName(dirName(dirName(dirName(__FILE__)))),"htslib","test","tabix","vcf_file.vcf.gz");
    auto reg = getIntervalFromString("1:3000151-3062916");

    auto vcf = VCFReader(fn);
    auto recs = vcf.query(reg.contig, reg.interval);
    assert(recs.count == 3);
    vcf.MAX_UNPACK = UnpackLevel.None;
    recs = vcf.query(reg.contig, reg.interval);
    VCFRecord rec = recs.front;

    assert(rec.chrom == "1");
    assert(rec.pos == 3000150);

    assert(approxEqual(rec.qual, 59.2));
    assert(rec.line.unpacked == UnpackLevel.None);

    assert(rec.id == ".");
    assert(rec.line.unpacked == UnpackLevel.AltAllele);

    assert(rec.refAllele == "C");
    assert(rec.line.unpacked == UnpackLevel.AltAllele);

    assert(rec.altAllelesAsArray == ["T"]);
    assert(rec.line.unpacked == UnpackLevel.AltAllele);

    assert(rec.filter == "PASS");
    assert(rec.line.unpacked == (UnpackLevel.Filter | UnpackLevel.AltAllele));

    assert(rec.getInfos["AN"].to!int == 4);
    assert(rec.line.unpacked == UnpackLevel.SharedFields);

    assert(rec.getFormats["DP"].to!int.array == [[32], [32]]);
    assert(rec.line.unpacked == UnpackLevel.All);

    /// only one extra test for pulling only format fields
    ///
    /// bcf_unpack automatically promotes UnpackLevel.Filter to
    /// (UnpackLevel.Filter | UnpackLevel.AltAllele)
    ///
    /// bcf_unpack automatically promotes UnpackLevel.Info to
    /// (UnpackLevel.Info | UnpackLevel.SharedFields)

    /// since front calls bcf_dup, we get a fresh unpacked record
    rec = recs.front;
    assert(rec.getFormats["DP"].to!int.array == [[32], [32]]);
    assert(rec.line.unpacked == UnpackLevel.Format);
    
}