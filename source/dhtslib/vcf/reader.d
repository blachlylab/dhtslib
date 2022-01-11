module dhtslib.vcf.reader;

import std.string: fromStringz, toStringz;
import std.algorithm : map;
import std.utf: toUTFz;
import std.traits : ReturnType;

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
    HtslibFile     f;    /// HtslibFile
    VCFHeader   vcfhdr;    /// rc header wrapper
    UnpackLevel MAX_UNPACK;     /// see htslib.vcf

    HtslibIterator!Bcf1 range; /// HtslibIterator

    
    /// read existing VCF file
    /// MAX_UNPACK: setting alternate value could speed reading
    this(string fn, UnpackLevel MAX_UNPACK = UnpackLevel.None)
    {
        if (fn == "") throw new Exception("Empty filename passed to VCFReader constructor");
        this.f = HtslibFile(fn);
        if (!this.f.fp) throw new Exception("Could not hts_open file");
        
        this.f.setExtraThreads(1);    // extra decoding thread

        this.f.loadHeader;
        this.vcfhdr = VCFHeader(this.f.bcfHdr);
        this.range = this.f.byRecord!Bcf1;
    }

    invariant
    {
        assert(this.vcfhdr.hdr != null);
    }

    VCFHeader getHeader()
    {
        return this.vcfhdr;
    }

    /// InputRange interface: iterate over all records
    @property bool empty()
    {
        return this.range.empty;
        
    }
    /// ditto
    void popFront()
    {
        this.range.popFront;

    }
    /// ditto
    VCFRecord front()
    {
        // note that VCFRecord's constructor will be responsible for
        // * UnpackLeveling and
        // * destroying
        // its copy
        return VCFRecord(this.vcfhdr, bcf_dup(this.range.front), this.MAX_UNPACK);
    }

    typeof(this) save()
    {
        typeof(this) newRange = this;
        newRange.f = f;
        newRange.vcfhdr = vcfhdr;
        newRange.MAX_UNPACK = MAX_UNPACK;
        newRange.range = this.range.save;
        return newRange;
    }

    auto query(CoordSystem cs)(string chrom, Interval!cs coords)
    {
        auto newCoords = coords.to!ZBHO();

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
        auto tid = bcf_hdr_name2id(this.vcfhdr.hdr,toStringz(chrom));
        import std.stdio;

        return newF.query!Bcf1(tid, newCoords.start, newCoords.end)
                .map!(x => VCFRecord(this.vcfhdr, bcf_dup(x), this.MAX_UNPACK));
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
    assert(vcf.count == 14);
    vcf = VCFReader(buildPath(dirName(dirName(dirName(dirName(__FILE__)))),"htslib","test","tabix","vcf_file.vcf"));
    VCFRecord rec = vcf.front;
    assert(rec.chrom == "1");
    assert(rec.pos == 3000149);
    assert(rec.refAllele == "C");
    assert(rec.altAllelesAsArray == ["T"]);
    assert(rec.allelesAsArray == ["C","T"]);
    assert(approxEqual(rec.qual,59.2));
    assert(rec.filter == "PASS");
    
    vcf.popFront;
    rec = vcf.front;

    assert(rec.chrom == "1");
    assert(rec.pos == 3000150);
    assert(rec.refAllele == "C");
    assert(rec.altAllelesAsArray == ["T"]);
    assert(rec.allelesAsArray == ["C","T"]);
    assert(approxEqual(rec.qual,59.2));
    assert(rec.filter == "PASS");
    
    vcf = VCFReader(buildPath(dirName(dirName(dirName(dirName(__FILE__)))),"htslib","test","tabix","vcf_file.vcf"));
    auto range1 = vcf.save;
    vcf.popFront;
    auto range2 = vcf.save;
    vcf.popFront;
    vcf.popFront;
    auto range3 = vcf.save;
    vcf.popFront;
    vcf.popFront;
    vcf.popFront;
    auto range4 = vcf.save;
    assert(range1.array.length == 14);
    assert(range2.array.length == 13);
    assert(range3.array.length == 11);
    assert(range4.array.length == 8);
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
    auto reg = getIntervalFromString("1:3000151-3062916");
    
    auto vcf = VCFReader(fn);
    auto range = vcf.query(reg.contig, reg.interval);
    assert(range.count == 3);
    range = vcf.query(reg.contig, reg.interval);
    VCFRecord rec = range.front;
    
    assert(rec.chrom == "1");
    assert(rec.pos == 3000150);
    assert(rec.refAllele == "C");
    assert(rec.altAllelesAsArray == ["T"]);
    assert(rec.allelesAsArray == ["C","T"]);
    assert(approxEqual(rec.qual, 59.2));
    assert(rec.filter == "PASS");

    range.popFront;
    rec = range.front;

    assert(rec.chrom == "1");
    assert(rec.pos == 3062914);
    assert(rec.refAllele == "GTTT");
    assert(rec.altAllelesAsArray == ["G"]);
    assert(rec.allelesAsArray == ["GTTT","G"]);
    assert(approxEqual(rec.qual,12.9));
    assert(rec.filter == "q10");
    // assert(rec. == ["C","T"]);

    range = vcf.query(reg.contig, reg.interval);
    auto range1 = range.save;
    range.popFront;
    auto range2 = range.save;
    range.popFront;
    auto range3 = range.save;
    range.popFront;
    auto range4 = range.save;
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
    auto range = vcf.query(reg.contig, reg.interval);

    assert(range.count == 3);
    range = vcf.query(reg.contig, reg.interval);
    VCFRecord rec = range.front;

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
    rec = range.front;
    assert(rec.getFormats["DP"].to!int.array == [[32], [32]]);
    assert(rec.line.unpacked == UnpackLevel.Format);
    
}