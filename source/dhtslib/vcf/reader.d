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
import htslib.vcf;
import htslib.hts_log;
import htslib.kstring;

alias BCFReader = VCFReader;

/** Basic support for reading VCF, BCF files */
auto VCFReader(string fn, int extra_threads = -1, UnpackLevel MAX_UNPACK = UnpackLevel.None)
{
    return VCFReaderImpl!(CoordSystem.zbc, false)(fn, extra_threads, MAX_UNPACK);
}

/** Basic support for reading VCF, BCF files via tabix*/
auto VCFReader(CoordSystem cs)(TabixIndexedFile tbxFile, string chrom, Interval!cs coords, int extra_threads = -1, UnpackLevel MAX_UNPACK = UnpackLevel.None)
{
    return VCFReaderImpl!(cs,true)(tbxFile, chrom, coords, extra_threads, MAX_UNPACK);
}

/** Basic support for reading VCF, BCF files */
struct VCFReaderImpl(CoordSystem cs, bool isTabixFile)
{
    // htslib data structures
    VcfFile     fp;    /// rc htsFile wrapper
    VCFHeader   vcfhdr;    /// rc header wrapper
    Bcf1 b;          /// rc bcf1_t wrapper, record for use in iterator, will be recycled
    UnpackLevel MAX_UNPACK;     /// see htslib.vcf


    private bool initialized;
    private int success;

    @disable this();
    static if(isTabixFile){

        TabixIndexedFile tbx; /// For tabix use
        ReturnType!(this.initializeTabixRange) tbxRange; /// For tabix use

        /// TabixIndexedFile and coordinates ctor
        this(TabixIndexedFile tbxFile, string chrom, Interval!cs coords, int extra_threads = -1, UnpackLevel MAX_UNPACK = UnpackLevel.None)
        {
            this.tbx = tbxFile;
            this.tbxRange = initializeTabixRange(chrom, coords);
            this.tbxRange.empty;
            /// read the header from TabixIndexedFile
            bcf_hdr_t * hdrPtr = bcf_hdr_init(toUTFz!(char *)("r"));
            auto err = bcf_hdr_parse(hdrPtr, toUTFz!(char *)(this.tbx.header));
            this.vcfhdr = VCFHeader(hdrPtr);
            
            this.b = Bcf1(bcf_init1());
            this.b.max_unpack = MAX_UNPACK;
            this.MAX_UNPACK = MAX_UNPACK;
            this.empty;
        }

        /// copy the TabixIndexedFile.region range
        auto initializeTabixRange(string chrom, Interval!cs coords)
        {
            return this.tbx.region(chrom, coords);
        }
    }else{
        /// read existing VCF file
        /// MAX_UNPACK: setting alternate value could speed reading
        this(string fn, int extra_threads = -1, UnpackLevel MAX_UNPACK = UnpackLevel.None)
        {
            import htslib.hts : hts_set_threads;

            assert(!isTabixFile);
            if (fn == "") throw new Exception("Empty filename passed to VCFReader constructor");
            this.fp = VcfFile(vcf_open(toStringz(fn), "r"c.ptr));
            if (!this.fp) throw new Exception("Could not hts_open file");

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
            } else if (extra_threads > 0)
            {
                if ((extra_threads + 1) > totalCPUs)
                    hts_log_warning(__FUNCTION__, "More threads requested than CPU cores detected");
                hts_set_threads(this.fp, extra_threads);
            }
            else if (extra_threads == 0)
            {
                hts_log_debug(__FUNCTION__, "Zero extra threads requested");
            }

            this.vcfhdr = VCFHeader(bcf_hdr_read(this.fp));

            this.b = Bcf1(bcf_init1());
            this.b.max_unpack = MAX_UNPACK;
            this.MAX_UNPACK = MAX_UNPACK;
            this.empty;
        }
    }

    /// Explicit postblit to avoid 
    /// https://github.com/blachlylab/dhtslib/issues/122
    this(this)
    {
        this.fp = fp;
        this.vcfhdr = vcfhdr;
        this.b = b;
        this.MAX_UNPACK = MAX_UNPACK;
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
        static if(isTabixFile){
            if(this.tbxRange.empty) return true;
        }
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
            static if(isTabixFile)
                hts_log_error(__FUNCTION__, "*** ERROR vcf_parse < -1");
            else
                hts_log_error(__FUNCTION__, "*** ERROR bcf_read < -1");
            return true;
        }
        
    }
    /// ditto
    void popFront()
    {
        //     int bcf_read(htsFile *fp, const(bcf_hdr_t) *h, bcf1_t *v);
        // documentation claims returns -1 on critical errors, 0 otherwise
        // however it looks like -1 is EOF and -2 is critical errors?
        static if(isTabixFile){
            if (this.initialized) this.tbxRange.popFront;
            auto str = this.tbxRange.front;
            kstring_t kstr;
            kputs(toUTFz!(char *)(str), &kstr);
            this.success = vcf_parse(&kstr, this.vcfhdr.hdr, this.b);
            
        }else
            this.success = bcf_read(this.fp, this.vcfhdr.hdr, this.b);

    }
    /// ditto
    VCFRecord front()
    {
        // note that VCFRecord's constructor will be responsible for
        // * UnpackLeveling and
        // * destroying
        // its copy
        auto copiedBcf_1t = bcf_dup(this.b);
        return VCFRecord(this.vcfhdr, copiedBcf_1t, this.MAX_UNPACK);
    }

    typeof(this) save()
    {
        typeof(this) newRange = this;
        static if(isTabixFile){
            newRange.tbx = tbx; /// For tabix use
            newRange.tbxRange = tbxRange.save; /// For tabix use
        }else{
            newRange.fp = HtsFile(copyHtsFile(this.fp));
        }
        newRange.vcfhdr = vcfhdr;
        newRange.b = Bcf1(bcf_dup(this.b));
        newRange.initialized = initialized;
        newRange.success = success;
        return newRange;
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
    assert(vcf.count == 15);
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
    auto vcf = VCFReader(tbx, reg.contig, reg.interval);
    assert(vcf.count == 3);
    vcf = VCFReader(tbx, reg.contig, reg.interval);
    VCFRecord rec = vcf.front;
    assert(rec.chrom == "1");
    assert(rec.pos == 3000150);
    assert(rec.refAllele == "C");
    assert(rec.altAllelesAsArray == ["T"]);
    assert(rec.allelesAsArray == ["C","T"]);
    assert(approxEqual(rec.qual, 59.2));
    assert(rec.filter == "PASS");
    assert(vcf.vcfhdr.getSamples == ["A", "B"]);

    vcf.popFront;
    rec = vcf.front;

    assert(rec.chrom == "1");
    assert(rec.pos == 3062914);
    assert(rec.refAllele == "GTTT");
    assert(rec.altAllelesAsArray == ["G"]);
    assert(rec.allelesAsArray == ["GTTT","G"]);
    assert(approxEqual(rec.qual,12.9));
    assert(rec.filter == "q10");
    // assert(rec. == ["C","T"]);

    vcf = VCFReader(tbx, reg.contig, reg.interval);
    auto range1 = vcf.save;
    vcf.popFront;
    auto range2 = vcf.save;
    vcf.popFront;
    auto range3 = vcf.save;
    vcf.popFront;
    auto range4 = vcf.save;
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
    auto tbx = TabixIndexedFile(fn);
    auto reg = getIntervalFromString("1:3000151-3062916");

    auto vcf = VCFReader(tbx, reg.contig, reg.interval);
    assert(vcf.count == 3);
    vcf = VCFReader(tbx, reg.contig, reg.interval, UnpackLevel.None);
    VCFRecord rec = vcf.front;

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
    rec = vcf.front;
    assert(rec.getFormats["DP"].to!int.array == [[32], [32]]);
    assert(rec.line.unpacked == UnpackLevel.Format);
    
}