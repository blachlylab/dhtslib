module dhtslib.vcf.reader;

import std.string: fromStringz, toStringz;
import std.utf: toUTFz;
import std.traits : ReturnType;

import dhtslib.vcf;
import dhtslib.tabix;
import dhtslib.coordinates;
import htslib.vcf;
import htslib.hts_log;
import htslib.kstring;

alias BCFReader = VCFReader;

auto VCFReader(string fn, UnpackLevel MAX_UNPACK = UnpackLevel.None)
{
    return VCFReaderImpl!(CoordSystem.zbc, false)(fn, MAX_UNPACK);
}

auto VCFReader(CoordSystem cs)(TabixIndexedFile tbxFile, string chrom, Interval!cs coords, UnpackLevel MAX_UNPACK = UnpackLevel.None)
{
    return VCFReaderImpl!(cs,true)(tbxFile, chrom, coords, MAX_UNPACK);
}

/** Basic support for reading VCF, BCF files */
struct VCFReaderImpl(CoordSystem cs, bool isTabixFile)
{
    // htslib data structures
    vcfFile     *fp;    /// htsFile
    //bcf_hdr_t   *hdr;   /// header
    VCFHeader   *vcfhdr;    /// header wrapper -- no copies
    bcf1_t* b;          /// record for use in iterator, will be recycled
    UnpackLevel MAX_UNPACK;     /// see htslib.vcf

    private static int refct;

    private bool initialized;
    private int success;

    this(this)
    {
        this.refct++;
    }
    @disable this();
    static if(isTabixFile){

        TabixIndexedFile tbx; /// For tabix use
        ReturnType!(this.initializeTabixRange) tbxRange; /// For tabix use

        /// TabixIndexedFile and coordinates ctor
        this(TabixIndexedFile tbxFile, string chrom, Interval!cs coords, UnpackLevel MAX_UNPACK = UnpackLevel.None)
        {
            this.tbx = tbxFile;
            this.tbxRange = initializeTabixRange(chrom, coords);
            this.tbxRange.empty;
            /// read the header from TabixIndexedFile
            bcf_hdr_t * hdrPtr = bcf_hdr_init(toUTFz!(char *)("r"));
            auto err = bcf_hdr_parse(hdrPtr, toUTFz!(char *)(this.tbx.header));
            this.vcfhdr = new VCFHeader(hdrPtr);
            
            this.b = bcf_init1();
            this.b.max_unpack = MAX_UNPACK;
            this.MAX_UNPACK = MAX_UNPACK;
        }

        /// copy the TabixIndexedFile.region range
        auto initializeTabixRange(string chrom, Interval!cs coords)
        {
            return this.tbx.region(chrom, coords);
        }
    }else{
        /// read existing VCF file
        /// MAX_UNPACK: setting alternate value could speed reading
        this(string fn, UnpackLevel MAX_UNPACK = UnpackLevel.None)
        {
            import htslib.hts : hts_set_threads;

            assert(!isTabixFile);
            if (fn == "") throw new Exception("Empty filename passed to VCFReader constructor");
            this.fp = vcf_open(toStringz(fn), "r"c.ptr);
            if (!this.fp) throw new Exception("Could not hts_open file");
            
            hts_set_threads(this.fp, 1);    // extra decoding thread

            this.vcfhdr = new VCFHeader( bcf_hdr_read(this.fp));

            this.b = bcf_init1();
            this.b.max_unpack = MAX_UNPACK;
            this.MAX_UNPACK = MAX_UNPACK;
        }
    }

    /// dtor
    ~this()
    {
        this.refct--;

        // block file close and bcf1_t free() with reference counting
        // to allow VCFReader to implement Range interface
        if(!this.refct) {
            static if(!isTabixFile){
                const ret = vcf_close(this.fp);
                if (ret != 0) hts_log_error(__FUNCTION__, "couldn't close VCF after reading");
            }

            // Deallocate header
            //bcf_hdr_destroy(this.hdr);
            // 2018-09-15: Do not deallocate header; will be free'd by VCFHeader dtor

            bcf_destroy1(this.b);
        }
    }
    invariant
    {
        assert(this.vcfhdr != null);
    }

    VCFHeader* getHeader()
    {
        return this.vcfhdr;
    }

    /** VCF version, e.g. VCFv4.2 */
    @property string vcfVersion() { return cast(string) fromStringz( bcf_hdr_get_version(this.vcfhdr.hdr) ).idup; }

    /++
        bcf_hrec_t *bcf_hdr_get_hrec(const(bcf_hdr_t) *hdr, int type, const(char) *key, const(char) *value,
                                                                                    const(char) *str_class);
    +/
    // TODO: check key for "", fail if empty
    // TODO: check value, str_class for "" and replace with NULL
    // TODO: Memory leak. We are never freeing the bcf_hrec_t, but, it escapes as pointer inside key/value string
    // TODO: handle structured and general lines
    string[string] getTagByKV(string tagType, T)(string key, string value, string str_class)
    if((tagType == "FILTER" || tagType == "INFO" || tagType == "FORMAT" || tagType == "contig") &&
        (isIntegral!T || isSomeString!T))
    {
        // hlt : header line type
        static if (tagType == "FILTER")     const int hlt = BCF_HL_FLT;
        else static if (tagType == "INFO")  const int hlt = BCF_HL_INFO; // @suppress(dscanner.suspicious.label_var_same_name)
        else static if (tagType == "FORMAT") const int hlt= BCF_HL_FMT; // @suppress(dscanner.suspicious.label_var_same_name)
        else static if (tagType == "contig") const int hlt= BCF_HL_CTG; // @suppress(dscanner.suspicious.label_var_same_name)
        else assert(0);

        bcf_hrec_t *hrec = bcf_hdr_get_hrec(this.vcfhdr.hdr, hlt,   toStringz(key),
                                                                    toStringz(value),
                                                                    toStringz(str_class));

        const int nkeys = hrec.nkeys;
        string[string] kv;

        foreach(int i; 0 .. nkeys) {
            string k = cast(string) fromStringz(hrec.keys[i]);
            string v = cast(string) fromStringz(hrec.vals[i]); 
            kv[k] = v;
        }

        return kv;
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
        return VCFRecord(this.vcfhdr, bcf_dup(this.b), this.MAX_UNPACK);
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
    vcf.empty;
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
    // assert(rec. == ["C","T"]);
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
    hts_log_info(__FUNCTION__, "Testing VCFReader");
    hts_log_info(__FUNCTION__, "Loading test file");
    auto fn = buildPath(dirName(dirName(dirName(dirName(__FILE__)))),"htslib","test","tabix","vcf_file.vcf.gz");
    auto tbx = TabixIndexedFile(fn);
    auto reg = getIntervalFromString("1:3000151-3062916");
    auto vcf = VCFReader(tbx, reg.contig, reg.interval);
    assert(vcf.count == 3);
    vcf = VCFReader(tbx, reg.contig, reg.interval);
    vcf.empty;
    VCFRecord rec = vcf.front;
    assert(rec.chrom == "1");
    assert(rec.pos == 3000150);
    assert(rec.refAllele == "C");
    assert(rec.altAllelesAsArray == ["T"]);
    assert(rec.allelesAsArray == ["C","T"]);
    assert(approxEqual(rec.qual, 59.2));
    assert(rec.filter == "PASS");

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
    hts_log_info(__FUNCTION__, "Testing bcf1_t Unpacking");
    hts_log_info(__FUNCTION__, "Loading test file");
    auto fn = buildPath(dirName(dirName(dirName(dirName(__FILE__)))),"htslib","test","tabix","vcf_file.vcf.gz");
    auto tbx = TabixIndexedFile(fn);
    auto reg = getIntervalFromString("1:3000151-3062916");

    auto vcf = VCFReader(tbx, reg.contig, reg.interval);
    assert(vcf.count == 3);
    vcf = VCFReader(tbx, reg.contig, reg.interval, UnpackLevel.None);
    vcf.empty;
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