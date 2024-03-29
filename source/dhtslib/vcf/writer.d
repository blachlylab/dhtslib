module dhtslib.vcf.writer;

import std.datetime;
import std.string: fromStringz, toStringz;
import std.traits: isArray, isDynamicArray, isBoolean, isIntegral, isFloatingPoint, isNumeric, isSomeString;
import std.conv: to, ConvException;
import std.format: format;
import std.parallelism : totalCPUs;

import dhtslib.memory;
import dhtslib.vcf;
import htslib.vcf;
import htslib.hts_log;
import htslib.hfile;
import htslib.hts;

alias BCFWriter = VCFWriter;

/// VCF/BCF/Compressed VCF/ Uncompressed BCF on-disk format.
/// `DEDUCE` will attempt to auto-detect from filename or other means
enum VCFWriterTypes
{
    VCF, /// Regular VCF
    BCF, /// Compressed BCF
    UBCF, /// Uncompressed BCF
    CVCF, /// Compressed VCF (vcf.gz)
    DEDUCE /// Determine based on file extension
}

/** Basic support for writing VCF, BCF files */
struct VCFWriter
{
    /// filename; as usable from D
    string filename;

    /// filename \0-terminated C string; reference needed to avoid GC reaping result of toStringz when ctor goes out of scope
    private const(char)* fn;

    // htslib data structures
    VcfFile     fp;    /// rc htsFile wrapper
    VCFHeader   vcfhdr;    /// header wrapper
    Bcf1[]    rows;   /// individual records

    /// hFILE if required
    private hFILE* f;

    @disable this();
    /// open file or network resources for writing
    /// setup incl header allocation
    /// if mode==w:
    ///     bcf_hdr_init automatically sets version string (##fileformat=VCFv4.2)
    ///     bcf_hdr_init automatically adds the PASS filter
    this(T)(T f, VCFWriterTypes t=VCFWriterTypes.DEDUCE)
    if (is(T == string) || is(T == File))
    {
        auto hdr = VCFHeader( bcf_hdr_init("w\0"c.ptr));
        hdr.addFiledate();
        bcf_hdr_sync(hdr.hdr);
        this(f, hdr, t);
    }
    /// setup and copy a header from another BCF/VCF as template
    this(T, H)(T f, H header, VCFWriterTypes t=VCFWriterTypes.DEDUCE, int extra_threads = -1)
    if((is(H == VCFHeader) || is(H == bcf_hdr_t*)) && (is(T == string) || is(T == File)))
    {

        char[] mode;
        if(t == VCFWriterTypes.BCF) mode=['w','b','\0'];
        else if(t == VCFWriterTypes.UBCF) mode=['w','b','0','\0'];
        else if(t == VCFWriterTypes.VCF) mode=['w','\0'];
        else if(t == VCFWriterTypes.CVCF) mode=['w','z','\0'];
        // open file
        static if (is(T == string))
        {
            if(t == VCFWriterTypes.DEDUCE){
                import std.path:extension, stripExtension;
                auto ext=extension(f);
                if(ext==".bcf") mode=['w','b','\0'];
                else if(ext==".vcf") mode=['w','\0'];
                else if(ext==".cram") mode=['w','c','\0'];
                else if(ext==".gz") {
                    auto extRemoved = stripExtension(f);
                    if (extension(extRemoved) == ".vcf") mode=['w','z','\0'];
                    else {
                        hts_log_error(__FUNCTION__,"extension "~extension(extRemoved)~ext~" not valid");
                        throw new Exception("DEDUCE VCFWriterType used with non-valid extension");
                    }
                }
                else {
                    hts_log_error(__FUNCTION__,"extension "~ext~" not valid");
                    throw new Exception("DEDUCE VCFWriterType used with non-valid extension");
                }
            }
            this.filename = f;
            this.fn = toStringz(f);

            if (f == "") throw new Exception("Empty filename passed to VCFWriter constructor");
            this.fp = vcf_open(toStringz(f), mode.ptr);
            if (!this.fp) throw new Exception("Could not hts_open file");
        }
        else static if (is(T == File))
        {
            assert(t!=VCFWriterTypes.DEDUCE);
            this.filename = f.name();
            this.fn = toStringz(f.name);
            this.f = hdopen(f.fileno, cast(immutable(char)*) "w");
            this.fp = hts_hopen(this.f, this.fn, mode.ptr);
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

        static if(is(H == VCFHeader*)) { this.vcfhdr      = VCFHeader( bcf_hdr_dup(header.hdr) ); }
        else static if(is(H == VCFHeader)) { this.vcfhdr      = VCFHeader( bcf_hdr_dup(header.hdr)); }
        else static if(is(H == bcf_hdr_t*)) { this.vcfhdr = VCFHeader( bcf_hdr_dup(header) ); }
        else assert(0);
    }

    /// Explicit postblit to avoid 
    /// https://github.com/blachlylab/dhtslib/issues/122
    this(this)
    {
        this.fp = fp;
        this.vcfhdr = vcfhdr;
        this.rows = rows;
    }

    VCFHeader getHeader()
    {
        return this.vcfhdr;
    }

    /// Add sample to this VCF
    /// * int bcf_hdr_add_sample(bcf_hdr_t *hdr, const(char) *sample);
    deprecated("use VCFHeader methods instead")
    int addSample(string name)
    in { assert(name != ""); }
    do
    {
        return this.vcfhdr.addSample(name);
    }

    deprecated("use VCFHeader methods instead")
    /// Add a new header line
    int addHeaderLineKV(string key, string value)
    {
        return this.vcfhdr.addHeaderLineKV(key, value);
    }

    deprecated("use VCFHeader methods instead")
    /// Add a new header line -- must be formatted ##key=value
    int addHeaderLineRaw(string line)
    {
        return this.vcfhdr.addHeaderLineRaw(line);
    }

    deprecated("use VCFHeader methods instead")
    /// Add a filedate= headerline, which is not called out specifically in  the spec,
    /// but appears in the spec's example files. We could consider allowing a param here.
    int addFiledate()
    {
        return this.vcfhdr.addFiledate;
    }
    
    /** Add INFO (§1.2.2) or FORMAT (§1.2.4) tag

    The INFO tag describes row-specific keys used in the INFO column;
    The FORMAT tag describes sample-specific keys used in the last, and optional, genotype column.

    Template parameter: string; must be INFO or FORMAT

    The first four parameters are required; NUMBER and TYPE have specific allowable values.
    source and version are optional, but recommended (for INFO only).

    *   id:     ID tag
    *   number: NUMBER tag; here a string because it can also take special values {A,R,G,.} (see §1.2.2)
    *   type:   Integer, Float, Flag, Character, and String
    *   description: Text description; will be double quoted
    *   source:      Annotation source  (eg dbSNP)
    *   version:     Annotation version (eg 142)
    */
    deprecated("use VCFHeader methods instead")
    void addTag(HeaderRecordType tagType)( string id,
                                    HeaderLengths number,
                                    HeaderTypes type,
                                    string description,
                                    string source="",
                                    string _version="")
    if(tagType == HeaderRecordType.Info || tagType == HeaderRecordType.Format)
    {
        this.vcfhdr.addHeaderLine!(tagType)(id, number, type, description, source, _version);
    }

    /** Add FILTER tag (§1.2.3) */
    deprecated("use VCFHeader methods instead")
    void addTag(HeaderRecordType tagType)(string id, string description)
    if(tagType == HeaderRecordType.Filter)
    {
        this.vcfhdr.addHeaderLine!(tagType)(id, description);
    }

    /** Add FILTER tag (§1.2.3) */
    deprecated("use VCFHeader methods instead")
    void addFilterTag(string id, string description)
    {
        this.vcfhdr.addFilter(id, description);
    }

    /** Add contig definition (§1.2.7) to header meta-info 
    
        other: "url=...,md5=...,etc."
    */
    deprecated("use VCFHeader methods instead")
    auto addTag(HeaderRecordType tagType)(const(char)[] id, const int length = 0, string other = "")
    if(tagType == HeaderRecordType.Contig)
    {
        return this.vcfhdr.addTag!(tagType)(id, length, other);
    }
    
    /**
        Add a record

        alleles:    comma-separated string, including ref allele
        qual:       a float, but accepts any numeric
    */
    int addRecord(S, N)(S contig, int pos, S id, S alleles, N qual, S[] filters)
    if ( (isSomeString!S || is(S == char*) ) && isNumeric!N)
    {        
        Bcf1 line = Bcf1(bcf_init1());

        line.rid = bcf_hdr_name2id(this.vcfhdr.hdr, toStringz(contig));
        if (line.rid == -1) hts_log_error(__FUNCTION__, "contig not found");

        line.pos = pos;

        bcf_update_id(this.vcfhdr.hdr, line, (id == "" ? null : toStringz(id)) );  // TODO: could support >1 id with array as with the filters

        bcf_update_alleles_str(this.vcfhdr.hdr, line, toStringz(alleles));

        line.qual = qual;

        // Update filter(s); if/else blocks for speed
        if(filters.length == 0)
        {
            int pass = bcf_hdr_id2int(this.vcfhdr.hdr, BCF_DT_ID, toStringz("PASS"c));
            bcf_update_filter(this.vcfhdr.hdr, line, &pass, 1);
        }
        else if(filters.length == 1)
        {
            int fid = bcf_hdr_id2int(this.vcfhdr.hdr, BCF_DT_ID, toStringz(filters[0]));
            if(fid == -1) hts_log_error(__FUNCTION__, format("filter not found (ignoring): ", filters[0]) );
            bcf_update_filter(this.vcfhdr.hdr, line, &fid, 1);
        }
        else    // TODO: factor out the check for -1 into a safe_update_filter or something
        {
            int[] filter_ids;
            foreach(f; filters) {
                const int fid = bcf_hdr_id2int(this.vcfhdr.hdr, BCF_DT_ID, toStringz(f));
                if(fid == -1) hts_log_error(__FUNCTION__, format("filter not found (ignoring): ", f) );
                else filter_ids ~= fid;
            }
            bcf_update_filter(this.vcfhdr.hdr, line, filter_ids.ptr, cast(int)filter_ids.length );
        }

        // Add a record
        /+
        // .. FORMAT
        bcf_hdr_append(hdr, "##FORMAT=<ID=TF,Number=1,Type=Float,Description=\"Test Float\">");
        float[4] test;
        bcf_float_set_missing(test[0]);
        test[1] = 47.11f;
        bcf_float_set_vector_end(test[2]);
        writeln("pre update format float");
        bcf_update_format_float(this.vcfhdr.hdr, line, toStringz("TF"), &test[0], 4);
        +/
        int tmpi = 1;
        bcf_update_info_int32(this.vcfhdr.hdr, line, toStringz("NS"c), &tmpi, 1);

        // Add the actual sample
        int[4] dp = [ 9000, 1, 2, 3];
        bcf_update_format(this.vcfhdr.hdr, line, toStringz("DP"c), &dp, 1, BCF_HT_INT);
        //int dp = 9000;
        //bcf_update_format_int32(this.vcfhdr.hdr, line, toStringz("DP"c), &dp, 1);
        //auto f = new float;
        //*f = 1.0;
        //bcf_update_format_float(this.vcfhdr.hdr, line, toStringz("XF"c), f, 1);

        this.rows ~= line;

        return 0;
    }

    /// as expected
    int writeHeader()
    {
        return bcf_hdr_write(this.fp, this.vcfhdr.hdr);
    }
    /// as expected
    int writeRecord(ref VCFRecord r)
    {
        const ret = bcf_write(this.fp, this.vcfhdr.hdr, r.line);
        if (ret != 0) hts_log_error(__FUNCTION__, "bcf_write error");
        return ret;
    }
    /// as expected
    int writeRecord(bcf_hdr_t *hdr, bcf1_t *rec)
    {
        hts_log_warning(__FUNCTION__, "pre call");
        const ret = bcf_write(this.fp, hdr, rec);
        hts_log_warning(__FUNCTION__, "post call");
        if (ret != 0) hts_log_error(__FUNCTION__, "bcf_write error");
        return ret;
    }
}

///
debug(dhtslib_unittest)
unittest
{
    import std.exception: assertThrown;
    import std.stdio: writeln, writefln;
    import dhtslib.vcf.writer;

    hts_set_log_level(htsLogLevel.HTS_LOG_TRACE);


    auto vw = VCFWriter("/dev/null", VCFWriterTypes.VCF);

    vw.addHeaderLineRaw("##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">");
    vw.addHeaderLineKV("INFO", "<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">");
    // ##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
    vw.addTag!(HeaderRecordType.Info)("AF", HeaderLengths.OnePerAltAllele, HeaderTypes.Integer, "Number of Samples With Data");
    vw.addTag!(HeaderRecordType.Filter)("filt","test");
    vw.addFilterTag("filt2","test2");

    assert(vw.getHeader.getHeaderRecord(HeaderRecordType.Filter, "filt2").getDescription == "\"test2\"");
    vw.writeHeader();
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
    hts_log_info(__FUNCTION__, "Testing VCFWriterTypes DEDUCE");
    hts_log_info(__FUNCTION__, "Loading test file");
    {
        auto vcf = VCFReader(buildPath(dirName(dirName(dirName(dirName(__FILE__)))),"htslib","test","tabix","vcf_file.vcf"));
        auto vcfw = VCFWriter("/tmp/test_vcf.vcf", vcf.vcfhdr);
        
        vcfw.writeHeader;
        foreach(rec;vcf) {
            vcfw.writeRecord(rec);
        }
        destroy(vcfw);
    }
    {
        auto vcf = VCFReader("/tmp/test_vcf.vcf");
        assert(vcf.count == 15);
        vcf = VCFReader("/tmp/test_vcf.vcf");

        VCFRecord rec = vcf.front;
        assert(rec.chrom == "1");
        assert(rec.pos == 3000149);
        assert(rec.refAllele == "C");
        assert(rec.altAllelesAsArray == ["T"]);
        assert(rec.allelesAsArray == ["C","T"]);
        assert(approxEqual(rec.qual,59.2));
        assert(rec.filter == "PASS");
    }
    {
        auto vcf = VCFReader(buildPath(dirName(dirName(dirName(dirName(__FILE__)))),"htslib","test","tabix","vcf_file.vcf"));
        auto vcfw = VCFWriter("/tmp/test_vcf.bcf", vcf.vcfhdr);
        
        vcfw.writeHeader;
        foreach(rec;vcf) {
            vcfw.writeRecord(rec);
        }
        destroy(vcfw);
    }
    {
        auto vcf = VCFReader("/tmp/test_vcf.bcf");
        assert(vcf.count == 15);
        vcf = VCFReader("/tmp/test_vcf.bcf");
        
        VCFRecord rec = vcf.front;
        assert(rec.chrom == "1");
        assert(rec.pos == 3000149);
        assert(rec.refAllele == "C");
        assert(rec.altAllelesAsArray == ["T"]);
        assert(rec.allelesAsArray == ["C","T"]);
        assert(approxEqual(rec.qual,59.2));
        assert(rec.filter == "PASS");
    }
    {
        auto vcf = VCFReader(buildPath(dirName(dirName(dirName(dirName(__FILE__)))),"htslib","test","tabix","vcf_file.vcf"));
        auto vcfw = VCFWriter("/tmp/test_vcf.vcf.gz", vcf.vcfhdr);
        
        vcfw.writeHeader;
        foreach(rec;vcf) {
            vcfw.writeRecord(rec);
        }
        destroy(vcfw);
    }
    {
        auto vcf = VCFReader("/tmp/test_vcf.vcf.gz");
        assert(vcf.count == 15);
        vcf = VCFReader("/tmp/test_vcf.vcf.gz");
        
        VCFRecord rec = vcf.front;
        assert(rec.chrom == "1");
        assert(rec.pos == 3000149);
        assert(rec.refAllele == "C");
        assert(rec.altAllelesAsArray == ["T"]);
        assert(rec.allelesAsArray == ["C","T"]);
        assert(approxEqual(rec.qual,59.2));
        assert(rec.filter == "PASS");
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
    hts_log_info(__FUNCTION__, "Testing VCFWriterTypes");
    hts_log_info(__FUNCTION__, "Loading test file");
    {
        auto vcf = VCFReader(buildPath(dirName(dirName(dirName(dirName(__FILE__)))),"htslib","test","tabix","vcf_file.vcf"));
        auto vcfw = VCFWriter("/tmp/test_vcf.cvcf", vcf.vcfhdr, VCFWriterTypes.CVCF);
        
        vcfw.writeHeader;
        foreach(rec;vcf) {
            vcfw.writeRecord(rec);
        }
        destroy(vcfw);
    }
    {
        auto vcf = VCFReader("/tmp/test_vcf.cvcf");
        assert(vcf.count == 15);
        vcf = VCFReader("/tmp/test_vcf.cvcf");

        VCFRecord rec = vcf.front;
        assert(rec.chrom == "1");
        assert(rec.pos == 3000149);
        assert(rec.refAllele == "C");
        assert(rec.altAllelesAsArray == ["T"]);
        assert(rec.allelesAsArray == ["C","T"]);
        assert(approxEqual(rec.qual,59.2));
        assert(rec.filter == "PASS");
    }
    {
        auto vcf = VCFReader(buildPath(dirName(dirName(dirName(dirName(__FILE__)))),"htslib","test","tabix","vcf_file.vcf"));
        auto vcfw = VCFWriter("/tmp/test_vcf.ubcf", vcf.vcfhdr, VCFWriterTypes.UBCF);
        
        vcfw.writeHeader;
        foreach(rec;vcf) {
            vcfw.writeRecord(rec);
        }
        destroy(vcfw);
    }
    {
        auto vcf = VCFReader("/tmp/test_vcf.ubcf");
        assert(vcf.count == 15);
        vcf = VCFReader("/tmp/test_vcf.ubcf");
        
        VCFRecord rec = vcf.front;
        assert(rec.chrom == "1");
        assert(rec.pos == 3000149);
        assert(rec.refAllele == "C");
        assert(rec.altAllelesAsArray == ["T"]);
        assert(rec.allelesAsArray == ["C","T"]);
        assert(approxEqual(rec.qual,59.2));
        assert(rec.filter == "PASS");
    }
    {
        auto vcf = VCFReader(buildPath(dirName(dirName(dirName(dirName(__FILE__)))),"htslib","test","tabix","vcf_file.vcf"));
        auto vcfw = VCFWriter("/tmp/test_vcf.txt", vcf.vcfhdr, VCFWriterTypes.VCF);
        
        vcfw.writeHeader;
        foreach(rec;vcf) {
            vcfw.writeRecord(rec);
        }
        destroy(vcfw);
    }
    {
        auto vcf = VCFReader("/tmp/test_vcf.txt");
        assert(vcf.count == 15);
        vcf = VCFReader("/tmp/test_vcf.txt");
        
        VCFRecord rec = vcf.front;
        assert(rec.chrom == "1");
        assert(rec.pos == 3000149);
        assert(rec.refAllele == "C");
        assert(rec.altAllelesAsArray == ["T"]);
        assert(rec.allelesAsArray == ["C","T"]);
        assert(approxEqual(rec.qual,59.2));
        assert(rec.filter == "PASS");
    }
}
