module dhtslib.vcf.writer;

import std.datetime;
import std.string: fromStringz, toStringz;
import std.traits: isArray, isDynamicArray, isBoolean, isIntegral, isFloatingPoint, isNumeric, isSomeString;
import std.conv: to, ConvException;
import std.format: format;

import dhtslib.vcf;
import htslib.vcf;
import htslib.hts_log;

alias BCFWriter = VCFWriter;

/** Basic support for writing VCF, BCF files */
struct VCFWriter
{
    // htslib data structures
    vcfFile     *fp;    /// htsFile
    //bcf_hdr_t   *hdr;   /// header
    VCFHeader   *vcfhdr;    /// header wrapper -- no copies
    bcf1_t*[]    rows;   /// individual records

    @disable this();
    /// open file or network resources for writing
    /// setup incl header allocation
    /// if mode==w:
    ///     bcf_hdr_init automatically sets version string (##fileformat=VCFv4.2)
    ///     bcf_hdr_init automatically adds the PASS filter
    this(string fn)
    {
        if (fn == "") throw new Exception("Empty filename passed to VCFWriter constructor");
        this.fp = vcf_open(toStringz(fn), toStringz("w"c));
        if (!this.fp) throw new Exception("Could not hts_open file");

        this.vcfhdr = new VCFHeader( bcf_hdr_init("w\0"c.ptr));
        addFiledate();
        bcf_hdr_sync(this.vcfhdr.hdr);

        /+
        hts_log_trace(__FUNCTION__, "Displaying header at construction");
        //     int bcf_hdr_format(const(bcf_hdr_t) *hdr, int is_bcf, kstring_t *str);
        kstring_t ks;
        bcf_hdr_format(this.vcfhdr.hdr, 0, &ks);
        char[] hdr;
        hdr.length = ks.l;
        import core.stdc.string: memcpy;
        memcpy(hdr.ptr, ks.s, ks.l);
        import std.stdio: writeln;
        writeln(hdr);
        +/
    }
    /// setup and copy a header from another BCF/VCF as template
    this(T)(string fn, T other)
    if(is(T == VCFHeader*) || is(T == bcf_hdr_t*))
    {
        if (fn == "") throw new Exception("Empty filename passed to VCFWriter constructor");
        this.fp = vcf_open(toStringz(fn), toStringz("w"c));
        if (!this.fp) throw new Exception("Could not hts_open file");

        static if(is(T == VCFHeader*)) { this.vcfhdr      = new VCFHeader( bcf_hdr_dup(other.hdr) ); }
        else static if(is(T == bcf_hdr_t*)) { this.vcfhdr = new VCFHeader( bcf_hdr_dup(other) ); }
    }
    /// dtor
    ~this()
    {
        const ret = vcf_close(this.fp);
        if (ret != 0) hts_log_error(__FUNCTION__, "couldn't close VCF after writing");

        // Deallocate header
        //bcf_hdr_destroy(this.hdr);
        // 2018-09-15: Do not deallocate header; will be free'd by VCFHeader dtor
    }
    invariant
    {
        assert(this.vcfhdr != null);
    }

    VCFHeader* getHeader()
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
        bcf1_t *line = new bcf1_t;

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


    auto vw = VCFWriter("/dev/null");

    vw.addHeaderLineRaw("##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">");
    vw.addHeaderLineKV("INFO", "<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">");
    // ##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
    vw.addTag!(HeaderRecordType.Info)("AF", HeaderLengths.OnePerAltAllele, HeaderTypes.Integer, "Number of Samples With Data");
    vw.addTag!(HeaderRecordType.Filter)("filt","test");
    vw.addFilterTag("filt2","test2");

    writeln(vw.getHeader.toString);
    auto expected = "##fileformat=VCFv4.2\n" ~
    "##FILTER=<ID=PASS,Description=\"All filters passed\">\n" ~
    "##filedate=20210909\n" ~
    "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">\n" ~
    "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n" ~
    "##INFO=<ID=AF,Number=A,Type=Integer,Description=Number of Samples With Data>\n" ~
    "##FILTER=<ID=filt,Description=\"test\">\n" ~
    "##FILTER=<ID=filt2,Description=\"test2\">\n" ~
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

    assert(vw.getHeader.toString == expected);

    vw.writeHeader();
}