module dhtslib.vcf;

import std.conv: to, ConvException;
import std.datetime;
import std.format: format;
import std.stdio: writeln;
import std.string: fromStringz, toStringz;
import std.traits: isBoolean, isIntegral, isFloatingPoint, isNumeric, isSomeString;

import dhtslib.htslib.hts_log;
import dhtslib.htslib.vcf;

alias BCFRecord = VCFRecord;
alias BCFWriter = VCFWriter;

/** Wrapper around bcf_hdr_t

    In order to avoid double free()'ing an instance bcf_hdr_t,
    this wrapper will be the authoritative holder of of bcf_hdr_t ptrs,
    it shall be passed around by reference, and copies are disabled.
*/
struct VCFHeader
{
    bcf_hdr_t *hdr;

    // Copies have to be disabled to avoid double free()
    @disable this(this);

    ~this()
    {
        // Deallocate header
        if (this.hdr != null) bcf_hdr_destroy(this.hdr);
    }

    invariant
    {
        assert(this.hdr != null);
    }


    /// Number of samples in the header
    pragma(inline, true)
    @property int nsamples() { return bcf_hdr_nsamples(this.hdr); }

}

/** Wrapper around bcf1_t 

    Because it uses bcf1_t internally, it must conform to the BCF2 part
    of the VCFv4.2 specs, rather than the loosey-goosey VCF specs. i.e.,
    INFO, CONTIG, FILTER records must exist in the header.

    TODO: Does this need to be kept in a consistent state?
    Ideally, VCFWriter would reject invalid ones, but we are informed
    that it is invalid (e.g. if contig not found) while building this
    struct; bcf_write1 will actually segfault, unfortunately. I'd like
    to avoid expensive validate() calls for every record before writing
    if possible, which means keeping this consistent. However, not
    sure what to do if error occurs when using the setters herein?
*/
struct VCFRecord
{
    bcf1_t* line;   /// htslib structured record

    VCFHeader *vcfheader;   /// corresponding header (required);
                            /// is ptr to avoid copying struct containing ptr to bcf_hdr_t (leads to double free())

    /** VCFRecord

    Construct a bcf/vcf record, backed by bcf1_t, from: an existing bcf1_t, parameters, or a VCF line.

    Internal backing by bcf1_t means it must conform to the BCF2 rules -- i.e., header must contain
    appropriate INFO, CONTIG, and FILTER lines.

    */
    this(bcf_hdr_t *h, bcf1_t *b)
    {
        this.vcfheader = new VCFHeader(h);  // this looks like it will also lead to a double free() bug if we don't own bcf_hdr_t ...

        this.line = b;
    }
    /// ditto
    this(VCFHeader *h, bcf1_t *b)
    {
        this.vcfheader = h;

        this.line = b;
    }
    /// ditto
    this(SS)(VCFHeader *vcfhdr, string chrom, int pos, string id, string _ref, string alt, float qual, SS filter, )
    if (isSomeString!SS || is(SS == string[]))
    {
        this.line = bcf_init1();
        this.vcfheader = vcfhdr;
        
        this.chrom = chrom;
        this.pos = pos;
        this.updateID(id);

        // alleles
        immutable string alleles = _ref ~ "," ~ alt ~ "\0";
        bcf_update_alleles_str(this.vcfheader.hdr, this.line, alleles.ptr);

        this.qual = qual;
        this.filter = filter;
    }
    /// ditto
    /// From VCF line
    this(string line)
    {
        assert(0);
    }
    /// disable copying to prevent double-free (which should not come up except when writeln'ing)
    @disable this(this);
    /// dtor
    ~this()
    {
        if (this.line) bcf_destroy1(this.line);
    }

    //////// FIXED FIELDS ////////
    
    /* CHROM */
    @property
    const(char)[] chrom()
    {
        return fromStringz(bcf_hdr_id2name(this.vcfheader.hdr, this.line.rid));
    }
    @property
    void chrom(string c)
    {
        auto rid = bcf_hdr_name2id(this.vcfheader.hdr, toStringz(c));
        if (rid == -1) {
            hts_log_error(__FUNCTION__, format("contig not found: %s", c));
            throw new Exception("contig not found");
        }
        else line.rid = rid;
    }

    /* POS */
    @property
    int pos()
    {
        return this.line.pos;
    }
    @property
    void pos(int p)
    {
        this.line.pos = p;
    }

    /* ID */
    /// Sets new ID string; comma-separated list allowed but no dup checking performed
    int updateID(string id)
    {
        if(id == "") return 0;
        return bcf_update_id(this.vcfheader.hdr, this.line, toStringz(id));
    }
    /// Append an ID (htslib performs duplicate checking)
    int addID(string id)
    {
        if(id == "") return 0;
        return bcf_add_id(this.vcfheader.hdr, this.line, toStringz(id));
    }


    /* Alleles (REF, ALT) */


    /* Quality (QUAL) */
    @property float qual()      { return this.line.qual; }
    @property void qual(float q){ this.line.qual = q; }


    /* FILTER */
    /// get FILTER column (TODO -- nothing in htslib)
    @property const(char)[] filter()
    {
        //TODO
        return "";
    }
    /// Remove all entries in FILTER
    void removeFilters()
    {
        auto ret = bcf_update_filter(this.vcfheader.hdr, this.line, null, 0);
        if (!ret) hts_log_error(__FUNCTION__,"error removing filters in removeFilters");
    }
    /// Set the FILTER column to f
    @property void filter(string f)
    {
        this.filter([f]);
    }
    /// Set the FILTER column to f0,f1,f2...
    @property void filter(string[] fs)
    {
        int[] fids;
        foreach(f; fs) {
            int fid = bcf_hdr_id2int(this.vcfheader.hdr, BCF_DT_ID, toStringz(f));
            if (fid == -1) hts_log_warning(__FUNCTION__, format("filter not found in header (ignoring): %s", f) );
            else fids ~= fid;
        }
        bcf_update_filter(this.vcfheader.hdr, this.line, fids.ptr, cast(int)fids.length);
    }

    /// Add a filter; from htslib: 
    /// "If flt_id is PASS, all existing filters are removed first. If other than PASS, existing PASS is removed."
    int addFilter(string f)
    {
        return bcf_add_filter(this.vcfheader.hdr, this.line, 
            bcf_hdr_id2int(this.vcfheader.hdr, BCF_DT_ID, toStringz(f)));
    }
    /+
    /**
     *  bcf_remove_filter() - removes from the FILTER column
     *  @flt_id:   filter ID to remove, numeric ID returned by bcf_hdr_id2int(hdr, BCF_DT_ID, "PASS")
     *  @pass:     when set to 1 and no filters are present, set to PASS
     */
    int bcf_remove_filter(const bcf_hdr_t *hdr, bcf1_t *line, int flt_id, int pass);
    /**
     *  Returns 1 if present, 0 if absent, or -1 if filter does not exist. "PASS" and "." can be used interchangeably.
     */
    int bcf_has_filter(const bcf_hdr_t *hdr, bcf1_t *line, char *filter);
    +/

    /* INFO */
    /+
        int tmpi = 1;
        bcf_update_info_int32(this.vcfhdr.hdr, line, toStringz("NS"c), &tmpi, 1);
    +/
    /// Add a tag:value to the INFO column -- tag must already exist in the header
    void addInfo(T)(string tag, T data)
    {
        int ret = -1;

        static if(isIntegral!T)
            ret = bcf_update_info_int32(this.vcfheader.hdr, this.line, toStringz(tag), &data, 1);
        
        else static if(isFloatingPoint!T) {
            auto flt = cast(float) data;    // simply passing "2.0" (or whatever) => data is a double
            ret = bcf_update_info_float(this.vcfheader.hdr, this.line, toStringz(tag), &flt, 1);
        }

        else static if(isSomeString!T)
            ret = bcf_update_info_string(this.vcfheader.hdr, this.line, toStringz(tag), toStringz(data));
        
        else static if(isBoolean!T) {
            immutable int set = data ? 1 : 0; // if data == true, pass 1 to bcf_update_info_flag(n=); n=0 => clear flag 
            ret = bcf_update_info_flag(this.vcfheader.hdr, this.line, toStringz(tag), null, set);
        }
        
        if (ret == -1) hts_log_warning(__FUNCTION__, format("Couldn't add tag (ignoring): %s", data));
    }
    /// ditto
    /// This handles a vector of values for the tag
    void addInfo(T)(string tag, T[] data)
    if(!is(T==immutable(char)))             // otherwise string::immutable(char)[] will match the other template
    {

    }

    /* FORMAT (sample info) */

    string toString()
    {
        return "[VCFRecord]";
    }
}


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
    /// bcf_hdr_init automatically sets version string (##fileformat=VCFv4.2)
    /// bcf_hdr_init automatically adds the PASS filter
    this(string fn)
    {
        this.fp = vcf_open(toStringz(fn), toStringz("w"c));
        // TODO: if !fp abort

        this.vcfhdr = new VCFHeader( bcf_hdr_init(toStringz("w"c)));
        addHeaderLineKV("filedate", (cast(Date) Clock.currTime()).toISOString );

        bcf_hdr_append(this.vcfhdr.hdr, "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">");
        bcf_hdr_append(this.vcfhdr.hdr, "##FORMAT=<ID=XF,Number=1,Type=Float,Description=\"Test Float\">");
        bcf_hdr_append(this.vcfhdr.hdr, "##FILTER=<ID=triallelic,Description=\"Triallelic site\">");

    }
    /// setup and copy a header from another BCF/VCF as template
    this(T)(string fn, T other)
    if(is(T == VCFHeader*) || is(T == bcf_hdr_t*))
    {
        this.fp = vcf_open(toStringz(fn), toStringz("w"c));
        // TODO: if !fp abort

        static if(is(T == VCFHeader*)) { this.vcfhdr.hdr  = new VCFHeader( bcf_hdr_dup(other.hdr) ); }
        else static if(is(T == bcf_hdr_t*)) { this.vcfhdr = new VCFHeader( bcf_hdr_dup(other) ); }
    }
    /// dtor
    ~this()
    {
        const ret = vcf_close(this.fp);
        if (ret != 0) hts_log_error(__FUNCTION__,"couldn't close VCF after writing");

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

    /// copy header lines from a template without overwiting existing lines
    void copyHeaderLines(bcf_hdr_t *other)
    {
        assert(this.vcfhdr != null);
        assert(0);
        //    bcf_hdr_t *bcf_hdr_merge(bcf_hdr_t *dst, const(bcf_hdr_t) *src);
    }

    /// Add sample to this VCF
    /// * int bcf_hdr_add_sample(bcf_hdr_t *hdr, const(char) *sample);
    int addSample(string name)
    {
        assert(this.vcfhdr != null);

        bcf_hdr_add_sample(this.vcfhdr.hdr, toStringz(name));

        // AARRRRGGGHHHH
        // https://github.com/samtools/htslib/issues/767
        bcf_hdr_sync(this.vcfhdr.hdr);

        return 0;
    }

    /// Add a new header line
    int addHeaderLineKV(string key, string value)
    {
        // TODO check that key is not INFO, FILTER, FORMAT (or contig?)
        string line = format("##%s=%s", key, value);

        return bcf_hdr_append(this.vcfhdr.hdr, toStringz(line));
    }
    /// Add a new header line -- must be formatted ##key=value
    int addHeaderLineRaw(string line)
    {
        assert(this.vcfhdr != null);
        //    int bcf_hdr_append(bcf_hdr_t *h, const(char) *line);
        return bcf_hdr_append(this.vcfhdr.hdr, toStringz(line));
    }
    /** Add INFO tag (§1.2.2)

    The first four parameters are required; NUMBER and TYPE have specific allowable values.
    source and version are optional, but recommended.

    *   id:     ID tag
    *   number: NUMBER tag; here a string because it can also take special values {A,R,G,.} (see §1.2.2)
    *   type:   Integer, Float, Flag, Character, and String
    *   description: Text description; will be double quoted
    *   source:      Annotation source  (eg dbSNP)
    *   version:     Annotation version (eg 142)
    */
    void addInfoTag(string id, string number, string type, string description, string source="", string _version="")
    {
        string info;

        // check ID
        if (id == "") {
            hts_log_error(__FUNCTION__,"no ID");
            return;
        }

        // check Number is either a special code {A,R,G,.} or an integer
        if (number != "A" &&
            number != "R" &&
            number != "G" &&
            number != ".") {
                // not a special ; check if integer
                try {
                    number.to!int;  // don't need to store result, will use format/%s
                }
                catch (ConvException e) {
                    hts_log_error(__FUNCTION__,"Number not A/R/G/. nor an integer");
                    return;
                }
        }

        // check Type
        if (type != "Integer" &&
            type != "Float" &&
            type != "Flag" &&
            type != "Character" &&
            type != "String") {
                hts_log_error(__FUNCTION__,"unrecognized type");
                return;
        }

        // check Description
        if (description == "") hts_log_error(__FUNCTION__,"no description");

        // check Source and Version
        if (source == "" && _version != "") hts_log_error(__FUNCTION__,"version wo source");

        // Format params
        if (source != "" && _version != "")
            info = format("##INFO=<ID=%s,Number=%s,Type=%s,Description=\"%s\",Source=\"%s\",Version=\"%s\">\0",
                            id, number, type, description, source, _version);
        else if (source != "" && _version == "")
            info = format("##INFO=<ID=%s,Number=%s,Type=%s,Description=\"%s\",Source=\"%s\">\0",
                            id, number, type, description, source);
        else
            info = format("##INFO=<ID=%s,Number=%s,Type=%s,Description=\"%s\">\0",
                            id, number, type, description);

        bcf_hdr_append(this.vcfhdr.hdr, info.ptr);
    }

    /** Add FILTER tag (§1.2.3) */
    void addFilterTag(string id, string description)
    {
        string filter = format("##FILTER=<ID=%s,Description=\"%s\">\0", id, description);
        bcf_hdr_append(this.vcfhdr.hdr, filter.ptr);
    }
    alias addFilter = addFilterTag;

    
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
        if (line.rid == -1) hts_log_error(__FUNCTION__,"contig not found");

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
                int fid = bcf_hdr_id2int(this.vcfhdr.hdr, BCF_DT_ID, toStringz(f));
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
        if (ret != 0) hts_log_error(__FUNCTION__,"bcf_write error");
        return ret;
    }
    /// as expected
    int writeRecord(bcf_hdr_t *hdr, bcf1_t *rec)
    {
        hts_log_warning(__FUNCTION__, "pre call");
        const ret = bcf_write(this.fp, hdr, rec);
        hts_log_warning(__FUNCTION__, "post call");
        if (ret != 0) hts_log_error(__FUNCTION__,"bcf_write error");
        return ret;
    }
}

debug
{
    // module constructor
    static this()
    {
        // TRACE is highest log levle (above DEBUG)
        hts_set_log_level(htsLogLevel.HTS_LOG_TRACE);
    }
}