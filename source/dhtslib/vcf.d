module dhtslib.vcf;

import std.conv: to, convException;
import std.datetime;
import std.format: format;
import std.stdio: writeln;
import std.string: fromStringz, toStringz;
import std.traits: isNumeric, isSomeString;

import dhtslib.htslib.vcf;

alias BCFRecord = VCFRecord;
alias BCFWriter = VCFWriter;

/** Wrapper around bcf_hdr_t */
struct VCFHeader
{
    bcf_hdr_t *hdr;

    /// Number of samples in the header
    pragma(inline, true)
    @property int nsamples() { return bcf_hdr_nsamples(this.hdr); }

}

/** Wrapper around bcf1_t 

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
    /// whether this is a valid record that won't cause bcf_write to segfault
    bool valid;

    bcf1_t* line;   /// htslib structured record

    VCFHeader *vcfheader; /// corresponding header (required)

    /** VCFRecord

    Construct a bcf/vcf record, backed by bcf1_t, from: an existing bcf1_t, parameters, or a VCF line

    */
    this(VCFHeader *h, bcf1_t *b)
    {
        this.vcfheader = h;

        this.line = b;
    }
    /// ditto
    this(bcf_hdr_t *h, bcf1_t *b)
    {
        this.vcfheader = new VCFHeader;
        this.vcfheader.hdr = h;

        this.line = b;
    }
    /// ditto
    this(SS)(string chrom, int pos, string id, string _ref, string alt, float qual, SS filter)
    if (isSomeString!SS || is(SS == string[]))
    {
        this.chrom = chrom;
        this.pos = pos;
        this.id = id;
        // alleles
        this.qual = qual;
        this.filter = filter;
    }
    /// ditto
    /// From VCF line
    this(string line)
    {
        assert(0);
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
        if (rid == -1) writeln("*** ERROR: contig not found");
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
        
        return "";
    }
    /// Remove all entries in FILTER
    void removeFilters()
    {
        auto ret = bcf_update_filter(this.vcfheader.hdr, this.line, null, 0);
        if (!ret) writeln("*** ERROR removing filters in removeFilters");
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
            if (fid == -1) writeln("*** ERROR, filter not in header: ", f);
            fids ~= fid;
        }
        auto ret = bcf_update_filter(this.vcfheader.hdr, this.line, fids.ptr, cast(int)fids.length);
        if (!ret) writeln("*** ERROR setting filters in @property filter, fids=", fids);
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

    /* FORMAT (sample info) */
}


/** Basic support for writing VCF, BCF files */
struct VCFWriter
{
    // htslib data structures
    vcfFile     *fp;    /// htsFile
    bcf_hdr_t   *hdr;   /// header
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

        this.hdr = bcf_hdr_init(toStringz("w"c));
        addHeaderLineKV("filedate", (cast(Date) Clock.currTime()).toISOString );

        bcf_hdr_append(this.hdr, "##contig=<ID=chr3,length=999999,assembly=hg19>");
        //bcf_hdr_append(this.hdr, "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">");
        // WIll NS be written automatically after bcf_hdr_add_samples are complete?
        bcf_hdr_append(this.hdr, "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">");
        bcf_hdr_append(this.hdr, "##FORMAT=<ID=XF,Number=1,Type=Float,Description=\"Test Float\">");
        bcf_hdr_append(this.hdr, "##FILTER=<ID=triallelic,Description=\"Triallelic site\">");

    }
    /// setup and copy a header from another BCF/VCF as template
    this(T)(string fn, T other)
    if(is(T == VCFHeader*) || is(T == bcf_hdr_t*))
    {
        this.fp = vcf_open(toStringz(fn), toStringz("w"c));
        // TODO: if !fp abort

        static if(is(T == VCFHeader*)) { this.hdr = bcf_hdr_dup(other.hdr); }
        else static if(is(T == bcf_hdr_t*)) { this.hdr = bcf_hdr_dup(other); }
    }

    /// copy header lines from a template without overwiting existing lines
    void copyHeaderLines(bcf_hdr_t *other)
    {
        assert(this.hdr != null);
        assert(0);
        //    bcf_hdr_t *bcf_hdr_merge(bcf_hdr_t *dst, const(bcf_hdr_t) *src);
    }

    /// Add sample to this VCF
    /// * int bcf_hdr_add_sample(bcf_hdr_t *hdr, const(char) *sample);
    int addSample(string name)
    {
        bcf_hdr_add_sample(this.hdr, toStringz(name));

        // AARRRRGGGHHHH
        // https://github.com/samtools/htslib/issues/767
        bcf_hdr_sync(this.hdr);

        return 0;
    }

    /// Add a new header line
    int addHeaderLineKV(string key, string value)
    {
        assert(this.hdr != null);

        // TODO check that key is not INFO, FILTER, FORMAT (or contig?)
        string line = format("##%s=%s", key, value);

        return bcf_hdr_append(this.hdr, toStringz(line));
    }
    /// Add a new header line -- must be formatted ##key=value
    int addHeaderLineRaw(string line)
    {
        assert(this.hdr != null);
        //    int bcf_hdr_append(bcf_hdr_t *h, const(char) *line);
        return bcf_hdr_append(this.hdr, toStringz(line));
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
            writeln("*** ERROR -- no ID");
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
                catch (convException e) {
                    writeln("*** ERROR -- Number not A/R/G/. nor an integer");
                    return;
                }
        }

        // check Type
        if (type != "Integer" &&
            type != "Float" &&
            type != "Flag" &&
            type != "Character" &&
            type != "String") {
                writeln("*** ERROR -- unrecognized type");
                return;
        }

        // check Description
        if (description == "") writeln("*** ERROR -- no description");

        // check Source and Version
        if (source == "" && _version != "") writeln("*** ERROR -- version wo source");

        // Format params
        if (source != "" && _version != "")
            info = format("##INFO=<ID=%s,Number=%s,Type=%s,Description=\"%s\",Source=\"%s\",Version=\"%s\"\0",
                            id, number, type, description, source, _version);
        else if (source != "" && _version == "")
            info = format("##INFO=<ID=%s,Number=%s,Type=%s,Description=\"%s\",Source=\"%s\"\0",
                            id, number, type, description, source);
        else
            info = format("##INFO=<ID=%s,Number=%d,Type=%s,Description=\"%s\">\0",
                            id, number, type, description);

        bcf_hdr_append(this.hdr, info.ptr);
    }

    /// Add contig record
    /// example:    "##contig=<ID=chr3,length=999999,assembly=hg19>"
    void addContig(string id, int length, string assembly, string url)
    {

    }

    /** Add a (structured) record */
    int addRecord(bcf1_t *b)
    {
        this.rows ~= b;
        return 0;
    }
    /** Add a (structured) record */
    int addRecord(VCFRecord r)
    {
        this.rows ~= r.line;
        return 0;
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

        line.rid = bcf_hdr_name2id(this.hdr, toStringz(contig));
        if (line.rid == -1) writeln("*** ERROR: contig not found");

        line.pos = pos;

        bcf_update_id(this.hdr, line, (id == "" ? null : toStringz(id)) );  // TODO: could support >1 id with array as with the filters

        bcf_update_alleles_str(this.hdr, line, toStringz(alleles));

        line.qual = qual;

        // Update filter(s); if/else blocks for speed
        if(filters.length == 0)
        {
            int pass = bcf_hdr_id2int(this.hdr, BCF_DT_ID, toStringz("PASS"c));
            bcf_update_filter(this.hdr, line, &pass, 1);
        }
        else if(filters.length == 1)
        {
            int fid = bcf_hdr_id2int(this.hdr, BCF_DT_ID, toStringz(filters[0]));
            if(fid == -1) writeln("*** ERROR: filter not found (ignoring): ", filters[0]);
            bcf_update_filter(this.hdr, line, &fid, 1);
        }
        else    // TODO: factor out the check for -1 into a safe_update_filter or something
        {
            int[] filter_ids;
            foreach(f; filters) {
                int fid = bcf_hdr_id2int(this.hdr, BCF_DT_ID, toStringz(f));
                if(fid == -1) writeln("*** ERROR: filter not found (ignoring): ", f);
                else filter_ids ~= fid;
            }
            bcf_update_filter(this.hdr, line, filter_ids.ptr, cast(int)filter_ids.length );
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
        bcf_update_format_float(this.hdr, line, toStringz("TF"), &test[0], 4);
        +/
        int tmpi = 1;
        bcf_update_info_int32(this.hdr, line, toStringz("NS"c), &tmpi, 1);

        // Add the actual sample
        int[4] dp = [ 9000, 1, 2, 3];
        bcf_update_format(this.hdr, line, toStringz("DP"c), &dp, 1, BCF_HT_INT);
        //int dp = 9000;
        //bcf_update_format_int32(this.hdr, line, toStringz("DP"c), &dp, 1);
        //auto f = new float;
        //*f = 1.0;
        //bcf_update_format_float(this.hdr, line, toStringz("XF"c), f, 1);

        this.rows ~= line;

        return 0;
    }

    /// as expected
    int writeHeader()
    {
        return bcf_hdr_write(this.fp, this.hdr);
    }
    /// as expected
    int writeRecord(VCFRecord r)
    {
        debug { writeln("hdr, bcf1_t: ", this.hdr, r.line); }
        return bcf_write(this.fp, this.hdr, r.line);
    }
    /// as expected
    int writeFile()
    {
        assert(this.hdr != null);

        int ret;

        ret = bcf_hdr_write(this.fp, this.hdr);

        // for each record in this.records
        foreach(r; this.rows) {
            ret = bcf_write(this.fp, this.hdr, r);
            if (ret != 0) writeln("*** ERROR *** VCFWriter:toVCF");
            // May need bcf_hdr_destroy(&r)
        }



        ret = vcf_close(this.fp);
        if (ret != 0) writeln("*** ERROR *** couldn't close VCF after writing");

        // Deallocate header
        bcf_hdr_destroy(this.hdr);

        return ret;
    }
}