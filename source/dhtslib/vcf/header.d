module dhtslib.vcf.header;

import std.datetime;
import std.string: fromStringz, toStringz;
import std.format: format;
import std.traits : isArray, isIntegral, isSomeString;
import std.conv: to, ConvException;
import std.algorithm : map;
import std.array : array;
import std.utf : toUTFz;

import htslib.vcf;
import htslib.hts_log;

/// Replacement for htslib BCF_HL_*
enum HDR_LINE
{
    FILTER =    0, /// header line: FILTER
    INFO =      1, /// header line: INFO
    FORMAT =    2, /// header line: FORMAT
    CONTIG =    3, /// header line: contig
    STRUCT =    4, /// header line: structured header line TAG=<A=..,B=..>
    GENERIC =   5, /// header line: generic header line
}

/// Replacement for htslib BCF_HT_*
enum HDR_TYPE
{
    FLAG =  0, /// header type: FLAG// header type
    INT =   1, /// header type: INTEGER
    REAL =  2, /// header type: REAL
    STR =   3, /// header type: STRING
    LONG =  BCF_HT_INT | 0x100, // BCF_HT_INT, but for int64_t values; VCF only!
}

/// Replacement for htslib BCF_VL_*
enum HDR_LENGTH
{
    FIXED = 0, /// variable length: fixed (?)// variable length
    VAR =   1, /// variable length: variable
    A =     2, /// variable length: ?
    G =     3, /// variable length: ?
    R =     4, /// variable length: ?
}

/// Replacement for htslib BCF_DT_*
enum HDR_DICT_TYPE
{
    BCF_DT_ID =     0, /// dictionary type: ID
    BCF_DT_CTG =    1, /// dictionary type: CONTIG
    BCF_DT_SAMPLE = 2, /// dictionary type: SAMPLE
}

struct HeaderRecord
{
    HDR_LINE type;
    string key;
    string value;
    int nkeys;
    string[] keys;
    string[] vals;

    /// ctor from a bcf_hrec_t
    this(bcf_hrec_t * rec){
        this.type = cast(HDR_LINE) rec.type;
        this.key = fromStringz(rec.key).dup;
        this.value = fromStringz(rec.value).dup;
        this.nkeys = rec.nkeys;

        for(auto i=0; i < rec.nkeys; i++){
            keys ~= fromStringz(rec.keys[i]).dup;
            vals ~= fromStringz(rec.vals[i]).dup;
        }
    }

    /// get a value from the KV pairs
    /// if key isn't present thows exception
    ref auto opIndex(string index)
    {
        foreach (i, string key; keys)
        {
            if(key == index){
                return vals[i];
            }
        }
        throw new Exception("Key " ~ index ~" not found");
    }

    /// set a value from the KV pairs
    /// if key isn't present a new KV pair is 
    /// added
    void opIndexAssign(string value, string index)
    {
        foreach (i, string key; keys)
        {
            if(key == index){
                vals[i] = value;
                return;
            }
        }
        this.nkeys++;
        keys~=key;
        keys~=value;
    }

    bcf_hrec_t * convert()
    {
        bcf_hrec_t rec;
        rec.type = this.type;
        rec.key = toUTFz!(char *)(this.key);
        rec.value = toUTFz!(char *)(this.value);
        rec.nkeys = this.nkeys;
        rec.keys = this.keys.map!(x=> toUTFz!(char *)(x)).array.ptr;
        rec.vals = this.vals.map!(x=> toUTFz!(char *)(x)).array.ptr;
        return bcf_hrec_dup(&rec);
    }
}

/** Wrapper around bcf_hdr_t

    Q: Why do we have VCFHeader, but not SAMHeader?
    A: Most (all)? of the SAM (bam1_t) manipulation functions do not require a ptr to the header,
    whereas almost all of the VCF (bcf1_t) manipulation functions do. Therefore, we track bcf_hdr_t*
    inside each VCFRecord; this is wrapped by VCFHeader for future convenience (for example,
    now we have @property nsamples; may move the tag reader and writer functions here?)

    In order to avoid double free()'ing an instance bcf_hdr_t,
    this wrapper will be the authoritative holder of of bcf_hdr_t ptrs,
    it shall be passed around by reference, and copies are disabled.
*/
struct VCFHeader
{
    /// Pointer to htslib BCF/VCF header struct; will be freed from VCFHeader dtor 
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

    /// copy this header
    auto dup(){
        return VCFHeader(bcf_hdr_dup(this.hdr));
    }

    /// List of contigs in the header
    @property string[] sequences()
    {
        import core.stdc.stdlib : free;
        int nseqs;

        /** Creates a list of sequence names. It is up to the caller to free the list (but not the sequence names) */
        //const(char) **bcf_hdr_seqnames(const(bcf_hdr_t) *h, int *nseqs);
        const(char*)*ary = bcf_hdr_seqnames(this.hdr, &nseqs);
        if (!nseqs) return [];

        string[] ret;
        ret.reserve(nseqs);

        for(int i; i < nseqs; i++) {
            ret ~= fromStringz(ary[i]).idup;
        }

        free(cast(void*)ary);
        return ret;        
    }

    /// Number of samples in the header
    pragma(inline, true)
    @property int nsamples() { return bcf_hdr_nsamples(this.hdr); }

    // TODO
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
    in { assert(name != ""); }
    do
    {
        assert(this.hdr != null);

        bcf_hdr_add_sample(this.hdr, toStringz(name));

        // AARRRRGGGHHHH
        // https://github.com/samtools/htslib/issues/767
        bcf_hdr_sync(this.hdr);

        return 0;
    }

    /// Add a new header line
    int addHeaderLineKV(string key, string value)
    {
        // TODO check that key is not INFO, FILTER, FORMAT (or contig?)
        string line = format("##%s=%s", key, value);

        return bcf_hdr_append(this.hdr, toStringz(line));
    }
    /// Add a new header line -- must be formatted ##key=value
    int addHeaderLineRaw(string line)
    {
        assert(this.hdr != null);
        //    int bcf_hdr_append(bcf_hdr_t *h, const(char) *line);
        const auto ret = bcf_hdr_append(this.hdr, toStringz(line));
        bcf_hdr_sync(this.hdr);
        return ret;
    }

    /// Add a new header line -- must be formatted ##key=value
    int addHeaderRecord(HeaderRecord rec)
    {
        assert(this.hdr != null);
        auto ret = bcf_hdr_add_hrec(this.hdr, rec.convert);
        bcf_hdr_sync(this.hdr);
        return ret;
    }

    /// Remove all header lines of a particular type
    void removeHeaderLines(HDR_LINE linetype)
    {
        bcf_hdr_remove(this.hdr, linetype, null);
        bcf_hdr_sync(this.hdr);
    }

    /// Remove a header line of a particular type with the key
    void removeHeaderLines(HDR_LINE linetype, string key)
    {
        bcf_hdr_remove(this.hdr, linetype, toStringz(key));
        bcf_hdr_sync(this.hdr);
    }

    HeaderRecord getHeaderRecord(HDR_LINE linetype, string key, string value)
    {
        auto rec = bcf_hdr_get_hrec(this.hdr, linetype, toUTFz!(const(char) *)(key),toUTFz!(const(char) *)(value), null);
        if(!rec) throw new Exception("Record could not be found");
        auto ret = HeaderRecord(rec);
        // bcf_hrec_destroy(rec);
        return ret;
    }

    /// Add a filedate= headerline, which is not called out specifically in  the spec,
    /// but appears in the spec's example files. We could consider allowing a param here.
    int addFiledate()
    {
        return addHeaderLineKV("filedate", (cast(Date) Clock.currTime()).toISOString );
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
    void addTag(string tagType, T)( string id,
                                    T number,
                                    string type,
                                    string description,
                                    string source="",
                                    string _version="")
    if((tagType == "INFO" || tagType == "FORMAT") && (isIntegral!T || isSomeString!T))
    {
        string line;    //  we'll suffix with \0, don't worry

        // check ID
        if (id == "") {
            hts_log_error(__FUNCTION__, "no ID");
            return;
        }

        // check Number is either a special code {A,R,G,.} or an integer
        static if(isSomeString!T) {
        if (number != "A" &&
            number != "R" &&
            number != "G" &&
            number != ".") {
                // not a special ; check if integer
                try {
                    number.to!int;  // don't need to store result, will use format/%s
                }
                catch (ConvException e) {
                    hts_log_error(__FUNCTION__, "Number not A/R/G/. nor an integer");
                    return;
                }
        }
        }

        // check Type
        if (type != "Integer" &&
            type != "Float" &&
            type != "Flag" &&
            type != "Character" &&
            type != "String") {
                hts_log_error(__FUNCTION__, "unrecognized type");
                return;
        }

        // check Description
        if (description == "") hts_log_error(__FUNCTION__, "no description");

        // check Source and Version
        if (source == "" && _version != "") hts_log_error(__FUNCTION__, "version wo source");

        // Format params
        if (source != "" && _version != "")
            line = format("##%s=<ID=%s,Number=%s,Type=%s,Description=\"%s\",Source=\"%s\",Version=\"%s\">\0",
                            tagType, id, number, type, description, source, _version);
        else if (source != "" && _version == "")
            line = format("##%s=<ID=%s,Number=%s,Type=%s,Description=\"%s\",Source=\"%s\">\0",
                            tagType, id, number, type, description, source);
        else
            line = format("##%s=<ID=%s,Number=%s,Type=%s,Description=\"%s\">\0",
                            tagType, id, number, type, description);

        bcf_hdr_append(this.hdr, line.ptr);
    }

    /** Add FILTER tag (§1.2.3) */
    void addTag(string tagType)(string id, string description)
    if(tagType == "FILTER")
    {
        // check ID
        if (id == "") {
            hts_log_error(__FUNCTION__, "no ID");
            return;
        }
        // check Description
        if (description == "") hts_log_error(__FUNCTION__, "no description");

        string line = format("##FILTER=<ID=%s,Description=\"%s\">\0", id, description);
        bcf_hdr_append(this.hdr, line.ptr);
    }

    /** Add FILTER tag (§1.2.3) */
    deprecated void addFilterTag(string id, string description)
    {
        string filter = format("##FILTER=<ID=%s,Description=\"%s\">\0", id, description);
        bcf_hdr_append(this.hdr, filter.ptr);
    }

    /** Add contig definition (§1.2.7) to header meta-info 
    
        other: "url=...,md5=...,etc."
    */
    auto addTag(string tagType)(const(char)[] id, const int length = 0, string other = "")
    if(tagType == "contig" || tagType == "CONTIG")
    {
        const(char)[] contig = "##contig=<ID=" ~ id ~
            (length > 0  ? ",length=" ~ length.to!string : "") ~
            (other != "" ? "," ~ other : "") ~
            ">\0";
        
        return bcf_hdr_append(this.hdr, contig.ptr);
    }
}

///
debug(dhtslib_unittest)
unittest
{
    import std.exception: assertThrown;
    import std.stdio: writeln, writefln;

    hts_set_log_level(htsLogLevel.HTS_LOG_TRACE);


    auto hdr = VCFHeader(bcf_hdr_init("w\0"c.ptr));

    hdr.addHeaderLineRaw("##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">");
    hdr.addHeaderLineKV("INFO", "<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">");
    // ##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
    hdr.addTag!"INFO"("AF", "A", "Integer", "Number of Samples With Data");
    hdr.addHeaderLineRaw("##contig=<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species=\"Homo sapiens\",taxonomy=x>"); // @suppress(dscanner.style.long_line)
    hdr.addHeaderLineRaw("##FILTER=<ID=q10,Description=\"Quality below 10\">");

    // Exercise header
    assert(hdr.nsamples == 0);
    hdr.addSample("NA12878");
    assert(hdr.nsamples == 1);
}

///
debug(dhtslib_unittest)
unittest
{
    import std.exception: assertThrown;
    import std.stdio: writeln, writefln;

    hts_set_log_level(htsLogLevel.HTS_LOG_TRACE);


    auto hdr = VCFHeader(bcf_hdr_init("w\0"c.ptr));

    hdr.addHeaderLineRaw("##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">");
    hdr.addHeaderLineKV("INFO", "<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">");

    auto rec = hdr.getHeaderRecord(HDR_LINE.INFO,"ID","NS");
    assert(rec.type == HDR_LINE.INFO);
    writeln(rec.key);
    writeln(rec.value);
    writeln(rec.vals);
    assert(rec.key == "INFO");
    assert(rec.nkeys == 5);
    assert(rec.keys == ["ID", "Number", "Type", "Description", "IDX"]);
    assert(rec.vals == ["NS", "1", "Integer", "\"Number of Samples With Data\"", "1"]);
    // // ##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
    // hdr.addTag!"INFO"("AF", "A", "Integer", "Number of Samples With Data");
    // hdr.addHeaderLineRaw("##contig=<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species=\"Homo sapiens\",taxonomy=x>"); // @suppress(dscanner.style.long_line)
    // hdr.addHeaderLineRaw("##FILTER=<ID=q10,Description=\"Quality below 10\">");

    // // Exercise header
    // assert(hdr.nsamples == 0);
    // hdr.addSample("NA12878");
    // assert(hdr.nsamples == 1);
}