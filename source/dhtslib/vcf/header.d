module dhtslib.vcf.header;

import std.datetime;
import std.string: fromStringz, toStringz;
import std.format: format;
import std.traits : isArray, isIntegral, isSomeString;
import std.conv: to, ConvException;
import std.algorithm : map;
import std.array : array;
import std.utf : toUTFz;

import dhtslib.vcf;
import htslib.vcf;
import htslib.hts_log;


/// Struct for easy setting and getting of hrec values for VCFheader
struct HeaderRecord
{
    /// HeaderRecordType type i.e INFO, contig, FORMAT
    HeaderRecordType recType = HeaderRecordType.NULL;

    /// string of HeaderRecordType type i.e INFO, contig, FORMAT ?
    /// or could be ##source=program
    ///               ======
    string key;

    /// mostly empty except for
    /// this ##source=program
    ///               =======
    string value;

    /// number kv pairs
    int nkeys;

    /// kv pair keys
    string[] keys;

    /// kv pair values
    string[] vals;

    /// HDR IDX value
    int idx = -1;

    /// HDR Length value A, R, G, ., FIXED
    HeaderLengths lenthType = HeaderLengths.NULL;

    /// if HDR Length value is FIXED
    /// this is the number
    int length = -1;

    /// HDR Length value INT, FLOAT, STRING
    HeaderTypes valueType = HeaderTypes.NULL;

    invariant
    {
        assert(this.keys.length == this.vals.length);
    }
    /// ctor from a bcf_hrec_t
    this(bcf_hrec_t * rec){

        /// Set the easy stuff
        this.recType = cast(HeaderRecordType) rec.type;
        this.key = fromStringz(rec.key).dup;
        this.value = fromStringz(rec.value).dup;
        this.nkeys = rec.nkeys;

        /// get the kv pairs
        /// special logic for Number and Type
        for(auto i=0; i < rec.nkeys; i++){
            keys ~= fromStringz(rec.keys[i]).dup;
            vals ~= fromStringz(rec.vals[i]).dup;
            if(keys[i] == "Number")
            {
                switch(vals[i]){
                    case "A":
                        this.lenthType = HeaderLengths.A;
                        break;
                    case "G":
                        this.lenthType = HeaderLengths.G;
                        break;
                    case "R":
                        this.lenthType = HeaderLengths.R;
                        break;
                    case ".":
                        this.lenthType = HeaderLengths.VAR;
                        break;
                    default:
                        this.lenthType = HeaderLengths.FIXED;
                        this.length = vals[i].to!int;
                        break;
                }
            }
            if(keys[i] == "Type")
            {
                switch(vals[i]){
                    case "Flag":
                        this.valueType = HeaderTypes.FLAG;
                        break;
                    case "Integer":
                        this.valueType = HeaderTypes.INT;
                        break;
                    case "Float":
                        this.valueType = HeaderTypes.REAL;
                        break;
                    case "Character":
                        this.valueType = HeaderTypes.CHAR;
                        break;
                    case "String":
                        this.valueType = HeaderTypes.STR;
                        break;
                    default:
                        throw new Exception(vals[i]~" is not a know BCF Header Type");
                }
            }
            if(keys[i] == "IDX"){
                this.nkeys--;
                this.idx = this.vals[$-1].to!int;
                this.keys = this.keys[0..$-1];
                this.vals = this.vals[0..$-1];
            }

        }
    }

    /// set Record Type i.e INFO, FORMAT ...
    void setHeaderRecordType(HeaderRecordType line)
    {
        this.recType = line;
        this.key = HeaderRecordTypeStrings[line];
    }

    /// get Record Type i.e INFO, FORMAT ...
    HeaderRecordType getHeaderRecordType()
    {
        return this.recType;
    }

    /// set Value Type length with integer
    void setLength(T)(T number)
    if(isIntegral!T)
    {
        this.lenthType = HeaderLengths.FIXED;
        this["Number"] = number.to!string;
    }

    /// set Value Type length i.e A, R, G, .
    void setLength(HeaderLengths number)
    {
        this.lenthType = number;
        this["Number"] = HeaderLengthsStrings[number];
    }

    /// get Value Type length
    string getLength()
    {
        return this["Number"];
    }

    /// set Value Type i.e Integer, String, Float
    void setValueType(HeaderTypes type)
    {
        this.valueType = type;
        this["Type"] = HeaderTypesStrings[type];
    }

    /// get Value Type i.e Integer, String, Float
    HeaderTypes getValueType()
    {
        return this.valueType;
    }

    /// set ID field
    void setID(string id)
    {
        this["ID"] = id;
    }

    /// get ID field
    string getID()
    {
        return this["ID"];
    }

    /// set Description field
    void setDescription(string des)
    {
        this["Description"] = des;
    }

    /// get Description field
    string getDescription()
    {
        return this["Description"];
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
        keys~=index;
        vals~=value;
    }

    /// convert to bcf_hrec_t for use with htslib functions
    bcf_hrec_t * convert(bcf_hdr_t * hdr)
    {
        if(this.recType == HeaderRecordType.INFO || this.recType == HeaderRecordType.FORMAT){
            assert(this.valueType != HeaderTypes.NULL);
            assert(this.lenthType != HeaderLengths.NULL);
        }

        auto str = this.toString;
        int parsed;
        auto rec = bcf_hdr_parse_line(hdr, toUTFz!(char *)(str), &parsed);
        rec.type = this.recType;
        return rec;
    }

    /// print a string representation of the header record
    string toString()
    {
        string ret = "##" ~ this.key ~ "=" ~ this.value;
        if(this.nkeys > 0){
            ret ~= "<";
            for(auto i =0; i < this.nkeys - 1; i++)
            {
                ret ~= this.keys[i] ~ "=" ~ this.vals[i] ~ ", ";
            }
            ret ~= this.keys[$-1] ~ "=" ~ this.vals[$-1] ~ ">";    
        }
        
        return ret;
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

    int getSampleId(string sam){
        auto ret = bcf_hdr_id2int(this.hdr, HeaderDictTypes.SAMPLE, toUTFz!(char *)(sam));
        if(ret == -1) hts_log_error(__FUNCTION__, "Couldn't find sample in header: " ~ sam);
        return ret;
    }

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

    /// Add a new header line using HeaderRecord 
    int addHeaderRecord(HeaderRecord rec)
    {
        assert(this.hdr != null);
        auto ret = bcf_hdr_add_hrec(this.hdr, rec.convert(this.hdr));
        if(ret < 0)
            hts_log_error(__FUNCTION__, "Couldn't add HeaderRecord");
        auto notAdded = bcf_hdr_sync(this.hdr);
        if(notAdded != 0)
            hts_log_error(__FUNCTION__, "Couldn't add HeaderRecord");
        return ret;
    }

    /// Remove all header lines of a particular type
    void removeHeaderLines(HeaderRecordType linetype)
    {
        bcf_hdr_remove(this.hdr, linetype, null);
        bcf_hdr_sync(this.hdr);
    }

    /// Remove a header line of a particular type with the key
    void removeHeaderLines(HeaderRecordType linetype, string key)
    {
        bcf_hdr_remove(this.hdr, linetype, toStringz(key));
        bcf_hdr_sync(this.hdr);
    }

    HeaderRecord getHeaderRecord(HeaderRecordType linetype, string key, string value)
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

    void addHeaderLine(HeaderRecordType lineType, T)(string id, T number, HeaderTypes type,
                                    string description="",
                                    string source="",
                                    string _version="")
    if((isIntegral!T || is(T == HeaderLengths)) && lineType != HeaderRecordType.NULL )       
    {
        HeaderRecord rec;
        rec.setHeaderRecordType = lineType;
        rec.setID(id);
        rec.setLength(number);
        rec.setValueType(type);
        static if(lineType == HeaderRecordType.INFO || lineType == HeaderRecordType.FILTER || lineType == HeaderRecordType.FORMAT){
            if(description == ""){
                throw new Exception("description cannot be empty for " ~ HeaderRecordTypeStrings[lineType]);    
            }
        }
        rec.setDescription(description);
        if(source != "")
            rec["source"] = "\"%s\"".format(source);
        if(_version != "")
            rec["version"] = "\"%s\"".format(_version);

        this.addHeaderRecord(rec);
    }

    /** Add FILTER tag (§1.2.3) */
    void addHeaderLine(HeaderRecordType lineType)(string id, string description)
    if(lineType == HeaderRecordType.FILTER)
    {
        HeaderRecord rec;
        rec.setHeaderRecordType = lineType;
        rec.setID(id);
        rec.setDescription(description);

        this.addHeaderRecord(rec);
    }

    /** Add FILTER tag (§1.2.3) */
    deprecated void addFilter(string id, string description)
    {
        addHeaderLine!(HeaderRecordType.FILTER)(id, description);
    }

    /// string representation of header
    string toString(){
        import htslib.kstring;
        kstring_t s;

        const int ret = bcf_hdr_format(this.hdr, 0, &s);
        if (ret)
        {
            hts_log_error(__FUNCTION__,
                format("bcf_hdr_format returned nonzero (%d) (likely EINVAL, invalid bcf_hdr_t struct?)", ret));
            return "[VCFHeader bcf_hdr_format parse_error]";
        }

        return cast(string) s.s[0 .. s.l];
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
    hdr.addHeaderLine!(HeaderRecordType.INFO)("AF", HeaderLengths.A, HeaderTypes.INT, "Number of Samples With Data");
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

    auto rec = hdr.getHeaderRecord(HeaderRecordType.INFO,"ID","NS");
    assert(rec.recType == HeaderRecordType.INFO);
    assert(rec.key == "INFO");
    assert(rec.nkeys == 4);
    assert(rec.keys == ["ID", "Number", "Type", "Description"]);
    assert(rec.vals == ["NS", "1", "Integer", "\"Number of Samples With Data\""]);
    assert(rec["ID"] == "NS");

    assert(rec.idx == 1);

    writeln(rec.toString);


    rec = HeaderRecord(rec.convert(hdr.hdr));

    assert(rec.recType == HeaderRecordType.INFO);
    assert(rec.key == "INFO");
    assert(rec.nkeys == 4);
    assert(rec.keys == ["ID", "Number", "Type", "Description"]);
    assert(rec.vals == ["NS", "1", "Integer", "\"Number of Samples With Data\""]);
    assert(rec["ID"] == "NS");
    // assert(rec["IDX"] == "1");
    // assert(rec.idx == 1);

    rec = hdr.getHeaderRecord(HeaderRecordType.INFO,"ID","NS");

    assert(rec.recType == HeaderRecordType.INFO);
    assert(rec.getLength == "1");
    assert(rec.getValueType == HeaderTypes.INT);
    
    rec.idx = -1;

    rec["ID"] = "NS2";

    hdr.addHeaderRecord(rec);
    auto hdr2 = hdr.dup;
    // writeln(hdr2.toString);

    rec = hdr2.getHeaderRecord(HeaderRecordType.INFO,"ID","NS2");
    assert(rec.recType == HeaderRecordType.INFO);
    assert(rec.key == "INFO");
    assert(rec.nkeys == 4);
    assert(rec.keys == ["ID", "Number", "Type", "Description"]);
    assert(rec.vals == ["NS2", "1", "Integer", "\"Number of Samples With Data\""]);
    assert(rec["ID"] == "NS2");

    assert(rec.idx == 3);

    rec = HeaderRecord.init;
    rec.setHeaderRecordType(HeaderRecordType.GENERIC);
    rec.key = "source";
    rec.value = "hello";
    hdr.addHeaderRecord(rec);

    rec = hdr.getHeaderRecord(HeaderRecordType.GENERIC,"source","hello");
    assert(rec.recType == HeaderRecordType.GENERIC);
    assert(rec.key == "source");
    assert(rec.value == "hello");
    assert(rec.nkeys == 0);

    hdr.addHeaderLine!(HeaderRecordType.FILTER)("nonsense","filter");

    rec = hdr.getHeaderRecord(HeaderRecordType.FILTER,"ID","nonsense");
    assert(rec.recType == HeaderRecordType.FILTER);
    assert(rec.key == "FILTER");
    assert(rec.value == "");
    assert(rec.getID == "nonsense");
    writeln(rec.idx);
    assert(rec.idx == 4);

    hdr.removeHeaderLines(HeaderRecordType.FILTER);

    auto expected = "##fileformat=VCFv4.2\n" ~ 
        "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">\n"~
        "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n"~
        "##INFO=<ID=NS2,Number=1,Type=Integer,Description=\"Number of Samples With Data\">\n"~
        "##source=hello\n"~
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
    assert(hdr.toString == expected);

    rec = rec.init;
    rec.setHeaderRecordType(HeaderRecordType.CONTIG);
    rec.setID("test");
    rec["length"] = "5";

    hdr.addHeaderRecord(rec);

    assert(hdr.sequences == ["test"]);


}