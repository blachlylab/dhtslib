module dhtslib.vcf.reader;

import std.string: fromStringz, toStringz;

import dhtslib.vcf.header;
import dhtslib.vcf.record;
import htslib.vcf;
import htslib.hts_log;

alias BCFReader = VCFReader;

/** Basic support for reading VCF, BCF files */
struct VCFReader
{
    // htslib data structures
    vcfFile     *fp;    /// htsFile
    //bcf_hdr_t   *hdr;   /// header
    VCFHeader   *vcfhdr;    /// header wrapper -- no copies
    bcf1_t* b;          /// record for use in iterator, will be recycled

    int MAX_UNPACK;     /// see htslib.vcf

    private static int refct;

    this(this)
    {
        this.refct++;
    }
    @disable this();
    /// read existing VCF file
    /// MAX_UNPACK: setting alternate value could speed reading
    this(string fn, int MAX_UNPACK = BCF_UN_ALL)
    {
        import htslib.hts : hts_set_threads;

        if (fn == "") throw new Exception("Empty filename passed to VCFReader constructor");
        this.fp = vcf_open(toStringz(fn), "r"c.ptr);
        if (!this.fp) throw new Exception("Could not hts_open file");
        
        hts_set_threads(this.fp, 1);    // extra decoding thread

        this.vcfhdr = new VCFHeader( bcf_hdr_read(this.fp));

        this.b = bcf_init1();
        this.b.max_unpack = MAX_UNPACK;
        this.MAX_UNPACK = MAX_UNPACK;
    }
    /// dtor
    ~this()
    {
        this.refct--;

        // block file close and bcf1_t free() with reference counting
        // to allow VCFReader to implement Range interface
        if(!this.refct) {
            const ret = vcf_close(this.fp);
            if (ret != 0) hts_log_error(__FUNCTION__, "couldn't close VCF after reading");

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
        //     int bcf_read(htsFile *fp, const(bcf_hdr_t) *h, bcf1_t *v);
        // documentation claims returns -1 on critical errors, 0 otherwise
        // however it looks like -1 is EOF and -2 is critical errors?
        immutable success = bcf_read(this.fp, this.vcfhdr.hdr, this.b);
        if (success >= 0) return false;
        else if (success == -1) {
            // EOF
            // see htslib my comments https://github.com/samtools/htslib/issues/246
            // and commit 9845bc9a947350d0f34e6ce69e79ab81b6339bd2
            return true;
        }
        else {
            hts_log_error(__FUNCTION__, "*** CRITICAL ERROR bcf_read < -1");
            // TODO: check b->errcode
            return true;
        }
    }
    /// ditto
    void popFront()
    {
        // noop? 
        // free this.b ?
        //bam_destroy1(this.b);

        // TODO: clear (not destroy) the bcf1_t for reuse?
    }
    /// ditto
    VCFRecord front()
    {
        // note that VCFRecord's constructor will be responsible for
        // * unpacking and
        // * destroying
        // its copy
        return VCFRecord(this.vcfhdr, bcf_dup(this.b), this.MAX_UNPACK);
    }
}