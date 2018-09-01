module dhtslib.sam;

import core.stdc.stdlib: malloc, free;
import std.format;
import std.parallelism: totalCPUs;
import std.stdio: writeln, writefln;
import std.string: fromStringz, toStringz;

import dhtslib.htslib.hts: htsFile, hts_open, hts_close;
import dhtslib.htslib.hts: hts_itr_t;
import dhtslib.htslib.hts: seq_nt16_str;
import dhtslib.htslib.sam;

/**
Encapsulates a SAM/BAM record,
using the bam1_t type for memory efficiency,
and the htslib helper functions for speed.
**/
class Record {
    ///
    bam1_t *b;

    private char[] s;
    private char[] q;

    ///
    this()
    {
        debug(dhtslib_debug) writeln("Record ctor");
        this.b = bam_init1();
    }
    ~this()
    {
        debug(dhtslib_debug) writeln("Record dtor");
        bam_destroy1(this.b);
    }

    /// bool bam_is_rev(bam1_t *b) { return ( ((*b).core.flag & BAM_FREVERSE) != 0 ); }
    @property bool isReversed() { return bam_is_rev(this.b); }
    /// bool bam_is_mrev(bam1_t *b) { return( ((*b).core.flag & BAM_FMREVERSE) != 0); }
    @property bool mateReversed() { return bam_is_mrev(this.b); }
    /// auto bam_get_qname(bam1_t *b) { return (cast(char*)(*b).data); }
    @property string queryName() { return fromStringz(bam_get_qname(this.b)).idup; }
    /// query (and quality string) length
    @property int length() { return this.b.core.l_qseq; }
    ///
    @property char[] sequence()
    {
        // auto bam_get_seq(bam1_t *b) { return ((*b).data + ((*b).core.n_cigar<<2) + (*b).core.l_qname); }
        auto seqdata = bam_get_seq(this.b);

        this.s.length = this.length;
/+        foreach(i; 0 .. this.length) {
            this.s[i] = seq_nt16_str[ bam_seqi(seqdata, i) ];
        }+/
        for(int i; i < this.b.core.l_qseq; i++)
            this.s[i] = seq_nt16_str[ bam_seqi(seqdata, i) ];

        return this.s;
    }
    ///
    @property char[] qscores()
    {
        // auto bam_get_qual(bam1_t *b) { return (*b).data + ((*b).core.n_cigar<<2) + (*b).core.l_qname + (((*b).core.l_qseq + 1)>>1); }
        char * qualdata = cast(char *) bam_get_qual(this.b);

        this.q.length = this.length;
        foreach(i; 0 .. this.length) {
            this.q[i] = cast(char) (qualdata[i] + 33);
        }
        return this.q;
    }
}

/**
Encapsulates a SAM/BAM file.
Implements InputRange interface using htslib calls to ().
*/
struct SAMFile {

    /// filename; reference needed to avoid GC reaping result of toStringz when ctor goes out of scope
    private immutable(char)* fn;

    /// htsFile
    private htsFile *fp;

    /// header struct
    bam_hdr_t *header = null;

    /// SAM/BAM index 
    private hts_idx_t* idx;

    /// htslib data structure representing the BGZF compressed file/stream fp
    //private BGZF* bgzf;

    private kstring_t line;

    /// disallow copying
    @disable this(this);

    ///
    this(string fn)
    {
        debug(dhtslib_debug) { writeln("SAMFile ctor"); }

        // open file
        this.fn = toStringz(fn);
        this.fp = hts_open(this.fn, cast(immutable(char)*)"r");

        // read header
        this.header = sam_hdr_read(this.fp);

    }
    ~this()
    {
        debug(dhtslib_debug) { writeln("SAMFile dtor" ); }

        bam_hdr_destroy(this.header);
        const auto ret = hts_close(fp);
        if (!ret) writeln("There was an error closing %s", this.fn);
    }

    /// number of reference sequences; from bam_hdr_t
    @property int n_targets() const { return this.header.n_targets; }

    /// length of specific reference sequence (by number)
    uint target_len(int target) const
    {
        return this.header.target_len[target];
    }

    /// lengths of the reference sequences
    @property uint[] target_lens() const
    {
        return this.header.target_len[0 .. this.n_targets].dup;
    }

    /// names of the reference sequences
    @property string[] target_names() const
    {
        string[] names;
        names.length = this.n_targets;
        foreach(i; 0 .. this.n_targets) {
            names[i] = fromStringz(this.header.target_name[i]).idup;
        }
        return names;
    }

    /// reference contig name to integer id
    /// Calls int bam_name2id(bam_hdr_t *h, const char *_ref);
    int target_id(string name) 
    {
        return bam_name2id(this.header, toStringz(name));
    }

    /** Query a region and return matching alignments as an InputRange */
    /// Query by chr, start, end
    void query(string chrom, int start, int end)
    {
        string q = format("%s:%d-%d", chrom, start, end);
        //return query(q);
    }
    /// Query by string chr:start-end
    void query(string q)
    {

    }
    /// Query by contig id, start, end
    void query(int tid, int start, int end)
    {
        // hts_itr_t *sam_itr_queryi(const hts_idx_t *idx, int tid, int beg, int end);
        //hts_itr_t* qiter = sam_itr_queryi(this.idx, tid, beg, end);

    }

    /// Iterate through all records in the SAM/BAM/CRAM
    struct AllRecordsRange
    {
        private htsFile     *fp;        // belongs to parent; shared
        private bam_hdr_t   *header;    // belongs to parent; shared
        private bam1_t      *b;

        /// InputRange
        @property bool empty()
        {
            //    int sam_read1(samFile *fp, bam_hdr_t *h, bam1_t *b);
            immutable success = sam_read1(this.fp, this.header, this.b);
            if (success >= 0) return false;
            else if (success == -1) return true;
            else {
                writeln("*** ERROR in sam::SAMFile::AllRecordsRange:empty");
                return true;
            }
        }
        /// ditto
        void popFront()
        {
            // noop? 
            // free this.b ?
            bam_destroy1(this.b);
        }
        /// ditto

    }
    struct RecordRange
    {
        private kstring_t line;     // shut up the compiler
        private htsFile *fp;        // this is shared/will point to a copy possibly used by other iterators
        private hts_itr_t *iter;

        private bam1_t *b;          /// This is the alignment object

        /+this()
        {
            debug { writeln("[sam:RecordRange:ctor]"); }
        }+/
        ~this()
        {
            debug(dhtslib_debug) { writeln("[sam:RecordRange:dtor]"); }
            sam_itr_destroy(this.iter);
        }
        /// InputRange interface; returned by query()
        @property bool empty()
        {
            // equivalent to htslib ks_release
            this.line.l = 0;
            this.line.m = 0;
            this.line.s = null;
            
            // int bgzf_getline(BGZF *fp, int delim, kstring_t *str);
            //immutable int res = bgzf_getline(this.bgzf, cast(int)'\n', &this.line);
            //return (res < 0 ? true : false);
            return true;
        }
        /// ditto
        void popFront()
        {

            free(this.line.s);

            // equivalent to htslib ks_release
            this.line.l = 0;
            this.line.m = 0;
            this.line.s = null;
            
        }
        /// ditto
        string front()
        {
            auto ret = fromStringz(this.line.s).idup;
            return ret;
        }
    }

}