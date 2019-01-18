module dhtslib.sam;

import core.stdc.stdlib: calloc, free;
import std.format;
import std.parallelism: totalCPUs;
import std.stdio: writeln, writefln,stderr;
import std.string: fromStringz, toStringz;

import dhtslib.htslib.hts: htsFile, hts_open, hts_close,hts_itr_next;
import dhtslib.htslib.hts: hts_itr_t;
import dhtslib.htslib.hts: seq_nt16_str;
import dhtslib.htslib.hts: hts_set_threads;

import dhtslib.htslib.hts_log;
import dhtslib.htslib.kstring;
import dhtslib.htslib.sam;

void test_log(string ctx, string msg)
{
    import std.stdio;
    stderr.writeln(ctx, msg);
}

/**
Encapsulates a SAM/BAM/CRAM record,
using the bam1_t type for memory efficiency,
and the htslib helper functions for speed.
**/
class SAMRecord {
    ///
    bam1_t *b;

    ///
    this()
    {
        //debug(dhtslib_debug) hts_log_debug(__FUNCTION__, "ctor()"); /// This line triggers memory error when __FUNCTION__, but not when "Other string"
        test_log(__FUNCTION__, "ctor()");   /// This line will also trigger the memory error when __FUNCTION__, but not other strings
        //writeln(__FUNCTION__);    // this will not trigger the memory error
        this.b = bam_init1();
        //assert(0);                // This will elide(?) the memory error
        //assert(1 == 2);           // This will elide(?) the memory error
    }
    ///
    this(bam1_t *b)
    {
        debug(dhtslib_debug) hts_log_debug(__FUNCTION__, "ctor(bam1_t *b)");
        this.b = b;
    }
    ~this()
    {
        writeln("Record dtor");
        debug(dhtslib_debug) hts_log_debug(__FUNCTION__, "dtor");
        //bam_destroy1(this.b); // we don't own it!
    }

    /// bool bam_is_rev(bam1_t *b) { return ( ((*b).core.flag & BAM_FREVERSE) != 0 ); }
    @property bool isReversed() { return bam_is_rev(this.b); }
    /// bool bam_is_mrev(bam1_t *b) { return( ((*b).core.flag & BAM_FMREVERSE) != 0); }
    @property bool mateReversed() { return bam_is_mrev(this.b); }
    /// auto bam_get_qname(bam1_t *b) { return (cast(char*)(*b).data); }
    @property char[] queryName() { return fromStringz(bam_get_qname(this.b)); }
    /// query (and quality string) length
    @property int length() { return this.b.core.l_qseq; }
    /// see samtools/sam_view.c: get_read
    @property char* sequence()
    {
        // calloc fills with \0; +1 len for Cstring
        char *s = cast(char *) calloc(1, this.b.core.l_qseq + 1);
        //char[] s;
        //s.length = this.b.core.l_qseq;

        // auto bam_get_seq(bam1_t *b) { return ((*b).data + ((*b).core.n_cigar<<2) + (*b).core.l_qname); }
        auto seqdata = bam_get_seq(this.b);

        for(int i; i < this.b.core.l_qseq; i++) {
            if (this.b.core.flag & BAM_FREVERSE)    s[i] = seq_nt16_str[seq_comp_table[bam_seqi(seqdata, i)]];
            else                                    s[i] = seq_nt16_str[bam_seqi(seqdata, i)];
        }
        if (this.b.core.flag & BAM_FREVERSE) reverse(s);
        return s;
    }
    /// see samtools/sam_view.c: get_quality
    @property char* qscores()
    {
        // calloc fills with \0; +1 len for Cstring
        char *q = cast(char *) calloc(1, this.b.core.l_qseq + 1);
        //char[] q;
        //q.length = this.b.core.l_qseq;

        // auto bam_get_qual(bam1_t *b) { return (*b).data + ((*b).core.n_cigar<<2) + (*b).core.l_qname + (((*b).core.l_qseq + 1)>>1); }
        char * qualdata = cast(char *) bam_get_qual(this.b);

        for(int i; i < this.b.core.l_qseq; i++)
            q[i] = cast(char) (qualdata[i] + 33);
        if (this.b.core.flag & BAM_FREVERSE) reverse(q);

        return q;
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
        import std.parallelism: totalCPUs;

        debug(dhtslib_debug) { writeln("SAMFile ctor"); }

        // open file
        this.fn = toStringz(fn);
        this.fp = hts_open(this.fn, cast(immutable(char)*)"r");

        if (totalCPUs > 1) {
            // TODO: make optional
            // TODO hts_log
            hts_log_info(__FUNCTION__, format("%d CPU cores detected; enabling multithreading", totalCPUs));
            //stderr.writefln("dhtslib: %d CPU cores detected; enabling multithreading.", totalCPUs);
            // hts_set_threads adds N _EXTRA_ threads, so totalCPUs - 1 seemed reasonable,
            // but overcomitting by 1 thread (i.e., passing totalCPUs) buys an extra 3% on my 2-core 2013 Mac
            hts_set_threads(this.fp, totalCPUs );
        }

        // read header
        this.header = sam_hdr_read(this.fp);
        this.idx=sam_index_load(this.fp,this.fn);
        debug(dhtslib_debug) { writeln("SAM index",this.idx); }
    }
    ~this()
    {
        debug(dhtslib_debug) { writeln("SAMFile dtor" ); }

        bam_hdr_destroy(this.header);
        //TODO:hts_close segfaults
        //const auto ret = hts_close(fp);
        //if (ret < 0) writefln("There was an error closing %s", fromStringz(this.fn));
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
    auto query(string chrom, int start, int end)
    {
        string q = format("%s:%d-%d", chrom, start, end);
        return query(q);
    }
    /// Query by string chr:start-end
    auto query(string q)
    {
        return RecordRange(this.fp,this.fn,this.idx,this.header,q);
    }
    /// Query by contig id, start, end
    auto query(int tid, int start, int end)
    {
        return RecordRange(this.fp,this.fn,this.idx,tid,start,end);
         // *sam_itr_queryi(const hts_idx_t *idx, int tid, int beg, int end);
        //hts_itr_t* qiter = sam_itr_queryi(this.idx, tid, beg, end);
    }

    AllRecordsRange all_records()
    {
        return AllRecordsRange(this.fp, this.header, bam_init1());
    }

    /// Iterate through all records in the SAM/BAM/CRAM
    struct AllRecordsRange
    {
        private htsFile     *fp;        // belongs to parent; shared
        private bam_hdr_t   *header;    // belongs to parent; shared
        private bam1_t      *b;

        ~this()
        {
            debug(dhtslib_debug) hts_log_debug(__FUNCTION__, "dtor");
        }

        /// InputRange
        @property bool empty()
        {
            hts_log_debug(__FUNCTION__, "Top of function");

            //    int sam_read1(samFile *fp, bam_hdr_t *h, bam1_t *b);
            immutable success = sam_read1(this.fp, this.header, this.b);
            if (success >= 0) return false;
            else if (success == -1) return true;
            else {
                writeln("*** ERROR in sam::SAMFile::AllRecordsRange:empty");
                hts_log_error(__FUNCTION__, "*** ERROR in sam::SAMFile::AllRecordsRange:empty, sam_read1 < -1");
                return true;
            }
        }
        /// ditto
        void popFront()
        {
            // noop? 
            // free this.b ?
            
            //bam_destroy1(this.b);
        }
        /// ditto
        SAMRecord front()
        {
            return new SAMRecord(bam_dup1(this.b));
        }

    }
    struct RecordRange
    {
        htsFile * fp;
        hts_itr_t * itr;
        bam1_t * b;
        int r;
        immutable(char)* q;
        this(htsFile * fp,immutable(char)* fn, hts_idx_t* idx, int tid, int start, int end){
            this.itr=sam_itr_queryi(idx,tid,start,end);
            import dhtslib.htslib.bgzf:bgzf_open;
            ////For whatever reason hts_open doesn't initialize the BGZF pointer
            ////in the the htsFile so we do it ourselves
            //TODO: this shouldn't be required
            fp.fp.bgzf=bgzf_open(fn,cast(immutable(char)*)"r");
            this.fp=fp;
            b=bam_init1();
            debug(dhtslib_debug) { writeln("sam_itr null? ",(cast(int)itr)==0); }
            popFront();
        }
        this(htsFile * fp,immutable(char)* fn, hts_idx_t* idx, bam_hdr_t* header, string query){
            this.itr=sam_itr_querys(idx,header,toStringz(query));
            import dhtslib.htslib.bgzf:bgzf_open;
            ////For whatever reason hts_open doesn't initialize the BGZF pointer
            ////in the the htsFile so we do it ourselves
            //TODO: this shouldn't be required
            fp.fp.bgzf=bgzf_open(fn,cast(immutable(char)*)"r");
            this.fp=fp;
            b=bam_init1();
            debug(dhtslib_debug) { writeln("sam_itr null? ",(cast(int)itr)==0); }
            popFront();
        }
        SAMRecord front(){
            return new SAMRecord(bam_dup1(b));
        }
        void popFront(){
            r=sam_itr_next(fp,itr,b);
        }
        @property bool empty(){
            return (r <=0 && itr.finished) ? true:false;
        }
    }

}

/// Nucleotide complement table; from samtools/sam_view.c
private const(char)[16] seq_comp_table = [ 0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15 ];

/// Reverse a string in place; from samtools/sam_view.c
//  Could be sped up, and potentially made safer? by passing strlen since already known
private char *reverse(char *str)
{
    import core.stdc.string: strlen;

    ulong i = strlen(str)-1,j=0;
    char ch;
    while (i>j) {
        ch = str[i];
        str[i]= str[j];
        str[j] = ch;
        i--;
        j++;
    }
    return str;
}