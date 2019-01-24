module dhtslib.sam;

import core.stdc.stdlib: calloc, free;
import std.format;
import std.parallelism: totalCPUs;
import std.stdio: writeln, writefln,stderr,File;
import std.string: fromStringz, toStringz;

import dhtslib.htslib.hts: htsFile, hts_open, hts_close, hdopen,hclose, hFILE, hts_hopen;
import dhtslib.htslib.hts: hts_itr_t,hts_itr_multi_t,hts_reglist_t,hts_pair32_t;
import dhtslib.htslib.hts: seq_nt16_str;
import dhtslib.htslib.hts: hts_set_threads;

import dhtslib.htslib.hts_log;
import dhtslib.htslib.kstring;
import dhtslib.htslib.sam;
import dhtslib.cigar;
import dhtslib.tagvalue;

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
        //test_log(__FUNCTION__, "ctor()");   /// This line will also trigger the memory error when __FUNCTION__, but not other strings
        //writeln(__FUNCTION__);    // this will not trigger the memory error
        this.b = bam_init1();
        //assert(0);                // This will elide(?) the memory error
        //assert(1 == 2);           // This will elide(?) the memory error
    }
    ///
    this(bam1_t *b)
    {
        //debug(dhtslib_debug) hts_log_debug(__FUNCTION__, "ctor(bam1_t *b)");
        this.b = b;
    }
    ~this()
    {
        //debug(dhtslib_debug) hts_log_debug(__FUNCTION__, "dtor");
        bam_destroy1(this.b); // we created our own in default ctor, or received copy via bam_dup1
    }

    /// bool bam_is_rev(bam1_t *b) { return ( ((*b).core.flag & BAM_FREVERSE) != 0 ); }
    @property bool isReversed() { return bam_is_rev(this.b); }
    /// is read mapped?
    @property bool isMapped() { return (b.core.flag & BAM_FUNMAP) ==0; }
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
    /// Create cigar from bam read data
    @property Cigar cigar()
    {
        return Cigar(bam_get_cigar(b),(*b).core.n_cigar);
    }
    /// Get aux tag from bam record and return TagValue
    TagValue opIndex(string val)
    {
        char[2] arr;
        TagValue t;
        if(val.length==2){
            arr=val.dup[0..2];
            t=TagValue(b,arr);
        }
        return t;
    }
}

/**
Encapsulates a SAM/BAM file.
Implements InputRange interface using htslib calls.
*/
struct SAMFile {

    /// filename; as usable from D
    string filename;

    /// filename \0-terminated C string; reference needed to avoid GC reaping result of toStringz when ctor goes out of scope
    private immutable(char)* fn;

    /// htsFile
    private htsFile *fp;

    /// hFILE if required
    private hFILE *f;

    /// header struct
    bam_hdr_t *header = null;

    /// SAM/BAM/CRAM index 
    private hts_idx_t* idx;

    private kstring_t line;

    /// disallow copying
    @disable this(this);

    /** Create a representation of SAM/BAM/CRAM file from given filename

        fn:             string filename (complete path passed to htslib; may support S3:// https:// etc.)
        extra_threads:  extra threads for compression/decompression
            -1 use all available cores
            0  use no extra threads
            >1 add indicated number of threads (to a default of 1)
    */
    this(string fn, int extra_threads = -1)
    {
        import std.parallelism: totalCPUs;

        debug(dhtslib_debug) { writeln("SAMFile ctor"); }

        // open file
        this.filename = fn;
        this.fn = toStringz(fn);
        this.fp = hts_open(this.fn, cast(immutable(char)*)"r");

        if (extra_threads == -1 && totalCPUs > 1) {
            hts_log_info(__FUNCTION__, format("%d CPU cores detected; enabling multithreading", totalCPUs));
            // hts_set_threads adds N _EXTRA_ threads, so totalCPUs - 1 seemed reasonable,
            // but overcomitting by 1 thread (i.e., passing totalCPUs) buys an extra 3% on my 2-core 2013 Mac
            hts_set_threads(this.fp, totalCPUs );
        } else if (extra_threads > 0 ) {
            if ((extra_threads+1) > totalCPUs)
                hts_log_warning(__FUNCTION__, "More threads requested than CPU cores detected");
            hts_set_threads(this.fp, extra_threads);
        } else if (extra_threads == 0) {
            hts_log_debug(__FUNCTION__, "Zero extra threads requested");
        } else {
            hts_log_warning(__FUNCTION__, "Invalid negative number of extra threads requested");
        }

        // read header
        this.header = sam_hdr_read(this.fp);
        this.idx = sam_index_load(this.fp, this.fn);
        if (this.idx == null) {
            hts_log_info(__FUNCTION__, "SAM index not found");
            // TODO: attempt to build
            // TODO: edit range to return empty immediately if no idx
        }
        hts_log_debug( __FUNCTION__, format("SAM index: %s", this.idx));
    }
    this(File f, int extra_threads = -1)
    {
        import std.parallelism: totalCPUs;

        debug(dhtslib_debug) { writeln("SAMFile ctor"); }

        // open file
        this.filename = f.name();
        this.fn = toStringz(f.name);
        this.f = hdopen(f.fileno,cast(immutable(char)*)"r");
        this.fp = hts_hopen(this.f,this.fn,cast(immutable(char)*)"r");

        if (extra_threads == -1 && totalCPUs > 1) {
            hts_log_info(__FUNCTION__, format("%d CPU cores detected; enabling multithreading", totalCPUs));
            // hts_set_threads adds N _EXTRA_ threads, so totalCPUs - 1 seemed reasonable,
            // but overcomitting by 1 thread (i.e., passing totalCPUs) buys an extra 3% on my 2-core 2013 Mac
            hts_set_threads(this.fp, totalCPUs );
        } else if (extra_threads > 0 ) {
            if ((extra_threads+1) > totalCPUs)
                hts_log_warning(__FUNCTION__, "More threads requested than CPU cores detected");
            hts_set_threads(this.fp, extra_threads);
        } else if (extra_threads == 0) {
            hts_log_debug(__FUNCTION__, "Zero extra threads requested");
        } else {
            hts_log_warning(__FUNCTION__, "Invalid negative number of extra threads requested");
        }

        // read header
        this.header = sam_hdr_read(this.fp);
        this.idx = sam_index_load(this.fp, this.fn);
        if (this.idx == null) {
            hts_log_info(__FUNCTION__, "SAM index not found");
            // TODO: attempt to build
            // TODO: edit range to return empty immediately if no idx
        }
        hts_log_debug( __FUNCTION__, format("SAM index: %s", this.idx));
    }
    ~this()
    {
        debug(dhtslib_debug) { writeln("SAMFile dtor" ); }

        bam_hdr_destroy(this.header);

        //TODO:hts_close segfaults
        //We cant close the file pointer if another program has given us the pointer
        if(this.f !is null){
            //Causes double free
            //auto close=hclose(this.f);
            //if (ret != 0) writeln("There was an error hfile");
        }else{
            const auto ret = hts_close(fp);
            if (ret < 0) writefln("There was an error closing %s", fromStringz(this.fn));
        }
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
        auto itr = sam_itr_querys(this.idx, this.header, toStringz(q));
        return RecordRange(this.fp, itr);
    }

    /// Query by contig id, start, end
    auto query(int tid, int start, int end)
    {
        auto itr = sam_itr_queryi(this.idx, tid, start, end);
        return RecordRange(this.fp, itr);
    }

    /// Query by ["chr1:1-2","chr1:1000-1001"]
    auto query(string[] regions)
    {
        return RecordRangeMulti(this.fp,this.idx,this.header,&(this),regions);
    }

    /// bam["chr1:1-2"]
    auto opIndex(string q)
    {
        return query(q);
    }

    /// bam["chr1",1..2]
    auto opIndex(string tid, int[2] pos)
    {
        return query(tid,pos[0],pos[1]);
    }

    /// bam["chr1",1]
    auto opIndex(string tid, int pos)
    {
        return query(tid,pos,pos+1);
    }

    /// bam["chr1",1,2]
    deprecated
    auto opIndex(string tid, int pos1, int pos2)
    {
        return query(tid,pos1,pos2);
    }

    /// Integer-based chr below
    /// bam[0,1..2]
    auto opIndex(int tid, int[2] pos)
    {
        return query(tid,pos[0],pos[1]);
    }

    /// bam[0,1]
    auto opIndex(int tid, int pos)
    {
        return query(tid,pos,pos+1);
    }

    /// bam[0,1,2]
    deprecated
    auto opIndex(int tid, int pos1, int pos2)
    {
        return query(tid,pos1,pos2);
    }

    int[2] opSlice(size_t dim)(int start,int end)
        if(dim>=0 && dim<2)
    {
        return [start,end];
    }


    /// Return an InputRange representing all recods in the SAM/BAM/CRAM
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
        int success;
        ~this()
        {
            //debug(dhtslib_debug) hts_log_debug(__FUNCTION__, "dtor");
        }

        /// InputRange interface
        @property bool empty()
        {
            //    int sam_read1(samFile *fp, bam_hdr_t *h, bam1_t *b);
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
            success = sam_read1(this.fp, this.header, this.b);
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

    /// Iterate over records falling within a queried region (TODO: itr_multi_query)
    struct RecordRange
    {
        private htsFile     *fp;
        private hts_itr_t   *itr;
        private bam1_t      *b;

        private int r;

        ///
        this(htsFile * fp, hts_itr_t *itr)
        {
            this.itr = itr;
            this.fp = fp;
            b = bam_init1();
            //debug(dhtslib_debug) { writeln("sam_itr null? ",(cast(int)itr)==0); }
            hts_log_debug(__FUNCTION__, format("SAM itr null?: %s", cast(int)itr == 0));
            popFront();
        }

        /// InputRange interface
        @property bool empty()
        {
            return (r <= 0 && itr.finished) ? true : false;
        }
        /// ditto
        void popFront()
        {
            r = sam_itr_next(this.fp, this.itr, this.b);
        }
        /// ditto
        SAMRecord front()
        {
            return new SAMRecord(bam_dup1(b));
        }
    }

    /// Iterate over records falling within queried regions using a RegionList
    struct RecordRangeMulti
    {
        private htsFile             *fp;
        private hts_itr_multi_t     *itr;
        private bam1_t              *b;
        private hts_reglist_t[] rlist;
        private int r;

        ///
        this(htsFile * fp,hts_idx_t *idx, bam_hdr_t *header,SAMFile * sam,string[] regions)
        {
            rlist=RegionList(sam,regions).getRegList();
            this.fp = fp;
            b = bam_init1();
            itr =sam_itr_regions(idx,header,rlist.ptr,cast(uint) rlist.length);
            //debug(dhtslib_debug) { writeln("sam_itr null? ",(cast(int)itr)==0); }
            hts_log_debug(__FUNCTION__, format("SAM itr null?: %s", cast(int)itr == 0));
            popFront();
        }

        /// InputRange interface
        @property bool empty()
        {
            return (r <= 0 && itr.finished) ? true : false;
        }
        /// ditto
        void popFront()
        {
            r = sam_itr_multi_next(this.fp, this.itr, this.b);
        }
        /// ditto
        SAMRecord front()
        {
            return new SAMRecord(bam_dup1(b));
        }
    }

    /// List of regions based on sam/bam
    struct RegionList
    {
        import std.algorithm.iteration:splitter;
        import std.algorithm.sorting:sort;
        import std.range:drop,array;
        import std.conv:to;
        hts_reglist_t[string] rlist;
        SAMFile * sam;

        ///
        this(SAMFile * sam,string[] queries)
        {
            this.sam=sam;
            foreach(q;queries){
                addRegion(q);
            }
        }
        void addRegion(string reg)
        {
            //chr:1-2
            //split into chr and 1-2
            auto split=reg.splitter(":");
            auto chr=split.front;
            //split 1-2 into 1 and 2
            split=split.drop(1).front.splitter("-");
            if(sam.target_id(chr)<0){
                hts_log_error(__FUNCTION__,"tid not present in sam/bam");
            }
            addRegion(sam.target_id(chr),split.front.to!int,split.drop(1).front.to!int);
        }
        void addRegion(int tid, int beg, int end)
        {
            if(tid >this.sam.n_targets || tid<0)
                hts_log_error(__FUNCTION__,"tid not present in sam/bam");

            auto val =(this.sam.target_names[tid] in this.rlist);
            hts_pair32_t p;

            if(beg<0)
                hts_log_error(__FUNCTION__,"first coordinate < 0");

            if(beg>=this.sam.target_len(tid))
                hts_log_error(__FUNCTION__,"first coordinate larger than tid length");

            if(end<0)
                hts_log_error(__FUNCTION__,"second coordinate < 0");

            if(end>=this.sam.target_len(tid))
                hts_log_error(__FUNCTION__,"second coordinate larger than tid length");

            p.beg=beg;
            p.end=end;
            hts_pair32_t[] plist;
            if(val is null){
                hts_reglist_t r;

                //set tid
                r.tid=tid;

                //create intervals
                plist=plist~p;
                r.intervals=plist.ptr;
                r.count=cast(uint)plist.length;
                r.min_beg=p.beg;
                r.max_end=p.end;
                this.rlist[this.sam.target_names[tid]]=r;
            }else{
                plist=(val.intervals[0..val.count]~p).sort!(cmpInterval).array;
                val.intervals=plist.ptr;
                val.count=cast(uint)plist.length;
                val.min_beg=plist[0].beg;
                val.max_end=plist[$-1].end;
            }
        }
        hts_reglist_t[] getRegList()
        {
            return rlist.byValue.array.sort!(cmpRegList).array;
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

bool cmpInterval(hts_pair32_t a,hts_pair32_t b){
    if(a.beg < b.beg){
        return true;
    }
    if(a.end < b.end){
        return true;
    }
    return false;
}
bool cmpRegList(hts_reglist_t a,hts_reglist_t b){
    if(a.tid < b.tid){
        return true;
    }
    return false;
}

int parseSam(string line,bam_hdr_t * header,bam1_t * b){
    import dhtslib.htslib.kstring:kstring_t;
    import std.utf:toUTFz;
    kstring_t k;
    k.s=toUTFz!(char *)(line.dup);
    k.m=line.length+1;
    k.l=line.length+1;
    return sam_parse1(&k,header,b);
}
unittest{
    import std.stdio;
    import std.range:drop;
    import std.utf:toUTFz;
    import dhtslib.htslib.hts_log;
    hts_set_log_level(htsLogLevel.HTS_LOG_TRACE);
    hts_log_info(__FUNCTION__,"Loading sam file");
    auto range=File("htslib/test/realn01_exp-a.sam").byLineCopy();
    auto b=bam_init1();
    auto hdr=bam_hdr_init();
    string hdr_str;
    for(auto i=0;i<4;i++){
        hdr_str~=range.front~"\n";
        range.popFront;
    }
    hts_log_info(__FUNCTION__,"Header");
    writeln(hdr_str);
    hdr=sam_hdr_parse(cast(int)hdr_str.length,toUTFz!(char *)(hdr_str));
    hts_log_info(__FUNCTION__,"Read status:"~parseSam(range.front,hdr,b).to!string);
    auto r=new SAMRecord(b);
    hts_log_info(__FUNCTION__,"Cigar"~r.cigar.toString);
    assert(r.cigar.toString=="6M1D117M5D28M");
}