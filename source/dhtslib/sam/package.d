/**

SAMRecord and SAMFile are wrappers for htslib functions relating to SAM/BAM/CRAM* files

SAMRecord is a structured representation of a SAM/BAM/CRAM* record,
backed internally by htslib's bam1_t, but with convenient getters and setters
for record attributes, including functions not included in vanilla htslib 
like returning the sequence or the qscore string (NB: actually char*)

SAMFile is a structured representation of SAM/BAM/CRAM* file,
backed internally by htslib's htsFile and bam_hdr_t,
but with convenient getters and setters for file and header attributes,
as well as query functions accessible explicitly (`query("chr1:999-9999"`)
and by indexing (`samfile["chr1", 999 .. 9999]`).
The file object can be iterated as an InputRange to obtain every record in the file.

Authors: James S Blachly, MD <james.blachly@gmail.com> ; Thomas Gregory <charles.gregory@osumc.edu>

Bugs: SAMRecord and SAMFile function only as readers, rather than writers (i.e, cannot build SAMFile)
Bugs: (*CRAM functionality is limited and untested)

Date: 2020-09-12

License: Apache 2.0

Standards: Sequence Alignment/Map Format Specification v1 14 Dec 2018 http://samtools.github.io/hts-specs/

*/
module dhtslib.sam;

public import dhtslib.sam.header;

import core.stdc.stdlib : calloc, free,realloc;
import core.stdc.string:memset;
import std.format;
import std.parallelism : totalCPUs;
import std.stdio : writeln, writefln, stderr, File;
import std.string : fromStringz, toStringz;
import std.traits : ReturnType;
import std.algorithm : filter;

import htslib.hts : htsFile, hts_open, hts_close, hts_hopen;
import htslib.hts : hts_idx_t, hts_itr_t, hts_itr_multi_t, hts_reglist_t, hts_pair32_t;
import htslib.hts : seq_nt16_str,seq_nt16_table;
import htslib.hts : hts_set_threads;
import htslib.hfile : hdopen, hclose, hFILE;

import htslib.hts_log;
import htslib.kstring;
import htslib.sam;
import dhtslib.cigar;
import dhtslib.tagvalue;

/**
Encapsulates a SAM/BAM/CRAM record,
using the bam1_t type for memory efficiency,
and the htslib helper functions for speed.
**/
class SAMRecord
{
    ///
    bam1_t* b;

    /// Construct blank SAMRecord with empty backing bam1_t
    this()
    {
        //debug(dhtslib_debug) hts_log_debug(__FUNCTION__, "ctor()"); /// This line triggers memory error when __FUNCTION__, but not when "Other string"
        //test_log(__FUNCTION__, "ctor()");   /// This line will also trigger the memory error when __FUNCTION__, but not other strings
        //writeln(__FUNCTION__);    // this will not trigger the memory error
        this.b = bam_init1();
        //assert(0);                // This will elide(?) the memory error
        //assert(1 == 2);           // This will elide(?) the memory error
    }
    
    /// Construct SAMRecord from supplied bam1_t
    this(bam1_t* b)
    {
        //debug(dhtslib_debug) hts_log_debug(__FUNCTION__, "ctor(bam1_t *b)");
        this.b = b;
    }

    ~this()
    {
        //debug(dhtslib_debug) hts_log_debug(__FUNCTION__, "dtor");
        bam_destroy1(this.b); // we created our own in default ctor, or received copy via bam_dup1
    }


    /* bam1_core_t fields */

    /// chromosome ID, defined by bam_hdr_t
    pragma(inline, true)
    @nogc @safe nothrow
    @property int tid() { return this.b.core.tid; }
    /// ditto
    pragma(inline, true)
    @nogc @safe nothrow
    @property void tid(int tid) { this.b.core.tid = tid; }

    /// 0-based leftmost coordinate
    pragma(inline, true)
    @nogc @safe nothrow
    @property long pos() { return this.b.core.pos; }
    /// ditto
    pragma(inline, true)
    @nogc @safe nothrow
    @property void pos(long pos) { this.b.core.pos = pos; }

    // TODO: @field  bin     bin calculated by bam_reg2bin()

    /// mapping quality
    pragma(inline, true)
    @nogc @safe nothrow
    @property ubyte qual() { return this.b.core.qual; }
    /// ditto
    pragma(inline, true)
    @nogc @safe nothrow
    @property void qual(ubyte q) { this.b.core.qual = q; }

    /// length of query name. Terminated by 1-4 \0, see l_extranul. Never set directly.
    pragma(inline, true)
    @nogc @safe nothrow
    @property int l_qname() { return this.b.core.l_qname; }

    /// number of EXTRA trailing \0 in qname; 0 <= l_extranul <= 3, see l_qname. Never set directly.
    pragma(inline, true)
    @nogc @safe nothrow
    @property int l_extranul() { return this.b.core.l_extranul; }

    /// bitwise flag
    pragma(inline, true)
    @nogc @safe nothrow
    @property ushort flag() { return this.b.core.flag; }
    /// ditto
    pragma(inline, true)
    @nogc @safe nothrow
    @property void flag(ushort fl) { this.b.core.flag = fl; }

    /// is read paired?
    pragma(inline, true)
    @property bool isPaired()
    {
        return (b.core.flag & BAM_FPAIRED);
    }
    
    /// is read reversed?
    /// bool bam_is_rev(bam1_t *b) { return ( ((*b).core.flag & BAM_FREVERSE) != 0 ); }
    @property bool isReversed()
    {
        version(LDC) pragma(inline, true);
        return bam_is_rev(this.b);
    }

    /// is read mapped?
    @property bool isMapped()
    {
        version(LDC){
            pragma(inline, true);
        }
        return (b.core.flag & BAM_FUNMAP) == 0;
    }

    /// is read mate mapped?
    @property bool isMateMapped()
    {
        version(LDC){
            pragma(inline, true);
        }
        return (b.core.flag & BAM_FMUNMAP) == 0;
    }

    /// is mate reversed?
    /// bool bam_is_mrev(bam1_t *b) { return( ((*b).core.flag & BAM_FMREVERSE) != 0); }
    pragma(inline, true)
    @property bool mateReversed()
    {
        return bam_is_mrev(this.b);
    }
    
    /// is read secondary?
    @property bool isSecondary()
    {
        version(LDC){
            pragma(inline, true);
        }
        return cast(bool)(b.core.flag & BAM_FSECONDARY);
    }

    /// is read supplementary?
    @property bool isSupplementary()
    {
        version(LDC){
            pragma(inline, true);
        }
        return cast(bool)(b.core.flag & BAM_FSUPPLEMENTARY);
    }

    /// Return read strand + or - (as char)
    @property char strand(){
        return ['+','-'][isReversed()];
    }

    /// auto bam_get_qname(bam1_t *b) { return (cast(char*)(*b).data); }
    pragma(inline, true)
    @property char[] queryName()
    {
        return fromStringz(bam_get_qname(this.b));
    }
    
    /// Set query name
    pragma(inline, true)
    @property void queryName(string qname)
    {
        assert(qname.length<252);

        //make cstring
        auto qname_n=qname.dup~'\0';

        //append extra nulls to 32bit align cigar
        auto extranul=0;
        for (; qname_n.length % 4 != 0; extranul++) qname_n~='\0';
        assert(extranul<=3);

        //get length of rest of data
        auto l_rest=b.l_data-b.core.l_qname;

        //copy rest of data
        ubyte[] rest=b.data[b.core.l_qname..b.l_data].dup;
        
        //if not enough space
        if(qname_n.length+l_rest>b.m_data){
            auto ptr=cast(ubyte*)realloc(b.data,qname_n.length+l_rest);
            assert(ptr!=null);
            b.data=ptr;
            b.m_data=cast(uint)qname_n.length+l_rest;
        }

        //reassign q_name and rest of data
        b.data[0..qname_n.length]=(cast(ubyte[])qname_n);
        b.data[qname_n.length..qname_n.length+l_rest]=rest;

        //reset data length, qname length, extra nul length
        b.l_data=cast(uint)qname_n.length+l_rest;
        b.core.l_qname=cast(ubyte)(qname_n.length);
        b.core.l_extranul=cast(ubyte)extranul;
    }

    /// query (and quality string) length
    pragma(inline, true)
    @property int length()
    {
        return this.b.core.l_qseq;
    }

    /// Return char array of the sequence
    /// see samtools/sam_view.c: get_read
    @property const(char)[] sequence()
    {
        // calloc fills with \0; +1 len for Cstring
        // char* s = cast(char*) calloc(1, this.b.core.l_qseq + 1);
        char[] s;
        s.length = this.b.core.l_qseq;

        // auto bam_get_seq(bam1_t *b) { return ((*b).data + ((*b).core.n_cigar<<2) + (*b).core.l_qname); }
        auto seqdata = bam_get_seq(this.b);

        for (int i; i < this.b.core.l_qseq; i++)
        {
            s[i] = seq_nt16_str[bam_seqi(seqdata, i)];
        }
        return s;
    }

    /// Assigns sequence and resets quality score
    pragma(inline, true)
    @property void sequence(const(char)[] seq)
    {
        //nibble encode sequence
        ubyte[] en_seq;
        en_seq.length=(seq.length+1)>>1;
        for(auto i=0;i<seq.length;i++){
            en_seq[i>>1]|=seq_nt16_table[seq[i]]<< ((~i&1)<<2);
        }

        //get length of data before seq
        uint l_prev=b.core.l_qname + cast(uint)(b.core.n_cigar*uint.sizeof);

        //get starting point after seq
        uint start_rest=l_prev+((b.core.l_qseq+1)>>1)+b.core.l_qseq;

        //copy previous data
        ubyte[] prev=b.data[0..l_prev].dup;

        //copy rest of data
        ubyte[] rest=b.data[start_rest..b.l_data].dup;
        
        //if not enough space
        if(en_seq.length+seq.length+l_prev+(b.l_data-start_rest)>b.m_data){
            auto ptr=cast(ubyte*)realloc(b.data,en_seq.length+seq.length+l_prev+(b.l_data-start_rest));
            assert(ptr!=null);
            b.data=ptr;
            b.m_data=cast(uint)en_seq.length+cast(uint)seq.length+l_prev+(b.l_data-start_rest);
        }

        //reassign q_name and rest of data
        b.data[0..l_prev]=prev;
        b.data[l_prev..en_seq.length+l_prev]=en_seq;
        //set qscores to 255 
        memset(b.data+en_seq.length+l_prev,255,seq.length);
        b.data[l_prev+en_seq.length+seq.length..l_prev+en_seq.length+seq.length+(b.l_data-start_rest)]=rest;

        //reset data length, seq length
        b.l_data=cast(uint)en_seq.length+cast(uint)seq.length+l_prev+(b.l_data-start_rest);
        b.core.l_qseq=cast(int)(seq.length);
    }

    /// Return char array of the qscores
    /// see samtools/sam_view.c: get_quality
    const(char)[] qscores(bool phredScale=true)()
    {
        // calloc fills with \0; +1 len for Cstring
        // char* q = cast(char*) calloc(1, this.b.core.l_qseq + 1);
        char[] q;
        q.length = this.b.core.l_qseq;

        // auto bam_get_qual(bam1_t *b) { return (*b).data + ((*b).core.n_cigar<<2) + (*b).core.l_qname + (((*b).core.l_qseq + 1)>>1); }
        char* qualdata = cast(char*) bam_get_qual(this.b);

        for (int i; i < this.b.core.l_qseq; i++)
            static if(phredScale) q[i] = cast(char)(qualdata[i] + 33);
            else q[i] = cast(char)(qualdata[i]);
        return q;
    }

    /// Add qscore sequence given that it is the same length as the bam sequence
    pragma(inline, true)
    void q_scores(bool phredScale=true)(const(char)[] seq)
    {
        if(seq.length!=b.core.l_qseq){
            hts_log_error(__FUNCTION__,"qscore length does not match sequence length");
            return;
        }
        auto s= seq.dup;
        for(auto i=0;i<seq.length;i++){
            static if(phredScale) s[i]=seq[i]-33;
        }
        
        //get length of data before seq
        uint l_prev=b.core.l_qname + cast(uint)(b.core.n_cigar*uint.sizeof)+((b.core.l_qseq+1)>>1);
        b.data[l_prev..s.length+l_prev]=cast(ubyte[])s;
    }

    /// Add qscore sequence given that it is the same length as the bam sequence
    pragma(inline, true)
    void qscores(ubyte[] seq)
    {
        if(seq.length!=b.core.l_qseq){
            hts_log_error(__FUNCTION__,"qscore length does not match sequence length");
            return;
        }
        
        //get length of data before seq
        uint l_prev=b.core.l_qname + cast(uint)(b.core.n_cigar*uint.sizeof)+((b.core.l_qseq+1)>>1);
        b.data[l_prev..seq.length+l_prev]=seq;
    }

    /// Create cigar from bam1_t record
    @property Cigar cigar()
    {
        return Cigar(bam_get_cigar(b), (*b).core.n_cigar);
    }

    /// Assign a cigar string
    pragma(inline, true)
    @property void cigar(Cigar cigar)
    {
        
        //get length of data before seq
        uint l_prev=b.core.l_qname;

        //get starting point after seq
        uint start_rest=l_prev+cast(uint)(b.core.n_cigar*uint.sizeof);

        //copy previous data
        ubyte[] prev=b.data[0..l_prev].dup;

        //copy rest of data
        ubyte[] rest=b.data[start_rest..b.l_data].dup;
        
        //if not enough space
        if(uint.sizeof*cigar.ops.length+l_prev+(b.l_data-start_rest)>b.m_data){
            auto ptr=cast(ubyte*)realloc(b.data,uint.sizeof*cigar.ops.length+l_prev+(b.l_data-start_rest));
            assert(ptr!=null);
            b.data=ptr;
            b.m_data=cast(uint)(uint.sizeof*cigar.ops.length)+l_prev+(b.l_data-start_rest);
        }

        //reassign q_name and rest of data
        b.data[0..l_prev]=prev;
        //without cigar.ops.dup we get the overlapping array copy error
        b.data[l_prev..uint.sizeof*cigar.ops.length+l_prev]=cast(ubyte[])(cast(uint[])cigar.ops.dup);
        b.data[l_prev+uint.sizeof*cigar.ops.length..l_prev+uint.sizeof*cigar.ops.length+(b.l_data-start_rest)]=rest;

        //reset data length, seq length
        b.l_data=cast(uint)(uint.sizeof*cigar.ops.length)+l_prev+(b.l_data-start_rest);
        b.core.n_cigar=cast(int)(cigar.ops.length);
    }

    /// Get aux tag from bam1_t record and return a TagValue
    TagValue opIndex(string val)
    {
        return TagValue(b, val[0..2]);
    }

    /// Assign a tag key value pair
    void opIndexAssign(T)(T value,string index)
    {
        import std.traits:isIntegral,isSomeString;
        static if(isIntegral!T){
            bam_aux_update_int(b,index[0..2],value);
        }else static if(is(T==float)){
            bam_aux_update_float(b,index[0..2],value);
        }else static if(isSomeString!T){
            bam_aux_update_str(b,index[0..2],cast(int)value.length+1,toStringz(value));
        }
    }

    /// Assign a tag key array pair
    void opIndexAssign(T:T[])(T[] value,string index)
    {
        bam_aux_update_array(b,index[0..2],TypeChars[TypeIndex!T],value.length,value.ptr);
    }

    /// chromosome ID of next read in template, defined by bam_hdr_t
    pragma(inline, true)
    @nogc @safe nothrow
    @property int mateTID() { return this.b.core.mtid; }
    /// ditto
    pragma(inline, true)
    @nogc @safe nothrow
    @property void mateTID(int mtid) { this.b.core.mtid = mtid; }

    /// 0-based leftmost coordinate of next read in template
    pragma(inline, true)
    @nogc @safe nothrow
    @property long matePos() { return this.b.core.mpos; }
    /// ditto
    pragma(inline, true)
    @nogc @safe nothrow
    @property void matePos(long mpos) { this.b.core.mpos = mpos; }

    /// Presumably Insert size, but is undocumented.
    /// Per samtools source, is measured 5' to 5'
    /// https://github.com/samtools/samtools/blob/bd1a409aa750d25d70a093405d174db8709a3f2c/bam_mate.c#L320
    pragma(inline, true)
    @nogc @safe nothrow
    @property long insertSize() { return this.b.core.isize; }
    /// ditto
    pragma(inline, true)
    @nogc @safe nothrow
    @property void insertSize(long isize) { this.b.core.isize = isize; }
/+ BEGIN CHERRY-PICK: TODO: confirm below chunk with htslib-1.10 (mainly long coordinates +/
    /// get aligned coordinates per each cigar op
    auto getAlignedCoordinates(){
        return AlignedCoordinatesItr(this.cigar);
    }

    /// get aligned coordinates per each cigar op within range
    /// range is 0-based half open using chromosomal coordinates
    auto getAlignedCoordinates(long start, long end){
        assert(start >= this.pos);
        assert(end <= this.pos + this.cigar.ref_bases_covered);
        return this.getAlignedCoordinates.filter!(x=> this.pos + x.rpos >= start).filter!(x=> this.pos + x.rpos < end);
    }

    struct AlignedPair(bool refSeq)
    {
        long qpos, rpos;
        Ops cigar_op;
        char queryBase;
        static if(refSeq) char refBase;

        string toString(){
            import std.conv : to;
            static if(refSeq) return rpos.to!string ~ ":" ~ refBase ~ " > " ~ qpos.to!string ~ ":" ~ queryBase;
            else return rpos.to!string ~ ":" ~ qpos.to!string ~ ":" ~ queryBase;
        }
    }

    /// get a range of aligned read and reference positions
    /// this is meant to recreate functionality from pysam:
    /// https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.get_aligned_pairs
    /// range is 0-based half open using chromosomal coordinates
    auto getAlignedPairs(bool withRefSeq)(long start, long end)
    {
        import dhtslib.md : MDItr;
        import std.range : dropExactly;
        struct AlignedPairs(bool refSeq)
        {
            ReturnType!getAlignedCoordinates coords;
            static if(refSeq) MDItr mdItr;
            AlignedPair!refSeq current;
            const(char)[] qseq;

            auto getAlignedCoordinates(SAMRecord rec, long start, long end){
                return rec.getAlignedCoordinates(start,end);
            }
            
            this(SAMRecord rec, long start, long end){
                assert(start < end);
                assert(start >= rec.pos);
                assert(end <= rec.pos + rec.cigar.ref_bases_covered);
                
                coords = getAlignedCoordinates(rec, start, end);
                qseq = rec.sequence;
                static if(refSeq){
                    assert(rec["MD"].exists,"MD tag must exist for record");
                    mdItr = MDItr(rec);
                    mdItr = mdItr.dropExactly(start - rec.pos);
                } 
                current.qpos = coords.front.qpos;
                current.rpos = coords.front.rpos;
                current.cigar_op = coords.front.cigar_op;
                if(!CigarOp(0, current.cigar_op).is_query_consuming) current.queryBase = '-';
                else current.queryBase = qseq[current.qpos];
                
                static if(refSeq){
                    if(!CigarOp(0, current.cigar_op).is_reference_consuming) current.refBase = '-';
                    else if(mdItr.front == '=') current.refBase = current.queryBase;
                    else current.refBase = mdItr.front;
                }
            }

            AlignedPair!refSeq front(){
                return current;
            }

            void popFront(){
                coords.popFront;
                if(coords.empty) return;
                current.qpos = coords.front.qpos;
                current.rpos = coords.front.rpos;
                current.cigar_op = coords.front.cigar_op;
                if(!CigarOp(0, current.cigar_op).is_query_consuming) current.queryBase = '-';
                else current.queryBase = qseq[current.qpos];
                
                static if(refSeq){
                    if(CigarOp(0, current.cigar_op).is_reference_consuming) mdItr.popFront;
                    if(!CigarOp(0, current.cigar_op).is_reference_consuming) current.refBase = '-';
                    else if(mdItr.front == '=') current.refBase = current.queryBase;
                    else current.refBase = mdItr.front;
                }

            }

            bool empty(){
                return coords.empty;
            }

        }
        return AlignedPairs!withRefSeq(this,start,end);
    }
    /// get a range of aligned read and reference positions
    /// this is meant to recreate functionality from pysam:
    /// https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.get_aligned_pairs
    auto getAlignedPairs(bool withRefSeq)(){
        return getAlignedPairs!withRefSeq(this.pos, this.pos + this.cigar.ref_bases_covered);
    }
/+ END CHERRY-PICK +/
}

debug(dhtslib_unittest)
unittest{
    writeln();
    import dhtslib.sam;
    import htslib.hts_log : hts_log_info;
    import std.path:buildPath,dirName;
    hts_set_log_level(htsLogLevel.HTS_LOG_TRACE);
    hts_log_info(__FUNCTION__, "Testing SAMRecord mutation");
    hts_log_info(__FUNCTION__, "Loading test file");
    auto bam = SAMFile(buildPath(dirName(dirName(dirName(dirName(__FILE__)))),"htslib","test","range.bam"), 0);
    auto read=bam.all_records.front;

    //change queryname
    hts_log_info(__FUNCTION__, "Testing queryname");
    assert(read.queryName=="HS18_09653:4:1315:19857:61712");
    read.queryName="HS18_09653:4:1315:19857:6171";
    assert(read.queryName=="HS18_09653:4:1315:19857:6171");
    assert(read.cigar.toString() == "78M1D22M");
    read.queryName="HS18_09653:4:1315:19857:617122";
    assert(read.queryName=="HS18_09653:4:1315:19857:617122");
    assert(read.cigar.toString() == "78M1D22M");
    assert(read["RG"].check!string);
    assert(read["RG"].to!string=="1");

    //change sequence
    hts_log_info(__FUNCTION__, "Testing sequence");
    assert(read.sequence=="AGCTAGGGCACTTTTTGTCTGCCCAAATATAGGCAACCAAAAATAATTTCCAAGTTTTTAATGATTTGTTGCATATTGAAAAAACATTTTTTGGGTTTTT");
    read.sequence=read.sequence[0..10];
    assert(read.sequence=="AGCTAGGGCA");
    ubyte[] q=[255,255,255,255,255,255,255,255,255,255];
    assert(read.qscores!false==cast(char[])(q));
    q=[1,14,20,40,30,10,14,20,40,30];
    read.qscores(q);
    assert(read.qscores!false==cast(char[])(q));
    q[]+=33;
    assert(read.qscores==cast(char[])(q));
    assert(read.cigar.toString() == "78M1D22M");
    assert(read["RG"].check!string);
    assert(read["RG"].to!string=="1");

    //change cigar
    hts_log_info(__FUNCTION__, "Testing cigar");
    auto c=read.cigar;
    c.ops[$-1].length=21;
    read.cigar=c;
    assert(read.cigar.toString() == "78M1D21M");
    assert(read["RG"].check!string);
    assert(read["RG"].to!string=="1");
    
    //change tag
    hts_log_info(__FUNCTION__, "Testing tags");
    read["X0"]=2;
    assert(read["X0"].to!ubyte==2);
    assert(read["RG"].check!string);
    
    read["RG"]="test";
    assert(read["RG"].to!string=="test");
    hts_log_info(__FUNCTION__, "Cigar:" ~ read.cigar.toString());
}

debug(dhtslib_unittest)
unittest{
    writeln();
    import dhtslib.sam;
    import htslib.hts_log : hts_log_info;
    import std.path:buildPath,dirName;
    hts_set_log_level(htsLogLevel.HTS_LOG_TRACE);
    hts_log_info(__FUNCTION__, "Testing SAMRecord mutation w/ realloc");
    hts_log_info(__FUNCTION__, "Loading test file");
    auto bam = SAMFile(buildPath(dirName(dirName(dirName(dirName(__FILE__)))),"htslib","test","range.bam"), 0);
    auto read=bam.all_records.front;

    //change queryname
    hts_log_info(__FUNCTION__, "Testing queryname");
    assert(read.queryName=="HS18_09653:4:1315:19857:61712");

    //we do this only to force realloc
    read.b.m_data=read.b.l_data;


    auto prev_m=read.b.m_data;
    read.queryName="HS18_09653:4:1315:19857:61712HS18_09653:4:1315:19857:61712";
    assert(read.b.m_data >prev_m);
    assert(read.queryName=="HS18_09653:4:1315:19857:61712HS18_09653:4:1315:19857:61712");
    assert(read.cigar.toString() == "78M1D22M");
    assert(read["RG"].check!string);
    assert(read["RG"].to!string=="1");

    //change sequence
    hts_log_info(__FUNCTION__, "Testing sequence");
    assert(read.sequence=="AGCTAGGGCACTTTTTGTCTGCCCAAATATAGGCAACCAAAAATAATTTCCAAGTTTTTAATGATTTGTTGCATATTGAAAAAACATTTTTTGGGTTTTT");
    prev_m=read.b.m_data;
    read.sequence="AGCTAGGGCACTTTTTGTCTGCCCAAATATAGGCAACCAAAAATAATTTCCAAGTTTTTAATGATTTGTTGCATATTGAAAAAACATTTTTTGGGTTTTTAGCTAGGGCACTTTTTGTCTGCCCAAATATAGGCAACCAAAAATAATTTCCAAGTTTTTAATGATTTGTTGCATATTGAAAAAACATTTTTTGGGTTTTT";
    assert(read.b.m_data >prev_m);
    assert(read.sequence=="AGCTAGGGCACTTTTTGTCTGCCCAAATATAGGCAACCAAAAATAATTTCCAAGTTTTTAATGATTTGTTGCATATTGAAAAAACATTTTTTGGGTTTTTAGCTAGGGCACTTTTTGTCTGCCCAAATATAGGCAACCAAAAATAATTTCCAAGTTTTTAATGATTTGTTGCATATTGAAAAAACATTTTTTGGGTTTTT");
    assert(read.cigar.toString() == "78M1D22M");
    assert(read["RG"].check!string);
    assert(read["RG"].to!string=="1");

    //change cigar
    hts_log_info(__FUNCTION__, "Testing cigar");
    auto c=read.cigar;
    c.ops=c.ops~c.ops;
    prev_m=read.b.m_data;
    read.cigar=c;
    assert(read.b.m_data >prev_m);
    assert(read.cigar.toString() == "78M1D22M78M1D22M");
    assert(read["RG"].check!string);
    assert(read["RG"].to!string=="1");
    
    //change tag
    hts_log_info(__FUNCTION__, "Testing tags");
    read["X0"]=2;
    assert(read["X0"].to!ubyte==2);
    assert(read["RG"].check!string);
    
    read["RG"]="test";
    assert(read["RG"].to!string=="test");
    hts_log_info(__FUNCTION__, "Cigar:" ~ read.cigar.toString());
}

unittest
{
    import std.stdio;
    import dhtslib.sam;
    import dhtslib.md : MDItr;
    import std.algorithm: map;
    import std.array: array;
    import std.path:buildPath,dirName;
    auto bam = SAMFile(buildPath(dirName(dirName(dirName(dirName(__FILE__)))),"htslib","test","range.bam"), 0);
    auto read=bam.all_records.front;
    auto pairs = read.getAlignedPairs!true(read.pos + 77, read.pos + 77 + 5);
    assert(pairs.map!(x => x.qpos).array == [77,77,78,79,80]);
    assert(pairs.map!(x => x.rpos).array == [77,78,79,80,81]);
    assert(pairs.map!(x => x.refBase).array == "GAAAA");
    assert(pairs.map!(x => x.queryBase).array == "G-AAA");
}

alias SAMFile = SAMReader;
/**
Encapsulates a SAM/BAM file.

Implements InputRange interface using htslib calls.
If indexed, Random-access query via multidimensional slicing.
*/
struct SAMReader
{
    /// filename; as usable from D
    string filename;

    /// filename \0-terminated C string; reference needed to avoid GC reaping result of toStringz when ctor goes out of scope
    private immutable(char)* fn;

    /// htsFile
    private htsFile* fp;

    /// hFILE if required
    private hFILE* f;

    /// header struct
    bam_hdr_t* header = null;

    /// SAM/BAM/CRAM index 
    private hts_idx_t* idx;

    private kstring_t line;

    /// disallow copying
    @disable this(this);

    /** Create a representation of SAM/BAM/CRAM file from given filename or File

    Params:
        fn =            string filename (complete path passed to htslib; may support S3:// https:// etc.)
        extra_threads = extra threads for compression/decompression
                        -1 use all available cores (default)
                        0  use no extra threads
                        >1 add indicated number of threads (to a default of 1)
    */
    this(T)(T f, int extra_threads = -1)
    if (is(T == string) || is(T == File))
    {
        import std.parallelism : totalCPUs;

        // open file
        static if (is(T == string))
        {
            this.filename = f;
            this.fn = toStringz(f);
            this.fp = hts_open(this.fn, cast(immutable(char)*) "r");
        }
        else static if (is(T == File))
        {
            this.filename = f.name();
            this.fn = toStringz(f.name);
            this.f = hdopen(f.fileno, cast(immutable(char)*) "r");
            this.fp = hts_hopen(this.f, this.fn, cast(immutable(char)*) "r");
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
        }
        else if (extra_threads > 0)
        {
            if ((extra_threads + 1) > totalCPUs)
                hts_log_warning(__FUNCTION__, "More threads requested than CPU cores detected");
            hts_set_threads(this.fp, extra_threads);
        }
        else if (extra_threads == 0)
        {
            hts_log_debug(__FUNCTION__, "Zero extra threads requested");
        }
        else
        {
            hts_log_warning(__FUNCTION__, "Invalid negative number of extra threads requested");
        }

        // read header
        this.header = sam_hdr_read(this.fp);
        this.idx = sam_index_load(this.fp, this.fn);
        if (this.idx == null)
        {
            hts_log_info(__FUNCTION__, "SAM index not found");
            // TODO: attempt to build
            // TODO: edit range to return empty immediately if no idx
        }
        hts_log_debug(__FUNCTION__, format("SAM index: %s", this.idx));
    }

    ~this()
    {
        debug (dhtslib_debug)
        {
            writeln("SAMFile dtor");
        }

        bam_hdr_destroy(this.header);

        if((this.fp !is null) && (this.f is null))
        {
            const auto ret = hts_close(fp);
            if (ret < 0)
                writefln("There was an error closing %s", fromStringz(this.fn));
        }
    }

    /// number of reference sequences; from bam_hdr_t
    @property int n_targets() const
    {
        return this.header.n_targets;
    }

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
        foreach (i; 0 .. this.n_targets)
        {
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
    auto query(string chrom, long start, long end)
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
    auto query(int tid, long start, long end)
    {
        auto itr = sam_itr_queryi(this.idx, tid, start, end);
        return RecordRange(this.fp, itr);
    }

    /// Query by ["chr1:1-2","chr1:1000-1001"]
    auto query(string[] regions)
    {
        return RecordRangeMulti(this.fp, this.idx, this.header, &(this), regions);
    }

    /// bam["chr1:1-2"]
    auto opIndex(string q)
    {
        return query(q);
    }

    /// bam["chr1", 1..2]
    auto opIndex(string tid, int[2] pos)
    {
        return query(tid, pos[0], pos[1]);
    }

    /// bam["chr1", 1]
    auto opIndex(string tid, int pos)
    {
        return query(tid, pos, pos + 1);
    }

    /** Deprecated: use multidimensional slicing with second parameter as range (["chr1", 1 .. 2]) */
    /// bam["chr1", 1, 2]
    deprecated auto opIndex(string tid, int pos1, int pos2)
    {
        return query(tid, pos1, pos2);
    }

    /// Integer-based chr below
    /// bam[0, 1..2]
    auto opIndex(int tid, int[2] pos)
    {
        return query(tid, pos[0], pos[1]);
    }

    /// bam[0, 1]
    auto opIndex(int tid, int pos)
    {
        return query(tid, pos, pos + 1);
    }

    /// bam[0, 1, 2]
    deprecated auto opIndex(int tid, int pos1, int pos2)
    {
        return query(tid, pos1, pos2);
    }

    /// support bam["chr1", 1..2 ]
    int[2] opSlice(size_t dim)(int start, int end) if (dim  == 1)
    {
        return [start, end];
    }

    /// Return an InputRange representing all recods in the SAM/BAM/CRAM
    AllRecordsRange all_records()
    {
        auto range = AllRecordsRange(this.fp, this.header, bam_init1());
        range.popFront();
        return range;
    }
    alias allRecords = all_records;

    /// Iterate through all records in the SAM/BAM/CRAM
    struct AllRecordsRange
    {
        private htsFile*    fp;     // belongs to parent; shared
        private bam_hdr_t*  header; // belongs to parent; shared
        private bam1_t*     b;
        private int success;
        ~this()
        {
            //debug(dhtslib_debug) hts_log_debug(__FUNCTION__, "dtor");
        }

        /// InputRange interface
        @property bool empty() // @suppress(dscanner.suspicious.incorrect_infinite_range)
        {
            //    int sam_read1(samFile *fp, bam_hdr_t *h, bam1_t *b);
            if (success >= 0)
                return false;
            else if (success == -1)
                return true;
            else
            {
                hts_log_error(__FUNCTION__, "*** ERROR sam_read1 < -1");
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
        private htsFile* fp;
        private hts_itr_t* itr;
        private bam1_t* b;

        private int r;

        ///
        this(htsFile* fp, hts_itr_t* itr)
        {
            this.itr = itr;
            this.fp = fp;
            b = bam_init1();
            //debug(dhtslib_debug) { writeln("sam_itr null? ",(cast(int)itr)==0); }
            hts_log_debug(__FUNCTION__, format("SAM itr null?: %s", cast(int) itr == 0));
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
        private htsFile* fp;
        private hts_itr_multi_t* itr;
        private bam1_t* b;
        private hts_reglist_t[] rlist;
        private int r;

        ///
        this(htsFile* fp, hts_idx_t* idx, bam_hdr_t* header, SAMFile* sam, string[] regions)
        {
            rlist = RegionList(sam, regions).getRegList();
            this.fp = fp;
            b = bam_init1();
            itr = sam_itr_regions(idx, header, rlist.ptr, cast(uint) rlist.length);
            //debug(dhtslib_debug) { writeln("sam_itr null? ",(cast(int)itr)==0); }
            hts_log_debug(__FUNCTION__, format("SAM itr null?: %s", cast(int) itr == 0));
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
        import std.algorithm.iteration : splitter;
        import std.algorithm.sorting : sort;
        import std.range : drop, array;
        import std.conv : to;

        private hts_reglist_t[string] rlist;
        private SAMFile* sam;

        ///
        this(SAMFile* sam, string[] queries)
        {
            this.sam = sam;
            foreach (q; queries)
            {
                addRegion(q);
            }
        }

        /// Add a region in standard format (chr1:2000-3000) to the RegionList
        void addRegion(string reg)
        {
            //chr:1-2
            //split into chr and 1-2
            auto split = reg.splitter(":");
            auto chr = split.front;
            //split 1-2 into 1 and 2
            split = split.drop(1).front.splitter("-");
            if (sam.target_id(chr) < 0)
            {
                hts_log_error(__FUNCTION__, "tid not present in sam/bam");
            }
            addRegion(sam.target_id(chr), split.front.to!int, split.drop(1).front.to!int);
        }

        /// Add a region by {target/contig/chr id, start coord, end coord} to the RegionList
        void addRegion(int tid, int beg, int end)
        {
            if (tid > this.sam.n_targets || tid < 0)
                hts_log_error(__FUNCTION__, "tid not present in sam/bam");

            auto val = (this.sam.target_names[tid] in this.rlist);
            hts_pair32_t p;

            if (beg < 0)
                hts_log_error(__FUNCTION__, "first coordinate < 0");

            if (beg >= this.sam.target_len(tid))
                hts_log_error(__FUNCTION__, "first coordinate larger than tid length");

            if (end < 0)
                hts_log_error(__FUNCTION__, "second coordinate < 0");

            if (end >= this.sam.target_len(tid))
                hts_log_error(__FUNCTION__, "second coordinate larger than tid length");

            p.beg = beg;
            p.end = end;
            hts_pair32_t[] plist;
            if (val is null)
            {
                hts_reglist_t r;

                //set tid
                r.tid = tid;

                //create intervals
                plist = plist ~ p;
                r.intervals = plist.ptr;
                r.count = cast(uint) plist.length;
                r.min_beg = p.beg;
                r.max_end = p.end;
                this.rlist[this.sam.target_names[tid]] = r;
            }
            else
            {
                plist = (val.intervals[0 .. val.count] ~ p).sort!(cmpInterval).array;
                val.intervals = plist.ptr;
                val.count = cast(uint) plist.length;
                val.min_beg = plist[0].beg;
                val.max_end = plist[$ - 1].end;
            }
        }

        hts_reglist_t[] getRegList()
        {
            return rlist.byValue.array.sort!(cmpRegList).array;
        }
    }
}

/// Nucleotide complement table; from samtools/sam_view.c
private const(char)[16] seq_comp_table = [0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15];

/// Reverse a string in place; from samtools/sam_view.c
//  TODO: Could be sped up, and potentially made safer? by passing strlen since already known
private char* reverse(char* str)
{
    import core.stdc.string : strlen;

    auto i = strlen(str) - 1, j = 0;
    char ch;
    while (i > j)
    {
        ch = str[i];
        str[i] = str[j];
        str[j] = ch;
        i--;
        j++;
    }
    return str;
}

/// SAM/BAM/CRAM disk format
enum SAMWriterTypes{
    BAM,
    UBAM,
    SAM,
    CRAM,
    DEDUCE
}

/// Encapsulates a SAM/BAM/CRAM, but as write-only
struct SAMWriter
{
    /// filename; as usable from D
    string filename;

    /// filename \0-terminated C string; reference needed to avoid GC reaping result of toStringz when ctor goes out of scope
    private immutable(char)* fn;

    /// htsFile
    private htsFile* fp;

    /// hFILE if required
    private hFILE* f;

    /// header struct
    bam_hdr_t* header = null;

    private kstring_t line;

    /// disallow copying
    @disable this(this);

    /** Create a representation of SAM/BAM/CRAM file from given filename or File

    Params:
        fn =            string filename (complete path passed to htslib; may support S3:// https:// etc.)
        extra_threads = extra threads for compression/decompression
                        -1 use all available cores (default)
                        0  use no extra threads
                        >1 add indicated number of threads (to a default of 1)
    */
    this(T)(T f,bam_hdr_t * header, SAMWriterTypes t=SAMWriterTypes.DEDUCE,int extra_threads = -1)
    if (is(T == string) || is(T == File))
    {
        import std.parallelism : totalCPUs;
        char[] mode;
        if(t == SAMWriterTypes.BAM) mode=['w','b','\0'];
        else if(t == SAMWriterTypes.UBAM) mode=['w','b','u','\0'];
        else if(t == SAMWriterTypes.SAM) mode=['w','\0'];
        else if(t == SAMWriterTypes.CRAM) mode=['w','c','\0'];
        // open file
        static if (is(T == string))
        {
            if(t == SAMWriterTypes.DEDUCE){
                import std.path:extension;
                auto ext=extension(f);
                if(ext==".bam") mode=['w','b','\0'];
                else if(ext==".sam") mode=['w','\0'];
                else if(ext==".cram") mode=['w','c','\0'];
                else {
                    hts_log_error(__FUNCTION__,"extension "~ext~" not valid");
                    throw new Exception("DEDUCE SAMWriterType used with non-valid extension");
                }
            }
            this.filename = f;
            this.fn = toStringz(f);
            this.fp = hts_open(this.fn, mode.ptr);
        }
        else static if (is(T == File))
        {
            assert(t!=SAMWriterTypes.DEDUCE);
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
        }
        else if (extra_threads > 0)
        {
            if ((extra_threads + 1) > totalCPUs)
                hts_log_warning(__FUNCTION__, "More threads requested than CPU cores detected");
            hts_set_threads(this.fp, extra_threads);
        }
        else if (extra_threads == 0)
        {
            hts_log_debug(__FUNCTION__, "Zero extra threads requested");
        }
        else
        {
            hts_log_warning(__FUNCTION__, "Invalid negative number of extra threads requested");
        }

        // read header
        this.header = bam_hdr_dup(header);
        sam_hdr_write(this.fp,this.header);
    }

    ~this()
    {
        debug (dhtslib_debug)
        {
            writeln("SAMWriter dtor");
        }

        bam_hdr_destroy(this.header);

    }

    /// close file
    void close(){
        const auto ret = hts_close(this.fp);
        if (ret < 0)
            stderr.writefln("There was an error closing %s", fromStringz(this.fn));
    }

    /// write out to disk
    void write(SAMRecord * rec){
        const auto ret = sam_write1(this.fp,this.header,rec.b);
        assert(ret>=0);
    }
}
debug(dhtslib_unittest)
unittest{
    writeln();
    import dhtslib.sam;
    import htslib.hts_log : hts_log_info;
    import std.path:buildPath,dirName;
    import std.string:fromStringz;
    hts_set_log_level(htsLogLevel.HTS_LOG_TRACE);
    hts_log_info(__FUNCTION__, "Testing SAMWriter");
    hts_log_info(__FUNCTION__, "Loading test file");
    auto sam = SAMFile(buildPath(dirName(dirName(dirName(dirName(__FILE__)))),"htslib","test","auxf#values.sam"), 0);
    auto readrange = sam.allRecords;
    auto sam2 = SAMWriter("test.bam",sam.header);
    hts_log_info(__FUNCTION__, "Getting read 1");
    auto read = readrange.front();
    sam2.write(&read);
    sam2.close;
    destroy(sam2);
    sam = SAMFile("test.bam");
    readrange = sam.allRecords;
    read = readrange.front();
    writeln(read.sequence);
    assert(read.sequence=="GCTAGCTCAG");
    // destroy(sam2);
}

/// Used in sorting
bool cmpInterval(hts_pair32_t a, hts_pair32_t b)
{
    if (a.beg < b.beg)
    {
        return true;
    }
    if (a.end < b.end)
    {
        return true;
    }
    return false;
}

/// Used in sorting
bool cmpRegList(hts_reglist_t a, hts_reglist_t b)
{
    if (a.tid < b.tid)
    {
        return true;
    }
    return false;
}

/// Parse text line of SAM; Used in unittest
private int parseSam(string line, bam_hdr_t* header, bam1_t* b)
{
    import htslib.kstring : kstring_t;
    import std.utf : toUTFz;

    kstring_t k;
    k.s = toUTFz!(char*)(line.dup);
    k.m = line.length + 1;
    k.l = line.length + 1;
    return sam_parse1(&k, header, b);
}
debug(dhtslib_unittest)
unittest{
    writeln();
    import dhtslib.sam;
    import htslib.hts_log : hts_log_info;
    import std.path:buildPath,dirName;
    import std.string:fromStringz;
    hts_set_log_level(htsLogLevel.HTS_LOG_TRACE);
    hts_log_info(__FUNCTION__, "Testing SAMFile & SAMRecord");
    hts_log_info(__FUNCTION__, "Loading test file");
    auto sam = SAMFile(buildPath(dirName(dirName(dirName(dirName(__FILE__)))),"htslib","test","auxf#values.sam"), 0);
    auto readrange = sam.allRecords;
    hts_log_info(__FUNCTION__, "Getting read 1");
    auto read = readrange.front();
    writeln(read.sequence);
    assert(read.sequence=="GCTAGCTCAG");
}

debug(dhtslib_unittest)
unittest
{
    import std.stdio : writeln;
    import std.range : drop;
    import std.utf : toUTFz;
    import htslib.hts_log; // @suppress(dscanner.suspicious.local_imports)
    import std.path:buildPath,dirName;
    import std.conv:to;
    hts_set_log_level(htsLogLevel.HTS_LOG_TRACE);
    hts_log_info(__FUNCTION__, "Loading sam file");
    auto range = File(buildPath(dirName(dirName(dirName(dirName(__FILE__)))),"htslib","test","realn01_exp-a.sam")).byLineCopy();
    auto b = bam_init1();
    auto hdr = bam_hdr_init();
    string hdr_str;
    for (auto i = 0; i < 4; i++)
    {
        hdr_str ~= range.front ~ "\n";
        range.popFront;
    }
    hts_log_info(__FUNCTION__, "Header");
    writeln(hdr_str);
    hdr = sam_hdr_parse(cast(int) hdr_str.length, toUTFz!(char*)(hdr_str));
    hts_log_info(__FUNCTION__, "Read status:" ~ parseSam(range.front, hdr, b).to!string);
    auto r = new SAMRecord(b);
    hts_log_info(__FUNCTION__, "Cigar" ~ r.cigar.toString);
    assert(r.cigar.toString == "6M1D117M5D28M");
}
