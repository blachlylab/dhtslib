module dhtslib.cigar;

import std.stdio;
import std.bitmanip:bitfields;
import std.array:join;
import std.algorithm:map;
import std.algorithm.iteration:each;
import std.conv:to;
import std.range:array;

import dhtslib.htslib.hts_log;


/**
Represents a cigar string
*/
struct Cigar
{
    CigarOp[] ops;

    this(uint * cigar, int length)
    {
        hts_log_debug(__FUNCTION__,"Cigar Length:"~length.to!string);
        cigar[0..length].each!(x=>ops~=CigarOp(x));
    }

    string toString()
    {
        return ops.map!(x=>x.length.to!string~CIGAR_STR[x.op]).array.join;
    }
}

/**
Represents a cigar operation
*/
union CigarOp
{
    uint raw;

    mixin(bitfields!(
    //lower 4 bits store op
    Ops,"op",4,
    //higher 28 bits store length
    uint,"length",28
    ));

    this(uint raw)
    {
        this.raw=raw;
    }
}
/**
Represents all ops
#define BAM_CIGAR_STR   "MIDNSHP=XB"
#define BAM_CMATCH      0
#define BAM_CINS        1
#define BAM_CDEL        2
#define BAM_CREF_SKIP   3
#define BAM_CSOFT_CLIP  4
#define BAM_CHARD_CLIP  5
#define BAM_CPAD        6
#define BAM_CEQUAL      7
#define BAM_CDIFF       8
#define BAM_CBACK       9

*/
string CIGAR_STR="MIDNSHP=XB";
enum Ops {
    MATCH    =  0,
    INS      =  1,
    DEL      =  2,
    REF_SKIP =  3,
    SOFT_CLIP=  4,
    HARD_CLIP=  5,
    PAD      =  6,
    EQUAL    =  7,
    DIFF     =  8,
    BACK     =  9
}


unittest{
    import dhtslib.sam;
    import dhtslib.htslib.hts_log;
    hts_set_log_level(htsLogLevel.HTS_LOG_TRACE);
    hts_log_info(__FUNCTION__,"Testing cigar");
    hts_log_info(__FUNCTION__,"Loading test file");
    auto bam=SAMFile("htslib/test/range.bam",0);
    auto readrange=bam["CHROMOSOME_I",914];
    hts_log_info(__FUNCTION__,"Getting read 1");
    auto read=readrange.front();
    writeln(read.queryName);
    hts_log_info(__FUNCTION__,"Cigar:"~read.cigar.toString());
    assert(read.cigar.toString()=="78M1D22M");
}


