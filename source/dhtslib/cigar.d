/**
    This module simplifies working with CIGAR strings/ops in SAM/BAM/CRAM format.
*/
module dhtslib.cigar;

import std.stdio;
import std.bitmanip : bitfields;
import std.array : join;
import std.algorithm : map;
import std.algorithm.iteration : each;
import std.conv : to;
import std.range : array;

import dhtslib.htslib.hts_log;

/// Represents a CIGAR string
/// https://samtools.github.io/hts-specs/SAMv1.pdf §1.4.6
struct Cigar
{
    /// array of distinct CIGAR ops 
    CigarOp[] ops;

    /// Construct Cigar from raw data
    this(uint* cigar, int length)
    {
        ops=(cast(CigarOp *)cigar)[0..length];
    }

    /// Construct Cigar from an array of CIGAR ops
    this(CigarOp[] ops)
    {
        this.ops = ops;
    }

    bool is_null(){
        return ops.length == 1 && ops[0].raw == '*';
    }

    /// Format Cigar struct as CIGAR string in accordance with SAM spec
    string toString()
    {
        return ops.map!(x => x.length.to!string ~ CIGAR_STR[x.op]).array.join;
    }
    @property int ref_bases_covered(){
        int len;
        foreach(op;this.ops){
            switch (op.op)
            {
            case Ops.MATCH:
            case Ops.EQUAL:
            case Ops.DIFF:
            case Ops.DEL:
            case Ops.REF_SKIP:
                len += op.length;
                break;
            default:
                break;
            }
        }
        return len;
    }
    /// return the alignment length expressed by this Cigar
    int alignedLength()
    {
        int length = 0;
        foreach (op; ops)
        {
            switch (op.op)
            {
            default:
                {
                    break;
                }

            case Ops.MATCH:
            case Ops.EQUAL:
            case Ops.DIFF:
            case Ops.DEL:
            case Ops.INS:
                {
                    length += op.length;
                    break;
                }
            }
        }
        return length;
    }
}

// Each pair of bits has first bit set iff the operation is query consuming,
// and second bit set iff it is reference consuming.
//                                            X  =  P  H  S  N  D  I  M
private static immutable uint CIGAR_TYPE = 0b11_11_00_00_01_10_10_01_11;


/// Represents a distinct cigar operation
union CigarOp
{
    /// raw opcode
    uint raw;

    mixin(bitfields!(//lower 4 bits store op
            Ops, "op", 4,//higher 28 bits store length
            uint, "length", 28));

    /// construct Op from raw opcode
    this(uint raw)
    {
        this.raw = raw;
    }

    /// construct Op from an operator and operand (length)
    this(uint len, Ops op)
    {
        this.op = op;
        this.length = len;
    }
    /// Credit to Biod for this code below
    /// https://github.com/biod/BioD from their bam.cigar module
    /// True iff operation is one of M, =, X, I, S
    bool is_query_consuming() @property const nothrow @nogc {
        return ((CIGAR_TYPE >> ((raw & 0xF) * 2)) & 1) != 0;
    }

    /// True iff operation is one of M, =, X, D, N
    bool is_reference_consuming() @property const nothrow @nogc {
        return ((CIGAR_TYPE >> ((raw & 0xF) * 2)) & 2) != 0;
    }

    /// True iff operation is one of M, =, X
    bool is_match_or_mismatch() @property const nothrow @nogc {
        return ((CIGAR_TYPE >> ((raw & 0xF) * 2)) & 3) == 3;
    }

    /// True iff operation is one of 'S', 'H'
    bool is_clipping() @property const nothrow @nogc {
        return ((raw & 0xF) >> 1) == 2; // 4 or 5
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
string CIGAR_STR = "MIDNSHP=XB";
/// ditto
enum Ops
{
    MATCH = 0,
    INS = 1,
    DEL = 2,
    REF_SKIP = 3,
    SOFT_CLIP = 4,
    HARD_CLIP = 5,
    PAD = 6,
    EQUAL = 7,
    DIFF = 8,
    BACK = 9
}
debug(dhtslib_unittest)
unittest
{
    writeln();
    import dhtslib.sam;
    import dhtslib.htslib.hts_log;

    hts_set_log_level(htsLogLevel.HTS_LOG_TRACE);
    hts_log_info(__FUNCTION__, "Testing cigar");
    hts_log_info(__FUNCTION__, "Loading test file");
    auto bam = SAMFile("htslib/test/range.bam", 0);
    auto readrange = bam["CHROMOSOME_I", 914];
    hts_log_info(__FUNCTION__, "Getting read 1");
    auto read = readrange.front();
    writeln(read.queryName);
    hts_log_info(__FUNCTION__, "Cigar:" ~ read.cigar.toString());
    assert(read.cigar.toString() == "78M1D22M");
}

/// return Cigar struct for a given CIGAR string (e.g. from SAM line)
Cigar cigarFromString(string cigar)
{
    import std.regex;

    return Cigar(match(cigar, regex(`(\d+)([A-Z=])`, "g")).map!(m => CigarOp(m[1].to!uint,
            m[2].to!char.charToOp)).array);
}

/// Convert single char representing a CIGAR Op to Op union
Ops charToOp(char c)
{
    foreach (i, o; CIGAR_STR)
    {
        if (c == o)
        {
            return cast(Ops) i;
        }
    }
    return cast(Ops) 9;
}
debug(dhtslib_unittest)
unittest
{
    writeln();
    hts_log_info(__FUNCTION__, "Testing string2cigar");
    hts_log_info(__FUNCTION__, "Loading test file");
    string c = "130M2D40M";
    auto cig = cigarFromString(c);
    hts_log_info(__FUNCTION__, "Cigar:" ~ cig.toString());
    assert(cig.toString() == c);
    assert(cig.ops[0].is_query_consuming && cig.ops[0].is_reference_consuming);
    assert(!cig.ops[1].is_query_consuming && cig.ops[1].is_reference_consuming);
}
