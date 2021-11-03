/**
This module simplifies working with CIGAR strings/ops from SAM/BAM/CRAM alignment records.
*/
module dhtslib.sam.cigar;

import std.stdio;
import std.bitmanip : bitfields;
import std.array : join;
import std.algorithm : map, each, sum;
import std.conv : to;
import std.range : array;
import std.traits : isIntegral;

import htslib.hts_log;
import dhtslib.sam.record : SAMRecord;
import dhtslib.coordinates;
import dhtslib.memory;
import htslib.sam : bam_get_cigar;

/// Represents a CIGAR string
/// https://samtools.github.io/hts-specs/SAMv1.pdf §1.4.6
struct Cigar
{
    /// array of distinct CIGAR ops 
    /// private now to force fixing some reference issues
    private CigarOp[] ops;
    private Bam1 b;
    /// Construct Cigar from SAMRecord
    this(Bam1 b)
    {
        this.b = b;
        auto cigar = bam_get_cigar(this.b);
        this(cigar, this.b.core.n_cigar);
    }

    /// Construct Cigar from raw data
    this(uint* cigar, uint length)
    {
        this((cast(CigarOp*) cigar)[0 .. length]);
    }

    /// Construct Cigar from an array of CIGAR ops
    this(CigarOp[] ops)
    {
        this.ops = ops;
    }

    /// null CIGAR -- don't read!
    bool is_null()
    {
        return this.ops.length == 0;
    }

    /// Format Cigar struct as CIGAR string in accordance with SAM spec
    string toString()
    {
        return ops.map!(x => x.length.to!string ~ CIGAR_STR[x.op]).array.join;
    }

    /// get length of cigar
    auto length(){
        return this.ops.length;
    }
    
    /// set length of cigar
    void length(T)(T len)
    if(isIntegral!T)
    {
        this.ops.length = len;
    }

    /// copy cigar
    auto const dup(){
        return Cigar(this.ops.dup);
    }

    /// return the alignment length expressed by this Cigar
    @property int refBasesCovered()
    {
        int len;
        foreach (op; this.ops)
        {
            // look ma no branching (ignore the above foreach)
            // add op.length if reference consuming
            // mask is 0 if not reference consuming
            // mask is -1 if reference consuming
            len += op.length & (uint.init - (((CIGAR_TYPE >> ((op.op & 0xF) << 1)) & 2) >> 1));
        }
        return len;
    }

    deprecated("Use camelCase names instead")
    alias ref_bases_covered = refBasesCovered;
    /// previous alignedLength function had a bug and 
    /// it is just a duplicate of ref_bases_covered
    alias alignedLength = refBasesCovered;

    /// allow foreach on Cigar
    int opApply(int delegate(ref CigarOp) dg) {
        int result = 0;

        foreach (op; ops) {
            result = dg(op);

            if (result) {
                break;
            }
        }

        return result;
    }
    
    /// Get a Slice of cigar ops from a range in the Cigar string
    auto opSlice()
    {
        return ops;
    }

    /// Get a Slice of cigar ops from a range in the Cigar string
    auto opSlice(T)(T start, T end) if (isIntegral!T)
    {
        return ops[start .. end];
    }

    /// Get a cigar op from a single position in the Cigar string
    ref auto opIndex(ulong i)
    {
        return ops[i];
    }

    /// Assign a cigar op at a single position in the Cigar string
    auto opIndexAssign(CigarOp value, size_t index)
    {
        return ops[index] = value;
    }

    /// Assign a range cigar ops over a range in the Cigar string
    auto opSliceAssign(T)(CigarOp[] values, T start, T end)
    {
        assert(end - start == values.length);
        return ops[start .. end] = values;
    }

    auto opAssign(T: CigarOp[])(T value)
    {
        this.ops = value;
        return this;
    }

    auto opAssign(T: Cigar)(T value)
    {
        this.ops = value.ops;
        return this;
    }

    size_t opDollar()
    {
        return length;
    }

    auto opBinary(string op)(const Cigar rhs) const
    if(op == "~")
    {
        return Cigar(this.dup[] ~ rhs.dup[]);
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

    mixin(bitfields!( //lower 4 bits store op
            Ops, "op", 4, //higher 28 bits store length
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

    string toString(){
        return this.length.to!string ~ CIGAR_STR[this.op];
    }
}

/// Credit to Biod for this code below
/// https://github.com/biod/BioD from their bam.cigar module
/// True iff operation is one of M, =, X, I, S
bool isQueryConsuming(T)(T op) nothrow @nogc
{
    return ((CIGAR_TYPE >> ((op & 0xF) * 2)) & 1) != 0;
}

/// True iff operation is one of M, =, X, D, N
bool isReferenceConsuming(T)(T op)nothrow @nogc
{
    return ((CIGAR_TYPE >> ((op & 0xF) * 2)) & 2) != 0;
}

/// True iff operation is one of M, =, X
bool isMatchOrMismatch(T)(T op) nothrow @nogc
{
    return ((CIGAR_TYPE >> ((op & 0xF) * 2)) & 3) == 3;
}

/// True iff operation is one of 'S', 'H'
bool isClipping(T)(T op)nothrow @nogc
{
    return ((op & 0xF) >> 1) == 2; // 4 or 5
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

debug (dhtslib_unittest) unittest
{
    writeln();
    import dhtslib.sam;
    import htslib.hts_log;
    import std.path : buildPath, dirName;

    hts_set_log_level(htsLogLevel.HTS_LOG_TRACE);
    hts_log_info(__FUNCTION__, "Testing cigar");
    hts_log_info(__FUNCTION__, "Loading test file");
    auto bam = SAMFile(buildPath(dirName(dirName(dirName(dirName(__FILE__)))), "htslib",
            "test", "range.bam"), 0);
    auto readrange = bam["CHROMOSOME_I", ZB(914)];
    hts_log_info(__FUNCTION__, "Getting read 1");
    assert(readrange.empty == false);
    auto read = readrange.front();
    auto cigar = read.cigar;
    writeln(read.queryName);
    hts_log_info(__FUNCTION__, "Cigar:" ~ cigar.toString());
    writeln(cigar.toString());
    assert(cigar.toString() == "78M1D22M");
    assert(cigar[0] == CigarOp(78, Ops.MATCH));
    assert(Cigar(cigar[0 .. 2]).toString == "78M1D");
    cigar[0] = CigarOp(4,Ops.HARD_CLIP);
    cigar[1..3] = [CigarOp(2, Ops.INS), CigarOp(21, Ops.MATCH)];
    assert(cigar[1 .. $] == [CigarOp(2, Ops.INS), CigarOp(21, Ops.MATCH)]);
    assert(cigar.toString() == "4H2I21M");
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

debug (dhtslib_unittest) unittest
{
    writeln();
    hts_log_info(__FUNCTION__, "Testing is_query_consuming and is_reference_consuming");
    string c = "130M2D40M";
    auto cig = cigarFromString(c);
    hts_log_info(__FUNCTION__, "Cigar:" ~ cig.toString());
    assert(cig.toString() == c);
    assert(cig.ops[0].op.isQueryConsuming && cig.ops[0].op.isReferenceConsuming);
    assert(!cig.ops[1].op.isQueryConsuming && cig.ops[1].op.isReferenceConsuming);
}

/// Range-based iteration of a Cigar string
/// Returns a range of Ops that is the length of all op lengths
/// e.g. if the Cigar is 4M5D2I4M
/// CigarItr will return a range of MMMMDDDDDIIMMMM
/// Range is of the Ops enum not chars
struct CigarItr
{
    Cigar cigar;
    CigarOp[] ops;
    CigarOp current;

    this(Cigar c)
    {
        // Copy the cigar
        cigar = c;
        ops = cigar[];
        current = ops[0];
        current.length = current.length - 1;
    }

    Ops front()
    {
        return current.op;
    }

    void popFront()
    {
        // import std.stdio;
        // writeln(current);
        // writeln(cigar);
        if(ops.length == 1 && current.length == 0)
            ops = [];
        else if (current.length == 0)
        {
            ops = ops[1 .. $];
            current = ops[0];
            current.length = current.length - 1;
        }
        else if (current.length != 0)
            current.length = current.length - 1;
    }

    bool empty()
    {
        return ops.length == 0;
    }
}

debug (dhtslib_unittest) unittest
{
    import std.algorithm : map;

    writeln();
    hts_log_info(__FUNCTION__, "Testing CigarItr");
    string c = "7M2D4M";
    auto cig = cigarFromString(c);
    hts_log_info(__FUNCTION__, "Cigar:" ~ cig.toString());
    auto itr = CigarItr(cig);
    assert(itr.map!(x => CIGAR_STR[x]).array.idup == "MMMMMMMDDMMMM");
}

/// Coordinate pair representing aligned query and reference positions,
/// and the CIGAR op at or effecting their alignment.
struct AlignedCoordinate
{
    Coordinate!(Basis.zero) qpos, rpos;
    Ops cigar_op;
}

/// Iterator yielding all AlignedCoordinate pairs for a given CIGAR string
struct AlignedCoordinatesItr
{
    CigarItr itr;
    AlignedCoordinate current;

    this(Cigar cigar)
    {
        itr = CigarItr(cigar);
        current.qpos = current.rpos = Coordinate!(Basis.zero)(-1);
        current.cigar_op = itr.front;
        current.qpos += ((CIGAR_TYPE >> ((current.cigar_op & 0xF) << 1)) & 1);
        current.rpos += (((CIGAR_TYPE >> ((current.cigar_op & 0xF) << 1)) & 2) >> 1);
    }

    /// InputRange interface
    AlignedCoordinate front()
    {
        return current;
    }

    /// InputRange interface
    void popFront()
    {
        itr.popFront;
        current.cigar_op = itr.front;
        current.qpos += ((CIGAR_TYPE >> ((current.cigar_op & 0xF) << 1)) & 1);
        current.rpos += (((CIGAR_TYPE >> ((current.cigar_op & 0xF) << 1)) & 2) >> 1);
    }

    /// InputRange interface
    bool empty()
    {
        return itr.empty;
    }
}

debug (dhtslib_unittest) unittest
{
    writeln();
    import dhtslib.sam;
    import htslib.hts_log;
    import std.path : buildPath, dirName;
    import std.algorithm : map;
    import std.array : array;

    hts_set_log_level(htsLogLevel.HTS_LOG_TRACE);
    hts_log_info(__FUNCTION__, "Testing cigar");
    hts_log_info(__FUNCTION__, "Loading test file");
    auto bam = SAMFile(buildPath(dirName(dirName(dirName(dirName(__FILE__)))), "htslib",
            "test", "range.bam"), 0);
    auto readrange = bam["CHROMOSOME_I", ZB(914)];
    hts_log_info(__FUNCTION__, "Getting read 1");
    assert(readrange.empty == false);
    auto read = readrange.front();
    auto coords = read.getAlignedCoordinates(77, 82);
    assert(coords.map!(x => x.rpos).array == [77,78,79,80,81]);
}
