/**
Module provides a parser for SAM/BAM record auxillary tags.

Reference: https://samtools.github.io/hts-specs/SAMtags.pdf
*/
module dhtslib.tagvalue;

import std.stdio;
import std.meta : AliasSeq, staticIndexOf;
import std.string : fromStringz;
import htslib.sam : bam_aux_get, bam1_t, bam_aux2i;
import htslib.hts_log;
import std.conv : to;

alias Types = AliasSeq!(byte, ubyte, short, ushort, int, uint, float, string, char);
enum TypeIndex(T) = staticIndexOf!(T, Types);
/// See https://samtools.github.io/hts-specs/SAMv1.pdf sec 1.5
char[9] TypeChars = ['c', 'C', 's', 'S', 'i', 'I', 'f', 'Z', 'A'];

/**

This represents a SAM/BAM record tag value, as outlined in the SAM specs §1.5.

The struct itself stores only a pointer to the tag, and has member functions
to parse into any of the tag types (but only if the tag matches that type) (TODO: is this true?)

Primary Types:
A   Printable character
i   Signed integer (see specs §1.5 footnote on size)
f   Single-precision float
Z   Printable string, including space
H   Byte array in the Hex format (network byte order / big-endian)
B   Integer or numeric array

Byte-array (B) types:
c   byte
C   ubyte
s   short
S   ushort
i   int32
I   uint32
f   float (spec does not indicate precision)

Memory layout
pipes delimit byte boundaries in an array
8/9 are example values
2 is a count of the array
the ubyte * starts at the type char
c | 8|
s |  | 8|
i |  |  |  | 8|
B |i |  |  |  | 2|  |  |  | 8|  |  |  | 9|


Alias seq allows us to have an enum of types.
https://forum.dlang.org/post/kmdjfzpugudmwfrdgson@forum.dlang.org
Thanks Paul!

Usage: auto t = TagValue(b, 'XX') where b is bam1_t* BAM record and XX is tag
*/
struct TagValue
{
    private ubyte* data;

    /** Constructor

    Usage: auto t = TagValue(b, 'XX') where b is bam1_t* BAM record and XX is tag
    */
    this(bam1_t* b, char[2] tag)
    {
        data = bam_aux_get(b, tag);
        debug
        {
            if (data is null)
                hts_log_warning(__FUNCTION__, (tag ~ " doesn't exist for this record").idup);
        }
    }

    /// check if empty/exists/null
    @property
    bool exists()
    {
        if (this.data is null) return false;
        return true;
    }

    /// Convert tag value
    string to(T : string)()
    {
        assert(this.data !is null);
        return fromStringz(cast(char*)&data[1]).idup;
    }
    /// Convert tag value
    T to(T)()
    {
        assert(this.data !is null);
        return *cast(T*) data[1 .. T.sizeof + 1].ptr;
    }
    /// Convert tag value
    T[] to(T : T[])()
    {
        assert(this.data !is null);
        int n = *cast(int*) data[2 .. 6].ptr;
        return (cast(T*)(data[6 .. T.sizeof + 6].ptr))[0 .. n];
    }
    /// Check if tag type is type T
    bool check(T)()
    {
        assert(this.data !is null);
        return TypeChars[TypeIndex!T] == cast(char) data[0];
    }
    /// Check if tag type is type T
    bool check(T : string)()
    {
        assert(this.data !is null);
        return TypeChars[TypeIndex!T] == cast(char) data[0];
    }
    /// Check if tag type is type T
    bool check(T : T[])()
    {
        assert(this.data !is null);
        return (cast(char) data[0] == 'B') && (TypeChars[TypeIndex!T] == cast(char) data[1]);
    }
    /// Convert tag value to string
    string toString() const
    {
        if (data !is null && cast(char) data[0] == 'Z')
        {
            return fromStringz(cast(char*)&data[1]).idup;
        }
        return "";
    }
    /// Convert tag value to integer
    long toInt()
    {
        assert(this.data !is null);
        switch (cast(char) data[0])
        {
        case 'c':
            return cast(long)(to!byte);
        case 'C':
            return cast(long)(to!ubyte);
        case 's':
            return cast(long)(to!short);
        case 'S':
            return cast(long)(to!ushort);
        case 'i':
            return cast(long)(to!int);
        case 'I':
            return cast(long)(to!uint);
        default:
            return long.min;
        }
    }
    /// Convert tag value to integer array
    long[] toIntArray()
    {
        assert(this.data !is null);
        switch (cast(char) data[1])
        {
        case 'c':
            return (to!(byte[]).to!(long[]));
        case 'C':
            return (to!(ubyte[]).to!(long[]));
        case 's':
            return (to!(short[]).to!(long[]));
        case 'S':
            return (to!(ushort[]).to!(long[]));
        case 'i':
            return (to!(int[]).to!(long[]));
        case 'I':
            return (to!(uint[]).to!(long[]));
        default:
            return [];
        }
    }
    /// Convert tag value to float array
    float[] toFloatArray()
    {
        assert(this.data !is null);
        return to!(float[]);
    }
}

debug (dhtslib_unittest) unittest
{
    TagValue v;
    ubyte[12] testdata;
    testdata[0] = cast(ubyte) 'B';
    testdata[1] = cast(ubyte) 'C';
    *cast(int*) testdata[2 .. 6].ptr = 3;
    testdata[6] = 1;
    testdata[8] = 2;
    testdata[10] = 3;
    v.data = testdata.ptr;
    writeln("testing array");
    assert(v.to!(ushort[]) == [1, 2, 3]);
    ubyte[5] testdata2;
    testdata2[0] = cast(ubyte) 'i';
    *cast(int*) testdata2[1 .. 5].ptr = 3;
    v.data = testdata2.ptr;
    writeln("testing int");
    assert(v.to!int == 3);
}

debug (dhtslib_unittest) unittest
{
    import dhtslib.sam; // @suppress(dscanner.suspicious.local_imports)
    import htslib.hts_log : hts_log_info;
    import std.path : buildPath, dirName;

    hts_set_log_level(htsLogLevel.HTS_LOG_TRACE);
    hts_log_info(__FUNCTION__, "Testing tagvalue");
    hts_log_info(__FUNCTION__, "Loading test file");
    auto bam = SAMFile(buildPath(dirName(dirName(dirName(__FILE__))), "htslib",
            "test", "auxf#values.sam"), 0);
    hts_log_info(__FUNCTION__, "Getting read 1");
    auto readrange = bam.all_records(); // @suppress(dscanner.suspicious.unmodified)
    auto read = readrange.front;
    hts_log_info(__FUNCTION__, "Testing string");
    assert(read["RG"].to!string == "ID");
    hts_log_info(__FUNCTION__, "Testing char");
    assert(read["A!"].to!char == '!');
    assert(read["Ac"].to!char == 'c');
    assert(read["AC"].to!char == 'C');
    hts_log_info(__FUNCTION__, "Testing int");
    assert(read["I0"].to!ubyte == 0);
    assert(read["I1"].to!ubyte == 1);
    assert(read["I2"].to!ubyte == 127);
    assert(read["I3"].to!ubyte == 128);
    assert(read["I4"].to!ubyte == 255);
    assert(read["I5"].to!ushort == 256);
    assert(read["I6"].to!ushort == 32_767);
    assert(read["I7"].to!ushort == 32_768);
    assert(read["I8"].to!ushort == 65_535);
    assert(read["I9"].to!uint == 65_536);
    assert(read["IA"].to!uint == 2_147_483_647);
    assert(read["i1"].to!byte == -1);
    assert(read["i2"].to!byte == -127);
    assert(read["i3"].to!byte == -128);
    assert(read["i4"].to!short == -255);
    assert(read["i5"].to!short == -256);
    assert(read["i6"].to!short == -32_767);
    assert(read["i7"].to!short == -32_768);
    assert(read["i8"].to!int == -65_535);
    assert(read["i9"].to!int == -65_536);
    assert(read["iA"].to!int == -2_147_483_647);
    assert(read["iB"].to!int == -2_147_483_648);
    assert(read["I0"].toInt == 0);
    assert(read["I1"].toInt == 1);
    assert(read["I2"].toInt == 127);
    assert(read["I3"].toInt == 128);
    assert(read["I4"].toInt == 255);
    assert(read["I5"].toInt == 256);
    assert(read["I6"].toInt == 32_767);
    assert(read["I7"].toInt == 32_768);
    assert(read["I8"].toInt == 65_535);
    assert(read["I9"].toInt == 65_536);
    assert(read["IA"].toInt == 2_147_483_647);
    assert(read["i1"].toInt == -1);
    assert(read["i2"].toInt == -127);
    assert(read["i3"].toInt == -128);
    assert(read["i4"].toInt == -255);
    assert(read["i5"].toInt == -256);
    assert(read["i6"].toInt == -32_767);
    assert(read["i7"].toInt == -32_768);
    assert(read["i8"].toInt == -65_535);
    assert(read["i9"].toInt == -65_536);
    assert(read["iA"].toInt == -2_147_483_647);
    assert(read["iB"].toInt == -2_147_483_648);
    hts_log_info(__FUNCTION__, "Testing float");
    assert(read["F0"].to!float == -1.0);
    assert(read["F1"].to!float == 0.0);
    assert(read["F2"].to!float == 1.0);
    hts_log_info(__FUNCTION__, "Running tag checking");
    assert(read["I0"].check!ubyte == true);
    assert(read["I5"].check!ushort == true);
    assert(read["I9"].check!uint == true);
    assert(read["i1"].check!byte == true);
    assert(read["i4"].check!short == true);
    assert(read["i8"].check!int == true);
    assert(read["F0"].check!float == true);
    readrange.popFront;
    read = readrange.front;
    hts_log_info(__FUNCTION__, "Testing arrays");
    assert(read["Bs"].to!(short[]) == [-32_768, -32_767, 0, 32_767]);
    assert(read["Bi"].to!(int[]) == [
            -2_147_483_648, -2_147_483_647, 0, 2_147_483_647
            ]);
    assert(read["BS"].to!(ushort[]) == [0, 32_767, 32_768, 65_535]);
    assert(read["BI"].to!(uint[]) == [
            0, 2_147_483_647, 2_147_483_648, 4_294_967_295
            ]);
    writeln(read["Bs"].toIntArray);
    assert(read["Bs"].toIntArray == [-32_768, -32_767, 0, 32_767]);
    assert(read["Bi"].toIntArray == [
            -2_147_483_648, -2_147_483_647, 0, 2_147_483_647
            ]);
    assert(read["BS"].toIntArray == [0, 32_767, 32_768, 65_535]);
    assert(read["BI"].toIntArray == [
            0, 2_147_483_647, 2_147_483_648, 4_294_967_295
            ]);
    hts_log_info(__FUNCTION__, "Running tag checking");
    assert(read["Bs"].check!(short[]) == true);
    assert(read["Bi"].check!(int[]) == true);
    assert(read["BS"].check!(ushort[]) == true);
    assert(read["BI"].check!(uint[]) == true);
}
