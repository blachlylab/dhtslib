/**
Module provides a parser for SAM/BAM record auxillary tags.

Reference: https://samtools.github.io/hts-specs/SAMtags.pdf
*/
module dhtslib.sam.tagvalue;

import std.stdio;
import std.meta : AliasSeq, staticIndexOf;
import std.string : fromStringz;
import htslib.sam : bam_aux_get, bam1_t, bam_aux2i;
import htslib.hts_log;
import std.conv : to;
import std.exception : enforce, assertThrown;
import std.math : approxEqual;

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
H   Byte array in the Hex format (network byte order / big-endian) //unknown if still supported
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
    private int refct = 1;      // Postblit refcounting in case the object is passed around

    private ubyte* data;

    /** Constructor

    Usage: auto t = TagValue(b, 'XX') where b is bam1_t* BAM record and XX is tag
    */
    this(bam1_t* b, char[2] tag)
    {
        data = bam_aux_get(b, tag);
    }

    this(this)
    {
        refct++;
    }

    ~this()
    {
        --refct;
    }

    invariant(){
        assert(refct >= 0);
    }

    auto references(){
        return refct;
    }

    /// check if empty/exists/null
    @property
    bool exists()
    {
        return this.data is null ? false : true;
    }

    /* Tag type checking */

    /// Check if tag type is type T
    bool check(T)()
    {
        enforce(this.exists,"Tag doesn't exist");
        return TypeChars[TypeIndex!T] == cast(char) data[0];
    }
    /// Check if tag type is type T
    bool check(T : string)()
    {
        enforce(this.exists,"Tag doesn't exist");
        return TypeChars[TypeIndex!T] == cast(char) data[0];
    }
    /// Check if tag type is type T
    bool check(T : T[])()
    {
        enforce(this.exists,"Tag doesn't exist");
        return (cast(char) data[0] == 'B') && (TypeChars[TypeIndex!T] == cast(char) data[1]);
    }

    /// Check if tag type is type T
    bool checkArray()
    {
        enforce(this.exists,"Tag doesn't exist");
        return cast(char) data[0] == 'B';
    }

    /// Check if tag type is type T
    bool checkHexByteArray()
    {
        enforce(this.exists,"Tag doesn't exist");
        return cast(char) data[0] == 'H';
    }

    /* Tag conversion */

    /// Convert tag value to D string
    string to(T : string)()
    {
        enforce(this.check!string || this.checkHexByteArray,"Tag is not type Z or H");
        return fromStringz(cast(char*)&data[1]).idup;
    }
    /// Convert tag value to D type
    T to(T)()
    {
        enforce(this.check!T,"Tag is not type " ~ T.stringof);
        return *cast(T*) data[1 .. T.sizeof + 1].ptr;
    }
    /// Convert array tag value D array
    T[] to(T : T[])()
    {
        enforce(this.check!(T[]),"Tag is not type " ~ T.stringof);
        int n = *cast(int*) data[2 .. 6].ptr;
        return (cast(T*)(data[6 .. T.sizeof + 6].ptr))[0 .. n];
    }
    

    /// Convert any tag value to string
    string toString()
    {
        enforce(this.exists,"Tag doesn't exist");
        switch (cast(char) data[0])
        {
        case 'c':
            return to!byte.to!string;
        case 'C':
            return to!ubyte.to!string;
        case 's':
            return to!short.to!string;
        case 'S':
            return to!ushort.to!string;
        case 'i':
            return to!int.to!string;
        case 'I':
            return to!uint.to!string;
        case 'f':
            return to!float.to!string;
        case 'Z':
        case 'H':
            return to!string;
        case 'B':
            switch (cast(char) data[1])
            {
            case 'c':
                return to!(byte[]).to!string;
            case 'C':
                return to!(ubyte[]).to!string;
            case 's':
                return to!(short[]).to!string;
            case 'S':
                return to!(ushort[]).to!string;
            case 'i':
                return to!(int[]).to!string;
            case 'I':
                return to!(uint[]).to!string;
            case 'f':
                return to!(float[]).to!string;
            default:
                throw new Exception("Array Tag malformed");    
            }
        default:
            throw new Exception("Tag malformed");
        }
    }
    /// Convert tag value to integer
    long toInt()
    {
        enforce(this.exists,"Tag doesn't exist");
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
            throw new Exception("Tag is not numeric or is malformed");
        }
    }
    /// Convert tag value to integer array
    long[] toIntArray()
    {
        enforce(this.exists,"Tag doesn't exist");
        enforce(this.checkArray,"Tag is not a numeric array");
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
            throw new Exception("Tag is malformed");
        }
    }
    /// Convert tag value to float array
    float[] toFloatArray()
    {
        enforce(this.exists,"Tag doesn't exist");
        enforce(this.checkArray,"Tag is not an array");
        enforce(this.check!(float[]),"Tag is not a float array");
        return to!(float[]);
    }
}

debug (dhtslib_unittest) unittest
{
    TagValue v;
    assert(!v.exists);
    ubyte[12] testdata;
    assertThrown(v.toIntArray);
    assertThrown(v.toInt);
    assertThrown(v.toString);
    testdata[0] = cast(ubyte) 'B';
    testdata[1] = cast(ubyte) 'S';
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
    auto bam = SAMFile(buildPath(dirName(dirName(dirName(dirName(__FILE__)))), "htslib",
            "test", "auxf#values.sam"), 0);

    hts_log_info(__FUNCTION__, "Getting read 1");
    auto readrange = bam.allRecords(); // @suppress(dscanner.suspicious.unmodified)
    assert(readrange.empty == false);
    auto read = readrange.front;
    assert(read["RG"].references == 2);

    hts_log_info(__FUNCTION__, "Testing string");
    assert(read["RG"].to!string == "ID");

    hts_log_info(__FUNCTION__, "Testing char");
    assert(read["A!"].to!char == '!');
    assert(read["Ac"].to!char == 'c');
    assert(read["AC"].to!char == 'C');

    hts_log_info(__FUNCTION__, "Testing integral checks");
    assert(read["I0"].check!ubyte);
    assert(read["I1"].check!ubyte);
    assert(read["I2"].check!ubyte);
    assert(read["I3"].check!ubyte);
    assert(read["I4"].check!ubyte);
    assert(read["I5"].check!ushort);
    assert(read["I6"].check!ushort);
    assert(read["I7"].check!ushort);
    assert(read["I8"].check!ushort);
    assert(read["I9"].check!uint);
    assert(read["IA"].check!uint);
    assert(read["i1"].check!byte);
    assert(read["i2"].check!byte);
    assert(read["i3"].check!byte);
    assert(read["i4"].check!short);
    assert(read["i5"].check!short);
    assert(read["i6"].check!short);
    assert(read["i7"].check!short);
    assert(read["i8"].check!int);
    assert(read["i9"].check!int);
    assert(read["iA"].check!int);
    assert(read["iB"].check!int);

    hts_log_info(__FUNCTION__, "Testing integral conversion");
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

    hts_log_info(__FUNCTION__, "Testing integral toString");
    assert(read["I0"].toString == "0");
    assert(read["I1"].toString == "1");
    assert(read["I2"].toString == "127");
    assert(read["I3"].toString == "128");
    assert(read["I4"].toString == "255");
    assert(read["I5"].toString == "256");
    assert(read["I6"].toString == "32767");
    assert(read["I7"].toString == "32768");
    assert(read["I8"].toString == "65535");
    assert(read["I9"].toString == "65536");
    assert(read["IA"].toString == "2147483647");
    assert(read["i1"].toString == "-1");
    assert(read["i2"].toString == "-127");
    assert(read["i3"].toString == "-128");
    assert(read["i4"].toString == "-255");
    assert(read["i5"].toString == "-256");
    assert(read["i6"].toString == "-32767");
    assert(read["i7"].toString == "-32768");
    assert(read["i8"].toString == "-65535");
    assert(read["i9"].toString == "-65536");
    assert(read["iA"].toString == "-2147483647");
    assert(read["iB"].toString == "-2147483648");

    hts_log_info(__FUNCTION__, "Testing integral toInt");
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

    hts_log_info(__FUNCTION__, "Testing float checks");

    assert(read["F0"].check!float);
    assert(read["F1"].check!float);
    assert(read["F2"].check!float);

    hts_log_info(__FUNCTION__, "Testing float conversion");
    assert(read["F0"].to!float == -1.0);
    assert(read["F1"].to!float == 0.0);
    assert(read["F2"].to!float == 1.0);

    hts_log_info(__FUNCTION__, "Testing float toString");

    assert(approxEqual(read["F0"].toString.to!float, -1.0));
    assert(approxEqual(read["F1"].toString.to!float, 0.0));
    assert(approxEqual(read["F2"].toString.to!float, 1.0));

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
    
    hts_log_info(__FUNCTION__, "Testing array toString");
    assert(read["Bs"].toString == "[-32768, -32767, 0, 32767]");
    assert(read["Bi"].toString == "[-2147483648, -2147483647, 0, 2147483647]");
    assert(read["BS"].toString == "[0, 32767, 32768, 65535]");
    assert(read["BI"].toString == "[0, 2147483647, 2147483648, 4294967295]");

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

    hts_log_info(__FUNCTION__, "Testing float Array");
    float[] arr = [10.0,11.0,12.1];
    read["fA"] = arr;
    assert(read["fA"].to!(float[]) == arr);
    assert(read["fA"].toFloatArray == arr);
    assert(read["fA"].toString == "[10, 11, 12.1]");

    hts_log_info(__FUNCTION__, "Testing byte Array");
    byte[] arr2 = [10, -10]; 
    read["cA"] = arr2;
    assert(read["cA"].to!(byte[]) == arr2);
    assert(read["cA"].toIntArray == arr2.to!(long[]));
    assert(read["cA"].toString == "[10, -10]");

    hts_log_info(__FUNCTION__, "Testing ubyte Array");
    ubyte[] arr3 = [10, 11]; 
    read["CA"] = arr3;
    assert(read["CA"].to!(ubyte[]) == arr3);
    assert(read["CA"].toIntArray == arr3.to!(long[]));
    assert(read["CA"].toString == "[10, 11]");

}
