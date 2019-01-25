module dhtslib.tagvalue;

import std.stdio;
import std.string:fromStringz;
import dhtslib.htslib.sam:bam_aux_get,bam1_t,bam_aux2i;
import dhtslib.htslib.hts_log;

/**
This represents a tag value from a bam record
This should be to the bam specification.
It stores only a pointer to the tag and from there
can be parsed into any of the tag types but only if
the tag matches that type.
c byte
C ubyte
s short
S ushort
i int
I uint
f float
B? array of type ?
Z char array

Memory layout
8/9 are example values
2 is a count of the array
the ubyte * starts at the type char
c | 8|
s |  | 8|
i |  |  |  | 8|
B |i |  |  |  | 2|  |  |  | 8|  |  |  | 9|

*/

struct TagValue{
    ubyte* data;
    this(bam1_t * b,char[2] tag){
        data=bam_aux_get(b,tag);
        if(data==null)
            hts_log_warning(__FUNCTION__,(tag~" doesn't exist for this record").idup);
    }

    string to(T:string)(){
        return fromStringz(cast(char*)&data[1]).idup;
    }

    T to(T)(){
        return *cast(T*)data[1..T.sizeof+1].ptr;
    }

    T[] to(T:T[])(){
        int n=*cast(int*)data[2..6].ptr;
        return (cast(T*)(data[6..T.sizeof+6].ptr))[0..n];
    }

    //bool check(T){
    //    return
    //}

    string toString(){
        if(data !is null &&cast(char)data[0]=='Z'){
            return fromStringz(cast(char*)&data[1]).idup;
        }
        return "";
    }
}

unittest{
    TagValue v;
    ubyte[12] testdata;
    testdata[0]=cast(ubyte)'B';
    testdata[1]=cast(ubyte)'C';
    *cast(int*)testdata[2..6].ptr= 3;
    testdata[6]=1;
    testdata[8]=2;
    testdata[10]=3;
    v.data=testdata.ptr;
    writeln("testing array");
    assert(v.to!(ushort[])==[1,2,3]);
    ubyte[5] testdata2;
    testdata2[0]=cast(ubyte)'i';
    *cast(int*)testdata2[1..5].ptr= 3;
    v.data=testdata2.ptr;
    writeln("testing int");
    assert(v.to!int==3);
}
unittest{
    import dhtslib.sam;
    import dhtslib.htslib.hts_log;
    hts_set_log_level(htsLogLevel.HTS_LOG_TRACE);
    hts_log_info(__FUNCTION__,"Testing tagvalue");
    hts_log_info(__FUNCTION__,"Loading test file");
    auto bam=SAMFile("htslib/test/auxf#values.sam",0);
    hts_log_info(__FUNCTION__,"Getting read 1");
    auto readrange=bam.all_records();
    auto read=readrange.front;
    hts_log_info(__FUNCTION__,"Testing string");
    assert(read["RG"].to!string=="ID");
    hts_log_info(__FUNCTION__,"Testing char");
    assert(read["A!"].to!char=='!');
    assert(read["Ac"].to!char=='c');
    assert(read["AC"].to!char=='C');
    hts_log_info(__FUNCTION__,"Testing int");
    assert(read["I0"].to!byte==0);
    assert(read["I1"].to!byte==1);
    assert(read["I2"].to!byte==127);
    assert(read["I3"].to!ubyte==128);
    assert(read["I4"].to!ubyte==255);
    assert(read["I5"].to!short==256);
    assert(read["I6"].to!short==32767);
    assert(read["I7"].to!ushort==32768);
    assert(read["I8"].to!ushort==65535);
    assert(read["I9"].to!int==65536);
    assert(read["IA"].to!int==2147483647);
    assert(read["i1"].to!byte==-1);
    assert(read["i2"].to!byte==-127);
    assert(read["i3"].to!byte==-128);
    assert(read["i4"].to!short==-255);
    assert(read["i5"].to!short==-256);
    assert(read["i6"].to!short==-32767);
    assert(read["i7"].to!short==-32768);
    assert(read["i8"].to!int==-65535);
    assert(read["i9"].to!int==-65536);
    assert(read["iA"].to!int==-2147483647);
    assert(read["iB"].to!int==-2147483648);
    hts_log_info(__FUNCTION__,"Testing float");
    assert(read["F0"].to!float==-1.0);
    assert(read["F1"].to!float==0.0);
    assert(read["F2"].to!float==1.0);
    //assert(read["F3"].to!float==9.9e-19);
    //assert(read["F4"].to!float==-9.9e-19);
    readrange.popFront;
    read=readrange.front;
    hts_log_info(__FUNCTION__,"Testing arrays");
    assert(read["Bs"].to!(short[])==[-32768,-32767,0,32767]);
    assert(read["Bi"].to!(int[])==[-2147483648,-2147483647,0,2147483647]);
    assert(read["BS"].to!(ushort[])==[0,32767,32768,65535]);
    assert(read["BI"].to!(uint[])==[0,2147483647,2147483648,4294967295]);
}
