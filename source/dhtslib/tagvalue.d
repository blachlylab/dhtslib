module dhtslib.tagvalue;

import std.stdio;
import std.meta:AliasSeq,staticIndexOf;
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
Bc array of type byte
Z char array
H hex?

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
*/

alias Types = AliasSeq!(byte,ubyte,short,ushort,int,uint,float,string,char);
enum TypeIndex(T)=staticIndexOf!(T,Types);
char[9] TypeChars = ['c','C','s','S','i','I','f','Z','A'];

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

    bool check(T)(){
        return TypeChars[TypeIndex!T] == cast(char)data[0];
    }
    bool check(T:string)(){
        return TypeChars[TypeIndex!T] == cast(char)data[0];
    }

    bool check(T:T[])(){
        return (cast(char)data[0]=='B') && (TypeChars[TypeIndex!T] == cast(char)data[1]);
    }

    string toString(){
        if(data !is null &&cast(char)data[0]=='Z'){
            return fromStringz(cast(char*)&data[1]).idup;
        }
        return "";
    }
}

debug(dhtslib_unittest)
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
debug(dhtslib_unittest)
unittest{
    import dhtslib.sam;
    import dhtslib.htslib.hts_log;
    import std.path:buildPath,dirName;
    hts_set_log_level(htsLogLevel.HTS_LOG_TRACE);
    hts_log_info(__FUNCTION__,"Testing tagvalue");
    hts_log_info(__FUNCTION__,"Loading test file");
    auto bam=SAMFile(buildPath(dirName(dirName(dirName(__FILE__))),"htslib","test","auxf#values.sam"),0);
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
    assert(read["I0"].to!ubyte==0);
    assert(read["I1"].to!ubyte==1);
    assert(read["I2"].to!ubyte==127);
    assert(read["I3"].to!ubyte==128);
    assert(read["I4"].to!ubyte==255);
    assert(read["I5"].to!ushort==256);
    assert(read["I6"].to!ushort==32_767);
    assert(read["I7"].to!ushort==32_768);
    assert(read["I8"].to!ushort==65_535);
    assert(read["I9"].to!uint==65_536);
    assert(read["IA"].to!uint==2_147_483_647);
    assert(read["i1"].to!byte==-1);
    assert(read["i2"].to!byte==-127);
    assert(read["i3"].to!byte==-128);
    assert(read["i4"].to!short==-255);
    assert(read["i5"].to!short==-256);
    assert(read["i6"].to!short==-32_767);
    assert(read["i7"].to!short==-32_768);
    assert(read["i8"].to!int==-65_535);
    assert(read["i9"].to!int==-65_536);
    assert(read["iA"].to!int==-2_147_483_647);
    assert(read["iB"].to!int==-2_147_483_648);
    hts_log_info(__FUNCTION__,"Testing float");
    assert(read["F0"].to!float==-1.0);
    assert(read["F1"].to!float==0.0);
    assert(read["F2"].to!float==1.0);
    hts_log_info(__FUNCTION__,"Running tag checking");
    assert(read["I0"].check!ubyte==true);
    assert(read["I5"].check!ushort==true);
    assert(read["I9"].check!uint==true);
    assert(read["i1"].check!byte==true);
    assert(read["i4"].check!short==true);
    assert(read["i8"].check!int==true);
    assert(read["F0"].check!float==true);
    readrange.popFront;
    read=readrange.front;
    hts_log_info(__FUNCTION__,"Testing arrays");
    assert(read["Bs"].to!(short[])==[-32_768,-32_767,0,32_767]);
    assert(read["Bi"].to!(int[])==[-2_147_483_648,-2_147_483_647,0,2_147_483_647]);
    assert(read["BS"].to!(ushort[])==[0,32_767,32_768,65_535]);
    assert(read["BI"].to!(uint[])==[0,2_147_483_647,2_147_483_648,4_294_967_295]);
    hts_log_info(__FUNCTION__,"Running tag checking");
    assert(read["Bs"].check!(short[])==true);
    assert(read["Bi"].check!(int[])==true);
    assert(read["BS"].check!(ushort[])==true);
    assert(read["BI"].check!(uint[])==true);
}
