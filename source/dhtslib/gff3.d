module dhtslib.gff3;

import std.range : inputRangeObject, InputRangeObject;

import dhtslib.bgzf;
import dhtslib.gff3d;

struct GFFReader(RecType)
if(is(RecType == GTF_Record) || is(RecType == GFF3_Record))
{
    BGZFile file;
    InputRangeObject range;

    this(string fn)
    {
        this.file = BGZFile(fn);
        this.range = this.file.byLineCopy.inputRangeObject;
    }

    RecType front()
    {
        return RecType(this.range.front);
    }

    void popFront()
    {
        this.range.popFront;
    }

    auto empty()
    {
        return this.range.empty;
    }
}

alias GFF3Reader = GFFReader!GFF3_Record;
alias GTFReader = GFFReader!GTF_Record;