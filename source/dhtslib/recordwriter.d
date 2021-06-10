module dhtslib.recordwriter;

import std.stdio : File;
import std.range : ElementType;


struct RecordWriter(RecType)
if(__traits(hasMember, RecType, "toString"))
{
    File f;

    this(string fn, string header="")
    {
        this.f = File(fn, "w");
        if(header != "") this.f.writeln(header);
    }

    void write(RecType rec){
        this.f.writeln(rec.toString);
    }

    void writeRecords(Range)(Range range)
    if(is(ElementType!Range == RecType))
    {
        foreach (rec; range)
        {
            this.write(rec);
        }
    }
}