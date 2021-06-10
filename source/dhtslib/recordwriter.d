module dhtslib.recordwriter;

import std.stdio : File, stdout;
import std.range : ElementType;

/**
* Abstraction for writing string-based records to a file.
* Intended for use with BedRecord, GFF(2|3)Record.
*/
struct RecordWriter(RecType)
if(__traits(hasMember, RecType, "toString"))
{

    File f; /// File writer

    /// ctor with filename and header
    this(string fn, string header="")
    {
        /// open file and write header if any
        if(fn == "-"){
            this.f = stdout;
        }else{
            this.f = File(fn, "w");
        }
        if(header != "") this.f.writeln(header);
    }

    /// Write record to file
    void write(RecType rec){
        this.f.writeln(rec.toString);
    }

    /// Write a range of records to file
    void writeRecords(Range)(Range range)
    if(is(ElementType!Range == RecType))
    {
        foreach (rec; range)
        {
            this.write(rec);
        }
    }
}