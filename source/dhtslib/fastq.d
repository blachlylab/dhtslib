module dhtslib.fastq;
import dhtslib.bgzf;
import std.range.interfaces : InputRange, inputRangeObject; 

/*
 *  Represents a complete fastq read record.
 */ 
struct FastqRecord
{
    string id;
    string sequence;
    string extra;
    string qscores;
}

/*
 *  Parses a bgzipped, gzipped, or plain-text fastq file 
 *  into a range of FastqRecords.
 */ 
struct FastqFile
{
    BGZFile f;
    InputRange!string lines;
    FastqRecord rec;
    bool last;

    this(string fn)
    {
        f = BGZFile(fn);
        lines = f.byLineCopy.inputRangeObject;
        popFront;
    }

    /// Explicit postblit to avoid 
    /// https://github.com/blachlylab/dhtslib/issues/122
    this(this)
    {
        this.f = f;
        this.lines = lines;
        this.rec = rec;
        this.last = last;
    }

    FastqRecord front()
    {
        return rec;
    }

    void popFront()
    {
        //get id line
        rec.id = lines.front;
        lines.popFront;

        //get seq line
        rec.sequence = lines.front;
        lines.popFront;

        //get extra line
        rec.extra = lines.front;
        lines.popFront;

        //get qscore line
        rec.qscores = lines.front;
        if(!lines.empty)
            lines.popFront;
        else
            last = true;
    }

    bool empty()
    {
        return lines.empty && last;
    }

}

///
debug(dhtslib_unittest) unittest
{
    import std.stdio;
    import htslib.hts_log;
    import std.algorithm : map;
    import std.array : array;
    import std.path : buildPath,dirName;
    hts_set_log_level(htsLogLevel.HTS_LOG_INFO);
    hts_log_info(__FUNCTION__, "Testing FastqFile");
    hts_log_info(__FUNCTION__, "Loading test file");

    auto fqs = FastqFile(buildPath(dirName(dirName(dirName(__FILE__))),"htslib","test","fastqs.fq"));
    assert(fqs.array.length == 125);
    // assert(bg.array == ["122333444455555"]);
}