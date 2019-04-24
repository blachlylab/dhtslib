module dhtslib.fastq;
import dhtslib.bgzf;

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
    FastqRecord rec;

    this(string fn)
    {
        f=BGZFile(fn);
        popFront;
    }

    FastqRecord front()
    {
        return rec;
    }

    void popFront()
    {
        //get id line
        rec.id = f.front;
        f.popFront;

        //get seq line
        rec.sequence = f.front;
        f.popFront;

        //get extra line
        rec.extra = f.front;
        f.popFront;

        //get qscore line
        rec.qscores = f.front;
        f.popFront;
    }

    bool empty()
    {
        return f.empty;
    }

}